#' Marginal distribution of the positions of the objects
#'
#' Computes the marginal matrix of a set of rankings (Marden, 1995, Section
#' 2.3): for each object and each position, the number (and the proportion) of
#' judges who put the object in that position, together with the mean rank of
#' each object. Use \code{plot} to display it as a bubble plot.
#'
#' @param X A n by m data matrix (or data frame), in which there are n judges
#'   and m objects to be judged. Each row is a ranking of the objects, which
#'   are represented by the columns; ties are allowed, and so are missing
#'   values (\code{NA} or \code{NaN}) for objects not ranked by a judge (e.g.
#'   top-k partial rankings), which are treated as tied in the last position.
#' @param wk Optional vector of n non-negative weights (e.g. the frequencies of
#'   the distinct rankings in the rows of \code{X}). Default \code{NULL}: every
#'   judge has weight 1.
#'
#' @details The ranks are read as an order: the actual values do not matter, so
#'   \code{c(1, 1, 3)}, \code{c(1, 1, 2)} and \code{c(2, 2, 7)} are the same
#'   ranking. Objects tied in a bucket that occupies the positions p, ..., q
#'   share these positions: the judge contributes 1/(q - p + 1) to each of them
#'   (the average over the random orderings of the tied objects). Hence every
#'   row and every column of the matrix of counts sums to the number (or the
#'   total weight) of the judges, and the mean rank of an object is its mean
#'   rank with ties (the midrank of its bucket). Objects not ranked by a judge
#'   form the last bucket of that judge's ranking. For complete rankings without
#'   ties this is the marginal matrix of Marden (1995).
#'
#' @return An object of class \code{marginals}, a list with:
#' \tabular{ll}{
#' counts \tab the m by m matrix of counts (objects in rows, positions in columns)\cr
#' prop \tab the matrix of proportions (counts divided by \code{nj})\cr
#' meanrank \tab the mean rank of each object\cr
#' nj \tab the number of judges (the sum of the weights)\cr
#' ties \tab \code{TRUE} if some rankings have ties or missing values\cr
#' items \tab the names of the objects
#' }
#'
#' @references Marden, J.I. (1995). Analyzing and Modeling Rank Data. Chapman &
#'   Hall, London.
#'
#' @author Luca Coraggio \email{luca.coraggio@unina.it}
#'
#' @seealso \code{\link{plot.marginals}}, \code{\link{pairmatrix}}
#'
#' @examples
#' data(Univranks)
#' mg <- marginals(Univranks$rankings)
#' mg
#' plot(mg)
#' # top-2 partial rankings (unranked objects coded as NaN)
#' data(EVS)
#' plot(marginals(EVS$rankings))
#'
#' @export

marginals <- function(X, wk = NULL) {

  R <- rankbuckets(X, wk, "marginals")                 # rankings as buckets, with weights
  m <- length(R$items)                                 # number of objects
  counts <- matrix(0, m, m)                            # objects in rows, positions in columns
  for (i in seq_len(nrow(R$lo))) {                     # for each distinct ranking
    lo <- R$lo[i, ]                                    # first position of the bucket of each object
    hi <- R$hi[i, ]                                    # last position of the bucket of each object
    for (j in seq_len(m)) {                            # for each object
      pos <- lo[j]:hi[j]                               # the positions shared by its bucket
      counts[j, pos] <- counts[j, pos] + R$w[i] / length(pos) # equal shares of the judge's weight
    }
  }
  dimnames(counts) <- list(R$items, paste0("r", seq_len(m))) # objects and positions
  nj <- sum(R$w)                                       # number of judges (total weight)
  meanrank <- as.vector(counts %*% seq_len(m)) / nj    # mean position of each object
  names(meanrank) <- R$items                           # named by object

  out <- list(counts = counts, prop = counts / nj, meanrank = meanrank,
              nj = nj, ties = R$ties, items = R$items)  # the results
  class(out) <- c("marginals", "ConsRankClass")        # class for print and plot
  return(out)
}

#-------------------------------------------------------------------------------
## Rankings as buckets, used by marginals and pairmatrix. Checks X and wk,
## codes the missing values as tied in the last position, merges identical
## rankings (summing their weights) and returns, for each distinct ranking and
## each object, the first (lo) and last (hi) position of the bucket of the
## object, the rankings in standard form (Y: the rank of the bucket, 1, 2, ...),
## the weights, the distinct ranking of each row of X (idx, NA for rows with
## zero weight), the weights of the rows (wk), the names of the objects and
## whether there are ties.

rankbuckets <- function(X, wk, fname) {

  if (is.data.frame(X)) X <- as.matrix(X)              # data frames are allowed
  if (is.null(dim(X))) X <- matrix(X, nrow = 1, dimnames = list(NULL, names(X))) # a single ranking
  if (!is.numeric(X))
    stop(fname, "() needs a numeric matrix of rankings: judges in rows, objects in columns.")
  m <- ncol(X)                                         # number of objects
  if (m < 2)
    stop(fname, "() needs at least two objects (columns).")
  items <- colnames(X)                                 # names of the objects
  if (is.null(items)) items <- paste0("o", seq_len(m)) # default names
  if (is.null(wk)) wk <- rep(1, nrow(X))               # default: weight 1 for each judge
  if (!is.numeric(wk) || length(wk) != nrow(X) || anyNA(wk) || any(wk < 0))
    stop("wk must be a vector of non-negative weights, one for each row of X.")
  allna <- rowSums(!is.na(X)) == 0                     # judges who ranked no object
  if (any(allna))
    stop(fname, "(): row(s) ", paste(utils::head(which(allna), 5), collapse = ", "),
         " contain only missing values: remove them.")

  Y <- t(apply(X, 1, function(x) {                     # standard form of each ranking
    x[is.na(x)] <- max(x, na.rm = TRUE) + 1            # unranked objects: tied in the last position
    match(x, sort(unique(x)))                          # rank of the bucket: 1, 2, ...
  }))
  ties <- any(apply(Y, 1, function(y) anyDuplicated(y) > 0)) # some tied objects

  key <- apply(Y, 1, paste, collapse = ",")            # identical rankings have the same key
  first <- !duplicated(key)                            # the first row of each distinct ranking
  w <- as.vector(tapply(wk, factor(key, levels = key[first]), sum)) # total weight of each ranking
  Y <- Y[first, , drop = FALSE]                        # the distinct rankings
  keep <- w > 0                                        # rankings with positive weight
  idx <- match(key, key[first])                        # distinct ranking of each row of X
  idx <- match(idx, which(keep))                       # renumbered after dropping zero weights (NA: dropped)
  Y <- Y[keep, , drop = FALSE]                         # only them
  w <- w[keep]                                         # and their weights
  if (length(w) == 0)
    stop("all the weights wk are zero.")

  lo <- t(apply(Y, 1, function(y) rank(y, ties.method = "min"))) # first position of each bucket
  hi <- t(apply(Y, 1, function(y) rank(y, ties.method = "max"))) # last position of each bucket
  dimnames(Y) <- dimnames(lo) <- dimnames(hi) <- list(NULL, items) # named by object

  return(list(Y = Y, lo = lo, hi = hi, w = w, idx = idx, wk = wk, items = items, ties = ties))
}

#-------------------------------------------------------------------------------
#' Print a marginal matrix
#'
#' Print method for objects of class \code{marginals}.
#'
#' @param x An object of class \code{marginals}.
#' @param digits Number of decimal digits of the proportions. Default 3.
#' @param \dots Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @export

print.marginals <- function(x, digits = 3, ...) {

  cat("\nMarginal matrix of", format(x$nj), "judges and", length(x$items), "objects") # header
  if (x$ties) cat(" (ties or unranked objects shared among their positions)") # how ties were treated
  cat("\n\nProportion of judges putting each object (rows) in each position (columns):\n\n")
  print(round(x$prop, digits))                         # the proportions
  cat("\nMean ranks:\n")
  print(round(sort(x$meanrank), digits))               # mean ranks, from the best object
  invisible(x)
}
