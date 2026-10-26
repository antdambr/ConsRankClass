#' Projection of rankings on the principal planes of the permutation polytope
#'
#' Projects a set of rankings onto the directions along which they are farthest
#' from the center of the permutation polytope (Marden, 1995, Section 2.2): the
#' biplot of Cohen and Mallows (1980), similar to principal components but
#' centred at the center of the polytope instead of at the mean ranking. The
#' poles of the objects are projected too. Use \code{plot} to display the
#' projection.
#'
#' @param X A n by m data matrix (or data frame), in which there are n judges
#'   and m objects to be judged. Each row is a ranking of the objects, which
#'   are represented by the columns; ties are allowed, and so are missing
#'   values (\code{NA} or \code{NaN}) for objects not ranked by a judge (e.g.
#'   top-k partial rankings), which are treated as tied in the last position.
#' @param wk Optional vector of n non-negative weights (e.g. the frequencies of
#'   the distinct rankings in the rows of \code{X}). Default \code{NULL}: every
#'   judge has weight 1.
#' @param k The number of dimensions returned. Default 2; at most m - 1.
#'
#' @details Let y be the vector of the ranks given by a judge (with ties, the
#'   midranks of the buckets, so that the ranks sum to m(m + 1)/2) and c the
#'   center of the polytope, whose elements are all (m + 1)/2. The axes are the
#'   first k eigenvectors of \eqn{W = \sum w (y - c)(y - c)' / \sum w}{W =
#'   sum w (y - c)(y - c)' / sum w}, and their eigenvalues are the mean squared
#'   distances of the rankings from the center along each axis. The scores are
#'   the projections of y - c divided by the distance of a permutation from the
#'   center, \eqn{\sqrt{m(m^2-1)/12}}{sqrt(m(m^2 - 1)/12)}, so that they lie in
#'   the unit ball: permutations on its surface, rankings with ties inside.
#'
#'   The pole of an object (Marden, 1995, Definition 2.2) is the half line from
#'   the center through the average of the permutations that rank the object
#'   first (1 for the object, 1 + m/2 for the others): the closer a ranking to
#'   the pole, the better the rank of the object. The poles are returned as unit
#'   vectors projected on the axes, so their length (at most 1) tells how well
#'   the object is represented in the projection. The sign of each axis is
#'   chosen so that its largest loading is positive.
#'
#' @return An object of class \code{rankproj}, a list with:
#' \tabular{ll}{
#' scores \tab the scores of the distinct rankings (one row each)\cr
#' rankings \tab the distinct rankings, in standard form (1, 2, ... for the buckets)\cr
#' weights \tab the number (total weight) of judges with each distinct ranking\cr
#' idx \tab the distinct ranking of each row of \code{X} (\code{NA} for rows with zero weight)\cr
#' wk \tab the weights of the rows of \code{X} (1 if not given)\cr
#' poles \tab the projected poles of the objects (one row each)\cr
#' loadings \tab the m by k matrix of the axes (eigenvectors)\cr
#' eigenvalues \tab all the non-zero eigenvalues of W\cr
#' prop \tab the proportion of the sum of the eigenvalues of each axis\cr
#' nj \tab the number of judges (the sum of the weights)\cr
#' items \tab the names of the objects
#' }
#'
#' @references Marden, J.I. (1995). Analyzing and Modeling Rank Data. Chapman &
#'   Hall, London.
#'
#'   Cohen, A., and Mallows, C.L. (1980). Analysis of ranking data. Technical
#'   memorandum, Bell Telephone Laboratories, Murray Hill, NJ.
#'
#' @author Luca Coraggio \email{luca.coraggio@unina.it}
#'
#' @seealso \code{\link{plot.rankproj}}, \code{\link{marginals}},
#'   \code{\link{pairmatrix}}
#'
#' @examples
#' data(sports, package = "ConsRank")
#' pr <- rankproj(sports)
#' pr
#' plot(pr)
#' data(Univranks)
#' plot(rankproj(Univranks$rankings), groups = Univranks$predictors$Stud)
#'
#' @export

rankproj <- function(X, wk = NULL, k = 2) {

  R <- rankbuckets(X, wk, "rankproj")                  # distinct rankings, with weights
  m <- length(R$items)                                 # number of objects
  if (length(k) != 1 || !is.numeric(k) || k != round(k) || k < 1 || k > m - 1)
    stop("k must be a whole number between 1 and the number of objects minus one (", m - 1, ").")
  Ymid <- (R$lo + R$hi) / 2                            # midranks: ties get the mean of their positions
  cen <- (m + 1) / 2                                   # center of the polytope (all elements)
  Z <- Ymid - cen                                      # rankings centred at the center of the polytope
  nj <- sum(R$w)                                       # number of judges (total weight)
  W <- crossprod(Z * sqrt(R$w)) / nj                   # W = sum w (y - c)(y - c)' / sum w
  ee <- eigen(W, symmetric = TRUE)                     # eigenvalues and eigenvectors
  G <- ee$vectors[, seq_len(k), drop = FALSE]          # the first k axes
  for (j in seq_len(k)) {                              # sign of each axis
    if (G[which.max(abs(G[, j])), j] < 0) G[, j] <- -G[, j] # largest loading positive
  }
  s <- sqrt(m * (m^2 - 1) / 12)                        # distance of a permutation from the center
  scores <- (Z %*% G) / s                              # scores in the unit ball

  zp <- matrix(1 + m / 2 - cen, m, m)                  # poles: 1 + m/2 for the other objects ...
  diag(zp) <- 1 - cen                                  # ... and 1 for the object (centred)
  zp <- zp / sqrt(rowSums(zp^2))                       # unit vectors (one row per object)
  poles <- zp %*% G                                    # projected poles

  ev <- ee$values[ee$values > 1e-12 * max(ee$values)]  # non-zero eigenvalues (at most m - 1)
  axn <- paste0("Dim", seq_len(k))                     # names of the axes
  dimnames(scores) <- list(NULL, axn)                  # scores: one row per distinct ranking
  dimnames(poles) <- list(R$items, axn)                # poles: one row per object
  dimnames(G) <- list(R$items, axn)                    # loadings
  out <- list(scores = scores, rankings = R$Y, weights = R$w, idx = R$idx, wk = R$wk,
              poles = poles, loadings = G, eigenvalues = ev,
              prop = ev[seq_len(k)] / sum(ev), nj = nj, items = R$items) # the results
  class(out) <- c("rankproj", "ConsRankClass")         # class for print and plot
  return(out)
}

#-------------------------------------------------------------------------------
#' Print a projection of rankings
#'
#' Print method for objects of class \code{rankproj}.
#'
#' @param x An object of class \code{rankproj}.
#' @param digits Number of decimal digits. Default 3.
#' @param \dots Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @export

print.rankproj <- function(x, digits = 3, ...) {

  k <- ncol(x$scores)                                  # number of dimensions
  cat("\nProjection of", format(x$nj), "rankings (", nrow(x$rankings), "distinct ) of",
      length(x$items), "objects on", k, "axes of the permutation polytope\n") # header
  cat("\nEigenvalues and proportions:\n")
  ev <- rbind(eigenvalue = x$eigenvalues[seq_len(k)], proportion = x$prop,
              cumulative = cumsum(x$prop))             # how much each axis explains
  colnames(ev) <- colnames(x$scores)                   # Dim1, Dim2, ...
  print(round(ev, digits))
  cat("\nProjected poles of the objects (length 1 = fully represented):\n")
  print(round(cbind(x$poles, length = sqrt(rowSums(x$poles^2))), digits)) # poles and their lengths
  invisible(x)
}
