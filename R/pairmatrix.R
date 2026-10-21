#' Matrix of pairwise preferences
#'
#' Computes the pairs matrix of a set of rankings (Marden, 1995, Section 2.4):
#' for each ordered pair of objects (a, b), the number of judges who prefer a to
#' b, and the number of judges who tie them. Use \code{plot} to display it as a
#' heatmap.
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
#' @details A smaller rank means a preferred object. The counts are kept
#'   separate: \code{wins[a, b]} judges prefer a to b, \code{wins[b, a]} prefer b
#'   to a, and \code{ties[a, b]} tie them, so \code{wins[a, b] + wins[b, a] +
#'   ties[a, b]} is the number (total weight) of the judges. An object ranked by
#'   a judge is preferred to the objects not ranked by that judge, and two
#'   unranked objects are tied.
#'
#'   Three matrices of proportions are returned. In \code{prop} each tie counts
#'   as half a preference for each object, so \code{prop[a, b] + prop[b, a] = 1}
#'   and, for complete rankings without ties, \code{prop} is the pairs matrix
#'   of Marden (1995); \code{strict} contains the proportions of strict
#'   preferences (ties excluded from the numerator) and \code{tieprop} the
#'   proportions of ties. The score of each object is the mean of its row of
#'   \code{prop}, i.e. the normalised Borda score: the proportion of the other
#'   objects it beats, ties counting one half. The combined input matrix of
#'   Emond and Mason (2002), used by \code{consrank}, is \code{wins + ties -
#'   t(wins)} (with zero diagonal).
#'
#' @return An object of class \code{pairmatrix}, a list with:
#' \tabular{ll}{
#' wins \tab the m by m matrix of the number of judges preferring the row object to the column object\cr
#' ties \tab the m by m symmetric matrix of the number of judges tying the two objects\cr
#' prop \tab the proportions of preferences, ties counting one half (\code{NA} on the diagonal)\cr
#' strict \tab the proportions of strict preferences (\code{NA} on the diagonal)\cr
#' tieprop \tab the proportions of ties (\code{NA} on the diagonal)\cr
#' score \tab the normalised Borda score of each object (mean of its row of \code{prop})\cr
#' nj \tab the number of judges (the sum of the weights)\cr
#' items \tab the names of the objects
#' }
#'
#' @references Marden, J.I. (1995). Analyzing and Modeling Rank Data. Chapman &
#'   Hall, London.
#'
#'   Emond, E.J., and Mason, D.W. (2002). A new rank correlation coefficient with
#'   application to the consensus ranking problem. Journal of Multi-Criteria
#'   Decision Analysis, 11, 17-28.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#'
#' @seealso \code{\link{plot.pairmatrix}}, \code{\link{marginals}}
#'
#' @examples
#' data(Univranks)
#' pm <- pairmatrix(Univranks$rankings)
#' pm
#' plot(pm)
#' plot(pm, type = "ties")
#'
#' @export

pairmatrix <- function(X, wk = NULL) {

  R <- rankbuckets(X, wk, "pairmatrix")                # distinct rankings in standard form, with weights
  m <- length(R$items)                                 # number of objects
  pr <- which(upper.tri(diag(m)), arr.ind = TRUE)      # pairs (a, b), a < b
  D <- R$Y[, pr[, 2], drop = FALSE] - R$Y[, pr[, 1], drop = FALSE] # > 0: a preferred; < 0: b preferred
  pa <- colSums(R$w * (D > 0))                         # judges preferring a to b
  pb <- colSums(R$w * (D < 0))                         # judges preferring b to a
  pt <- colSums(R$w * (D == 0))                        # judges tying a and b

  wins <- matrix(0, m, m, dimnames = list(R$items, R$items)) # row object preferred to column object
  wins[pr] <- pa                                       # a over b
  wins[pr[, 2:1, drop = FALSE]] <- pb                  # b over a
  ties <- matrix(0, m, m, dimnames = list(R$items, R$items)) # judges tying the two objects
  ties[pr] <- pt                                       # upper triangle
  ties[pr[, 2:1, drop = FALSE]] <- pt                  # symmetric

  nj <- sum(R$w)                                       # number of judges (total weight)
  prop <- (wins + ties / 2) / nj                       # ties count one half
  strict <- wins / nj                                  # strict preferences only
  tieprop <- ties / nj                                 # proportion of ties
  diag(prop) <- diag(strict) <- diag(tieprop) <- NA    # no comparison of an object with itself
  score <- rowSums(prop, na.rm = TRUE) / (m - 1)       # normalised Borda score

  out <- list(wins = wins, ties = ties, prop = prop, strict = strict,
              tieprop = tieprop, score = score, nj = nj, items = R$items) # the results
  class(out) <- c("pairmatrix", "ConsRankClass")       # class for print and plot
  return(out)
}

#-------------------------------------------------------------------------------
#' Print a matrix of pairwise preferences
#'
#' Print method for objects of class \code{pairmatrix}.
#'
#' @param x An object of class \code{pairmatrix}.
#' @param digits Number of decimal digits of the proportions. Default 3.
#' @param \dots Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @export

print.pairmatrix <- function(x, digits = 3, ...) {

  cat("\nPairwise preferences of", format(x$nj), "judges and", length(x$items), "objects\n") # header
  cat("\nProportion of judges preferring the row object to the column object (ties count one half):\n\n")
  print(round(x$prop, digits))                         # proportions with half ties
  if (any(x$ties > 0)) {                               # ties: show their proportions too
    cat("\nProportion of ties:\n\n")
    print(round(x$tieprop, digits))
  }
  cat("\nNormalised Borda scores:\n")
  print(round(sort(x$score, decreasing = TRUE), digits)) # from the best object
  invisible(x)
}
