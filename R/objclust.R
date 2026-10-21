#' Hierarchical clustering of the objects
#'
#' Clusters the objects of a set of rankings (Marden, 1995, Example 2.18): the
#' distance between two objects is the mean squared difference between the
#' ranks the judges give them, so objects are close when the judges assign them
#' similar ranks. Use \code{plot} to draw the dendrogram.
#'
#' @param X A n by m data matrix (or data frame), in which there are n judges
#'   and m objects to be judged. Each row is a ranking of the objects, which
#'   are represented by the columns; ties are allowed, and so are missing
#'   values (\code{NA} or \code{NaN}) for objects not ranked by a judge (e.g.
#'   top-k partial rankings), which are treated as tied in the last position.
#' @param wk Optional vector of n non-negative weights (e.g. the frequencies of
#'   the distinct rankings in the rows of \code{X}). Default \code{NULL}: every
#'   judge has weight 1.
#' @param method The agglomeration method of \code{\link[stats]{hclust}}.
#'   Default \code{"complete"} (complete linkage, as in Marden, 1995).
#'
#' @details The distance between objects k and l is
#'   \eqn{\sum w (y_k - y_l)^2 / \sum w}{sum w (y_k - y_l)^2 / sum w}, where
#'   \eqn{y_k} is the rank of object k given by a judge (with ties, the midrank
#'   of its bucket) and w the weight of the judge.
#'
#' @return An object of class \code{objclust}, a list with:
#' \tabular{ll}{
#' dist \tab the distances between the objects (an object of class \code{dist})\cr
#' hclust \tab the hierarchical clustering (an object of class \code{hclust})\cr
#' nj \tab the number of judges (the sum of the weights)\cr
#' items \tab the names of the objects
#' }
#' Use \code{\link[stats]{cutree}} on the component \code{hclust} to get the
#' groups of objects.
#'
#' @references Marden, J.I. (1995). Analyzing and Modeling Rank Data. Chapman &
#'   Hall, London.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#'
#' @seealso \code{\link{plot.objclust}}, \code{\link{rankproj}}
#'
#' @examples
#' data(sports, package = "ConsRank")
#' oc <- objclust(sports)
#' oc
#' plot(oc)
#' stats::cutree(oc$hclust, k = 3)
#'
#' @importFrom stats hclust as.dist
#' @export

objclust <- function(X, wk = NULL, method = "complete") {

  R <- rankbuckets(X, wk, "objclust")                  # distinct rankings, with weights
  m <- length(R$items)                                 # number of objects
  Ymid <- (R$lo + R$hi) / 2                            # midranks: ties get the mean of their positions
  nj <- sum(R$w)                                       # number of judges (total weight)
  S <- crossprod(Ymid * sqrt(R$w)) / nj                # mean products of ranks: sum w y_k y_l / sum w
  D <- outer(diag(S), diag(S), "+") - 2 * S            # mean squared differences: S_kk + S_ll - 2 S_kl
  D[D < 0] <- 0                                        # rounding errors
  dimnames(D) <- list(R$items, R$items)                # named by object
  d <- stats::as.dist(D)                               # distances between the objects
  hc <- stats::hclust(d, method = method)              # hierarchical clustering

  out <- list(dist = d, hclust = hc, nj = nj, items = R$items) # the results
  class(out) <- c("objclust", "ConsRankClass")         # class for print and plot
  return(out)
}

#-------------------------------------------------------------------------------
#' Print a clustering of the objects
#'
#' Print method for objects of class \code{objclust}.
#'
#' @param x An object of class \code{objclust}.
#' @param digits Number of decimal digits. Default 3.
#' @param \dots Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @export

print.objclust <- function(x, digits = 3, ...) {

  cat("\nHierarchical clustering (", x$hclust$method, " linkage) of ", length(x$items),
      " objects ranked by ", format(x$nj), " judges\n", sep = "") # header
  cat("\nMean squared differences of the ranks:\n\n")
  print(round(x$dist, digits))                         # the distances
  cat("\nMerging heights:", round(x$hclust$height, digits), "\n") # heights of the dendrogram
  invisible(x)
}

#-------------------------------------------------------------------------------
#' Plot a clustering of the objects
#'
#' Draws the dendrogram of the hierarchical clustering of the objects.
#'
#' @param x An object of class \code{objclust}, as returned by
#'   \code{\link{objclust}}.
#' @param k Optional number of groups: if given, the groups are enclosed in
#'   rectangles (see \code{\link[stats]{rect.hclust}}).
#' @param main Title of the plot. Default: none.
#' @param ylab Label of the vertical axis. Default "Mean squared rank difference".
#' @param \dots Other arguments passed to \code{plot} for \code{hclust} objects.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{objclust}}
#'
#' @examples
#' data(sports, package = "ConsRank")
#' plot(objclust(sports), k = 3)
#'
#' @importFrom stats rect.hclust
#' @importFrom graphics plot
#' @export

plot.objclust <- function(x, k = NULL, main = "", ylab = "Mean squared rank difference", ...) {

  graphics::plot(x$hclust, hang = -1, main = main, xlab = "", sub = "", ylab = ylab, ...) # dendrogram
  if (!is.null(k)) stats::rect.hclust(x$hclust, k = k, border = "grey40") # groups in rectangles
  invisible(x)
}
