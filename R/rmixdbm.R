## -----------------------------------------------------------------------------
## [IT] rmixdbm (esportata): simula n ranking da una mistura di modelli di
## Mallows (con eventuale componente uniforme se p ha G + 1 elementi). La
## componente di ogni ranking e' restituita nell'attributo "component".
## -----------------------------------------------------------------------------
#' Simulate from a mixture of distance-based models
#'
#' @param n Number of rankings.
#' @param centers A G by M matrix of central rankings (or a vector if G = 1).
#' @param lambda Vector of G precisions.
#' @param p Vector of G mixing proportions, or G + 1 with the last one for a
#'   uniform noise component.
#' @param dist The distance (see \code{\link{rankdist}}).
#'
#' @return A n by M matrix of rankings with attribute \code{"component"}.
#'
#' @examples
#' set.seed(4)
#' X <- rmixdbm(500, centers = rbind(1:5, 5:1), lambda = c(1, 1),
#'              p = c(0.5, 0.5), dist = "kendall")
#' table(attr(X, "component"))
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mixdbm}}, \code{\link{rmallows}}
#' @export
rmixdbm <- function(n, centers, lambda, p, dist = "kendall") {
  dist <- .matchdist(dist)
  centers <- .checkrank(centers, "centers")
  G <- nrow(centers); M <- ncol(centers)
  if (length(lambda) != G) stop("'lambda' must have one value per center")
  if (!length(p) %in% c(G, G + 1)) stop("'p' must have G or G + 1 elements")
  comp <- sample.int(length(p), n, replace = TRUE, prob = p)
  X <- matrix(0L, n, M)
  for (g in seq_along(p)) {
    id <- which(comp == g)
    if (!length(id)) next
    if (g > G) {
      X[id, ] <- t(replicate(length(id), sample.int(M)))
      next
    }
    X[id, order(centers[g, ])] <- .rdbmidentity(length(id), M, lambda[g], dist)
  }
  colnames(X) <- colnames(centers)
  attr(X, "component") <- comp
  X
}
