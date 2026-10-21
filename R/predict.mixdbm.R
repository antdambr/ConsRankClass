## [IT] predict.mixdbm: probabilita' a posteriori di appartenenza (e
## classe piu' probabile) per nuovi ranking.
#' Cluster membership of new rankings from a mixdbm fit
#'
#' Posterior membership probabilities of new rankings under a fitted
#' mixture of distance-based models.
#'
#' @param object An object of class \code{"mixdbm"}.
#' @param newdata A matrix of complete rankings of the same objects.
#' @param \dots Not used.
#' @return A list with \code{z} (matrix of posterior probabilities, one
#'   column per component, noise last) and \code{cluster} (column of
#'   \code{z} with the highest probability).
#' @examples
#' set.seed(1)
#' X <- rmixdbm(300, centers = rbind(1:5, 5:1), lambda = c(1, 1),
#'              p = c(0.5, 0.5), dist = "kendall")
#' fit <- mixdbm(X, G = 2, dist = "kendall", nstart = 3)
#' predict(fit, rbind(1:5, c(5, 4, 3, 1, 2)))
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mixdbm}}
#' @export
predict.mixdbm <- function(object, newdata, ...) {
  X <- .checkrank(newdata, "newdata")
  if (ncol(X) != object$M) stop("'newdata' must rank ", object$M, " objects")
  lp <- .mixlogdens(object, X)
  z <- exp(lp - .logsumexprows(lp))
  colnames(z) <- names(object$p)
  list(z = z, cluster = max.col(z, ties.method = "first"))
}
