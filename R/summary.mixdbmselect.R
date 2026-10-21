## [IT] summary.mixdbmselect: come print, piu' la stampa completa dei due
## modelli scelti.
#' Summary method for mixdbmselect objects
#'
#' As \code{print}, plus the full printout of the models chosen by BIC and
#' by ICL.
#'
#' @param object An object of class \code{"mixdbmselect"}.
#' @param top Number of models shown for each criterion.
#' @param \dots Not used.
#' @return The object, invisibly.
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mixdbmselect}}
#' @export
summary.mixdbmselect <- function(object, top = 10, ...) {
  print(object, top = top)
  cat("Best model by BIC:\n"); print(object$bestBIC)
  cat("Best model by ICL:\n"); print(object$bestICL)
  invisible(object)
}
