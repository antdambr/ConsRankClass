## [IT] logLik.mallows: log-verosimiglianza come oggetto "logLik" (cosi'
## funzionano AIC() e BIC() di stats).
#' Log-likelihood of a mallows object
#'
#' @param object An object of class \code{"mallows"}.
#' @param \dots Not used.
#' @return An object of class \code{"logLik"} (one degree of freedom: the
#'   center is not counted), usable with \code{AIC} and \code{BIC}.
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mallows}}
#' @export
logLik.mallows <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = object$n, class = "logLik")
