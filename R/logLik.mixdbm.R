## [IT] logLik.mixdbm: log-verosimiglianza come oggetto "logLik".
#' Log-likelihood of a mixdbm object
#'
#' @param object An object of class \code{"mixdbm"}.
#' @param \dots Not used.
#' @return An object of class \code{"logLik"} with the degrees of freedom of
#'   the model (see \code{\link{mixdbm}}), usable with \code{AIC} and \code{BIC}.
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mixdbm}}
#' @export
logLik.mixdbm <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = object$n, class = "logLik")
