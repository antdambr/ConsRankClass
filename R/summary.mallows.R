## [IT] summary.mallows: come print, piu' il confronto tra distanza media
## osservata e attesa sotto l'uniforme e un intervallo di confidenza per lambda.
#' Summary method for mallows objects
#'
#' As \code{print}, plus the expected distance under the uniform model, an
#' approximate 95\% confidence interval for the precision (conditional on
#' the estimated center) and the BIC.
#'
#' @param object An object of class \code{"mallows"}.
#' @param digits Number of significant digits.
#' @param \dots Not used.
#' @return The object, invisibly.
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mallows}}
#' @export
summary.mallows <- function(object, digits = 4, ...) {
  print(object, digits = digits)
  e0 <- .expectd(0, object$M, object$dist)
  cat("Expected distance under the uniform model: ", format(e0, digits = digits),
      "\n", sep = "")
  ci <- object$lambda + c(-1, 1) * stats::qnorm(0.975) * object$se
  cat("Approximate 95% CI for lambda (center fixed): [",
      paste(format(pmax(ci, 0), digits = digits), collapse = ", "), "]\n", sep = "")
  cat("BIC = ", format(object$BIC, digits = digits + 2),
      "  (df = ", object$df, ", center not counted)\n\n", sep = "")
  invisible(object)
}
