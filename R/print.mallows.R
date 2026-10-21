## [IT] print.mallows: stampa sintetica del modello di Mallows stimato.
#' Print method for mallows objects
#'
#' Prints the estimated center (as ranking and as ordering), the precision
#' with its standard error, the mean distance and the log-likelihood.
#'
#' @param x An object of class \code{"mallows"}.
#' @param digits Number of significant digits.
#' @param \dots Not used.
#' @return The object, invisibly.
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mallows}}
#' @export
print.mallows <- function(x, digits = 4, ...) {
  cat("\nMallows model,", .distlabel(x$dist), "distance\n")
  cat("Objects: ", x$M, ".  Rankings: ", format(x$n), ".\n\n", sep = "")
  cat("Central ranking", if (x$centerknown) "(given)" else "(estimated)", ":\n")
  print(x$center)
  cat("Ordering: ", paste(x$ocenter, collapse = " > "), "\n\n", sep = "")
  cat("lambda = ", format(x$lambda, digits = digits),
      "  (s.e. ", format(x$se, digits = 3), ")\n", sep = "")
  cat("Mean distance from the center: ", format(x$meandist, digits = digits), "\n", sep = "")
  cat("loglik = ", format(x$loglik, digits = digits + 2), "\n\n", sep = "")
  invisible(x)
}
