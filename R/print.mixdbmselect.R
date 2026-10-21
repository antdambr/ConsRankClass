## [IT] print.mixdbmselect: i migliori modelli per BIC e per ICL e i due
## modelli scelti.
#' Print method for mixdbmselect objects
#'
#' Prints the best models by BIC and by ICL.
#'
#' @param x An object of class \code{"mixdbmselect"}.
#' @param top Number of models shown for each criterion.
#' @param \dots Not used.
#' @return The object, invisibly.
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mixdbmselect}}
#' @export
print.mixdbmselect <- function(x, top = 10, ...) {
  tab <- x$table
  cat("\nTop models by BIC:\n")
  print(utils::head(tab[order(-tab$BIC), ], top), row.names = FALSE)
  cat("\nTop models by ICL:\n")
  print(utils::head(tab[order(-tab$ICL), ], top), row.names = FALSE)
  cat("\nBest model by BIC:", x$bestBIC$model, "\n")
  cat("Best model by ICL:", x$bestICL$model, "\n\n")
  invisible(x)
}
