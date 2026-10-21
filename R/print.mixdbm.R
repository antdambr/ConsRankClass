## [IT] print.mixdbm: stampa nello stile di print.cca (centri come
## ordinazioni, precisioni, proporzioni, criteri di informazione).
#' Print method for mixdbm objects
#'
#' Prints the cluster centers (as orderings and as rankings), precisions,
#' mixing proportions, log-likelihood, BIC and ICL.
#'
#' @param x An object of class \code{"mixdbm"}.
#' @param digits Number of significant digits.
#' @param \dots Not used.
#' @return The object, invisibly.
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mixdbm}}
#' @export
print.mixdbm <- function(x, digits = 3, ...) {
  cat("\nMixture of distance-based models (Murphy & Martin, 2003)\n")
  cat("Model: ", x$model, ".  Rankings: ", x$n, ".  Objects: ", x$M, ".\n\n", sep = "")
  if (x$G > 0) {
    cat("Cluster centers (orderings):\n")
    for (g in seq_len(x$G))
      cat(sprintf("  %-3s %s\n", rownames(x$centers)[g],
                  paste(x$ocenters[g, ], collapse = " > ")))
    cat("\n")
    tab <- data.frame(x$centers, check.names = FALSE)
    tab$lambda <- signif(x$lambda, digits)
    tab$prop <- signif(x$p[seq_len(x$G)], digits)
    print(tab)
  }
  if (x$noise) cat("Noise (uniform) component proportion: ",
                   format(x$p[length(x$p)], digits = digits), "\n", sep = "")
  cat(sprintf("\nloglik = %.3f   df = %d   BIC = %.3f   ICL = %.3f\n\n",
              x$loglik, as.integer(x$df), x$BIC, x$ICL))
  if (length(x$degenerate))
    cat("Warning: degenerate component(s) ", paste(x$degenerate, collapse = ", "),
        " (precision at lambdamax, few rankings): every start was degenerate.\n\n",
        sep = "")
  invisible(x)
}
