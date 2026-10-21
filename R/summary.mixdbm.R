## [IT] summary.mixdbm: come print, piu' numerosita' delle classi (partizione
## netta), probabilita' a posteriori media per classe, convergenza e quante
## partenze hanno raggiunto la log-verosimiglianza migliore (indicazione
## sulla presenza di massimi locali).
#' Summary method for mixdbm objects
#'
#' As \code{print}, plus the size of each cluster of the crisp partition,
#' the mean posterior probability within each cluster, the convergence of
#' the EM algorithm, how many random starts reached the best
#' log-likelihood (an indication of local maxima) and how many degenerate
#' starts were discarded.
#'
#' @param object An object of class \code{"mixdbm"}.
#' @param digits Number of significant digits.
#' @param \dots Not used.
#' @return The object, invisibly.
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mixdbm}}
#' @export
summary.mixdbm <- function(object, digits = 3, ...) {
  print(object, digits = digits)
  cn <- colnames(object$z)
  sz <- table(factor(cn[object$cluster], levels = cn))
  pm <- vapply(seq_along(cn), function(g)
    if (any(object$cluster == g)) mean(object$z[object$cluster == g, g]) else NA, 0)
  cat("Crisp partition (highest posterior probability):\n")
  tab <- rbind(size = as.vector(sz), mean.post = round(pm, digits))
  colnames(tab) <- cn
  print(tab)
  ll <- object$loglikstarts
  if (!is.null(object$degstarts) && !all(object$degstarts)) ll <- ll[!object$degstarts]
  nbest <- sum(ll > max(ll) - 1e-4)
  cat("\nEM iterations: ", object$iter, if (object$converged) " (converged)" else
        " (NOT converged)", ".  Starts reaching the best loglik: ", nbest, "/",
      length(object$loglikstarts), ".\n", sep = "")
  ndeg <- sum(object$degstarts)
  if (ndeg && !length(object$degenerate))
    cat("Degenerate starts discarded: ", ndeg, " (best discarded loglik = ",
        format(max(object$loglikstarts[object$degstarts]), nsmall = 3), ").\n", sep = "")
  cat("\n")
  invisible(object)
}
