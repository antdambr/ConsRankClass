## -----------------------------------------------------------------------------
## [IT] mixdbmselect (esportata): stima tutti i modelli di una griglia
## (distanze, vincoli sulle precisioni, componente uniforme si'/no, G) e li
## confronta con BIC e ICL, come nelle Tabelle 2-5 del paper. Con
## nested = TRUE il modello con G componenti parte anche dal migliore con
## G - 1 componenti piu' una nuova componente.
## -----------------------------------------------------------------------------
#' Model selection for mixtures of distance-based models
#'
#' Fits \code{\link{mixdbm}} over a grid of distances, precision
#' constraints, noise options and numbers of components, and compares the
#' models by BIC and ICL.
#'
#' @param X A n by M matrix of complete rankings.
#' @param G Vector of numbers of components (0 is used only with noise).
#' @param dist Vector of distances (see \code{\link{rankdist}}).
#' @param lambda Vector of precision constraints (\code{"free"}, \code{"equal"}).
#' @param noise Vector of noise options (\code{FALSE}, \code{TRUE}).
#' @param nested If \code{TRUE} each fit with G components also starts from
#'   the best fit with G - 1 components plus one new component.
#' @param trace If \code{TRUE} (default) prints one line per fitted model.
#' @param control,\dots Options, see \code{\link{mixdbmcontrol}}.
#'
#' @return An object of class \code{"mixdbmselect"}: \code{table} (one row per
#'   model with loglik, df, BIC, ICL and \code{degenerate}, TRUE when every
#'   start gave a degenerate solution, see \code{\link{mixdbm}}), \code{fits} (all the fitted models,
#'   named as in the table), \code{bestBIC}, \code{bestICL} (the best
#'   models among the non-degenerate ones).
#'
#' @examples
#' \donttest{
#' data(Irish)
#' ## complete rankings only (rankings with ties are not supported yet)
#' X <- Irish$rankings[apply(Irish$rankings, 1, function(r) all(sort(r) == 1:5)), ]
#' set.seed(1)
#' sel <- mixdbmselect(X, G = 1:3, dist = c("kendall", "cayley"),
#'                     lambda = "free", noise = FALSE, nstart = 5)
#' sel
#' }
#'
#' @references Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based
#'   models for ranking data. Computational Statistics & Data Analysis, 41,
#'   645-655.
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mixdbm}}
#' @export
mixdbmselect <- function(X, G = 1:5, dist = c("kendall", "spearman", "cayley"),
                         lambda = c("free", "equal"), noise = c(FALSE, TRUE),
                         nested = TRUE, trace = TRUE,
                         control = mixdbmcontrol(...), ...) {
  G <- sort(G)
  dist <- vapply(dist, .matchdist, "")
  lambda <- match.arg(lambda, several.ok = TRUE)
  X <- .checkrank(X, "X")
  ag <- .aggregaterank(X)
  labs <- .objnames(X)
  fits <- list(); rows <- list()
  for (d in dist) {
    ctx <- .makectx(ag$U, ncol(X), d, control$centermethod, control$lambdamax)
    for (nz in noise) for (lam in lambda) for (g in G) {
      if (g == 0 && !nz) next
      if (g <= 1 && lam == "equal") next           # identico a "free"
      prev <- fits[[.modelname(d, if (g - 1 <= 1) "free" else lam, g - 1, nz)]]
      f <- .mixdbmfit(ag, nrow(X), labs, g, lam, nz, ctx, control,
                      if (nested) prev)
      fits[[f$model]] <- f
      rows[[f$model]] <- data.frame(model = f$model, dist = d,
        lambda = if (g <= 1) "free" else lam, G = g, noise = nz,
        loglik = f$loglik, df = f$df, BIC = f$BIC, ICL = f$ICL,
        degenerate = length(f$degenerate) > 0)
      if (trace)
        cat(sprintf("%-30s loglik = %10.3f  BIC = %10.3f  ICL = %10.3f%s\n",
                    f$model, f$loglik, f$BIC, f$ICL,
                    if (length(f$degenerate)) "  (degenerate)" else ""))
    }
  }
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  if (any(tab$degenerate))
    warning("degenerate solutions (every start degenerate) for: ",
            paste(tab$model[tab$degenerate], collapse = "; "),
            ". They are not eligible as best models.", call. = FALSE)
  ## i modelli degeneri non possono essere scelti (salvo che lo siano tutti)
  ok <- if (all(tab$degenerate)) rep(TRUE, nrow(tab)) else !tab$degenerate
  pick <- function(crit) tab$model[ok][which.max(crit[ok])]
  structure(list(table = tab, fits = fits,
                 bestBIC = fits[[pick(tab$BIC)]],
                 bestICL = fits[[pick(tab$ICL)]]),
            class = "mixdbmselect")
}
