## -----------------------------------------------------------------------------
## [IT] mixdbm (esportata): stima di una mistura di G modelli di Mallows con
## la distanza scelta (Murphy & Martin, 2003). Interfaccia come cca():
## mixdbm(X, G, ..., control = mixdbmcontrol(...)) e le opzioni di controllo
## si possono anche passare direttamente (ad es. nstart = 10).
## -----------------------------------------------------------------------------
#' Mixtures of distance-based (Mallows) models for rankings
#'
#' Fits the mixture of distance-based models of Murphy and Martin (2003),
#' \deqn{f(r) = \sum_{g=1}^G p_g C(\lambda_g) \exp\{-\lambda_g d(r, R_g)\} \; (+ \, p_0 / M!)}{f(r) = sum_g p_g C(lambda_g) exp(-lambda_g d(r, R_g)) (+ p_0 / M!)}
#' by the EM algorithm from multiple random starts.
#'
#' @param X A n by M matrix of complete rankings (rank of object j in column j).
#' @param G Number of (non-uniform) components. \code{G = 0} with
#'   \code{noise = TRUE} gives the uniform model ("0+N" in the paper).
#' @param dist The distance: \code{"kendall"}, \code{"cayley"},
#'   \code{"hamming"}, \code{"ulam"}, \code{"footrule"}, \code{"spearman"}
#'   (as in the paper) or \code{"spearmansq"} (see \code{\link{rankdist}}).
#' @param lambda Constraint on the precisions: \code{"free"} or \code{"equal"}.
#' @param noise If \code{TRUE} a uniform (\eqn{\lambda = 0}{lambda = 0})
#'   noise component is added.
#' @param startfrom Optional \code{mixdbm} fit with G - 1 components (same
#'   distance and noise option): three extra starts are made of its solution
#'   plus a new component.
#' @param control A list of options, see \code{\link{mixdbmcontrol}}.
#' @param \dots Options passed to \code{\link{mixdbmcontrol}} bypassing
#'   \code{control}.
#'
#' @details Data are collapsed into distinct rankings with frequencies. The
#' M-step for the centers is exact for Kendall (exhaustive search, or the
#' branch-and-bound of \code{\link[ConsRank]{consrank}} when M > 8), footrule,
#' Hamming and squared Spearman, and for the other distances when M <= 7;
#' otherwise a local search is used (generalised EM). The normalising
#' constants are exact. BIC = 2 loglik - df log(n) and ICL = BIC + 2 sum z log z
#' as in the paper (larger is better). Components are ordered by decreasing
#' proportion; the noise component is the last one.
#'
#' A component is degenerate when its precision reaches \code{lambdamax}
#' and its expected size is below \code{minsize} (see
#' \code{\link{mixdbmcontrol}}): it describes a few copies of one ranking,
#' not a group of similar rankings, and the likelihood increases without
#' limit towards lambda = Inf. Such solutions can win by BIC and ICL, so the
#' best start without degenerate components is kept; if every start is
#' degenerate the best one is returned with a warning.
#'
#' @return An object of class \code{"mixdbm"}: \code{centers} (G by M matrix
#' of central rankings), \code{ocenters} (centers as orderings of the
#' objects), \code{lambda} (precisions), \code{p} (mixing proportions, noise
#' last), \code{z} (n by (G + noise) posterior membership probabilities),
#' \code{cluster} (crisp partition: column of \code{z} with the highest
#' probability), \code{loglik}, \code{df}, \code{BIC}, \code{ICL}, \code{n},
#' \code{M}, \code{G}, \code{noise}, \code{dist}, \code{lambdamodel},
#' \code{model}, \code{iter}, \code{converged}, \code{lltrace} (log-likelihood
#' at each EM iteration of the best start), \code{loglikstarts} (final
#' log-likelihood of every start), \code{degstarts} (logical: was each start
#' degenerate?), \code{degenerate} (names of the degenerate components of the
#' returned solution: empty unless every start was degenerate),
#' \code{control}.
#'
#' @examples
#' \donttest{
#' data(Irish)
#' ## complete rankings only (rankings with ties are not supported yet)
#' X <- Irish$rankings[apply(Irish$rankings, 1, function(r) all(sort(r) == 1:5)), ]
#' set.seed(1)
#' fit <- mixdbm(X, G = 2, dist = "kendall", nstart = 5)
#' fit
#' summary(fit)
#' }
#'
#' @references Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based
#'   models for ranking data. Computational Statistics & Data Analysis, 41,
#'   645-655.
#' @references Celeux, G. and Diebolt, J. (1985). The SEM algorithm: a
#'   probabilistic teacher algorithm derived from the EM algorithm for the
#'   mixture problem. Computational Statistics Quarterly, 2, 73-82.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mixdbmcontrol}}, \code{\link{mixdbmselect}},
#'   \code{\link{rmixdbm}}, \code{\link{mallows}}, \code{\link{cca}}
#' @export
mixdbm <- function(X, G, dist = "kendall", lambda = c("free", "equal"),
                   noise = FALSE, startfrom = NULL,
                   control = mixdbmcontrol(...), ...) {
  dist <- .matchdist(dist)
  lambdamodel <- match.arg(lambda)
  if (G < 0 || G != round(G)) stop("'G' must be a non-negative integer")
  if (G == 0 && !noise) stop("G = 0 requires noise = TRUE")
  X <- .checkrank(X, "X")
  ag <- .aggregaterank(X)
  ctx <- .makectx(ag$U, ncol(X), dist, control$centermethod, control$lambdamax)
  .logpsi(0, ncol(X), dist)              # errore subito se M e' troppo grande
  fit <- .mixdbmfit(ag, nrow(X), .objnames(X), G, lambdamodel, noise, ctx,
                    control, startfrom)
  fit$call <- match.call()
  if (length(fit$degenerate))
    warning(.degeneratemsg(fit), call. = FALSE)
  fit
}
