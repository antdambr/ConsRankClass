## -----------------------------------------------------------------------------
## [IT] mallows (esportata): stima di massima verosimiglianza del modello di
## Mallows con la distanza scelta (l'equivalente di PerMallows::lmm, con pesi
## e con tutte le distanze). Il centro stimato e' la mediana dei ranking
## (mallowscenter): infatti la log-verosimiglianza profilo e' decrescente
## nella distanza media. Poi lambda risolve E_lambda[d] = distanza media.
## L'errore standard di lambda e' 1 / sqrt(n Var_lambda[d]), condizionato al
## centro stimato. Se il centro e' dato, si stima solo lambda.
## -----------------------------------------------------------------------------
#' Mallows (distance-based) model for complete rankings
#'
#' Maximum likelihood estimation of the Mallows model
#' \eqn{P(r) = C(\lambda) \exp\{-\lambda d(r, R)\}}{P(r) = C(lambda) exp(-lambda d(r, R))}
#' with center (modal ranking) R and precision \eqn{\lambda}{lambda}, for any
#' of the distances of \code{\link{rankdist}}.
#'
#' @param X A n by M matrix of complete rankings (rank of object j in column j).
#' @param dist The distance: \code{"kendall"} (default, the classical Mallows
#'   model), \code{"cayley"}, \code{"hamming"}, \code{"ulam"},
#'   \code{"footrule"}, \code{"spearman"} or \code{"spearmansq"}.
#' @param center Optional known center: if given only \eqn{\lambda}{lambda}
#'   is estimated.
#' @param weights Optional non-negative case weights.
#' @param centermethod Method for the center (see \code{\link{mallowscenter}}).
#' @param lambdamax Upper bound for \eqn{\lambda}{lambda}.
#'
#' @details The MLE of the center is the median ranking with respect to the
#' chosen distance (see \code{\link{mallowscenter}}); given the center, the
#' MLE of \eqn{\lambda}{lambda} solves
#' \eqn{E_\lambda[d] = \bar d}{E_lambda[d] = mean distance}. The standard error
#' of \eqn{\lambda}{lambda} is computed from the Fisher information
#' \eqn{n \, Var_\lambda[d]}{n Var_lambda[d]}, conditionally on the
#' estimated center. The normalising constant is exact (closed forms for
#' Kendall and Cayley; for Ulam \eqn{M \le 50}, for footrule and Spearman
#' \eqn{M \le 14}). As in \code{\link{mixdbm}}, \code{BIC} is
#' 2 loglik - df log(n) (larger is better; \code{stats::BIC(fit)} gives the
#' same value with the opposite sign).
#'
#' @return An object of class \code{"mallows"}: \code{center} (the central
#' ranking), \code{ocenter} (the center as an ordering of the objects),
#' \code{lambda}, \code{se} (its standard error), \code{meandist} (mean
#' distance from the center), \code{loglik}, \code{df}, \code{BIC},
#' \code{n}, \code{M}, \code{dist}.
#'
#' @examples
#' data(Irish)
#' ## complete rankings only (rankings with ties are not supported yet)
#' X <- Irish$rankings[apply(Irish$rankings, 1, function(r) all(sort(r) == 1:5)), ]
#' fit <- mallows(X, dist = "kendall")
#' fit
#' summary(fit)
#' mallows(X, dist = "cayley")
#'
#' @references Mallows, C.L. (1957). Non-null ranking models. I.
#'   Biometrika, 44, 114-130.
#' @references Fligner, M.A. and Verducci, J.S. (1986). Distance based ranking
#'   models. Journal of the Royal Statistical Society B, 48, 359-369.
#' @references Irurozki, E., Calvo, B. and Lozano, J.A. (2016). PerMallows:
#'   An R package for Mallows and generalized Mallows models. Journal of
#'   Statistical Software, 71(12), 1-30.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mixdbm}}, \code{\link{rmallows}},
#'   \code{\link{mallowscenter}}, \code{\link{mallowsnorm}}
#' @export
mallows <- function(X, dist = "kendall", center = NULL, weights = NULL,
                    centermethod = c("auto", "exhaustive", "local", "bb"),
                    lambdamax = 50) {
  cl <- match.call()
  dist <- .matchdist(dist)
  centermethod <- match.arg(centermethod)
  X <- .checkrank(X, "X")
  M <- ncol(X)
  if (!is.null(weights) && (length(weights) != nrow(X) || any(weights < 0)))
    stop("'weights' must be a non-negative vector of length nrow(X)")
  ag <- .aggregaterank(X, weights)
  .logpsi(0, M, dist)                    # errore subito se M e' troppo grande
  if (is.null(center)) {
    ctx <- .makectx(ag$U, M, dist, centermethod, lambdamax)
    R <- .centerestimate(ag$U, ag$w, dist, ctx)
  } else {
    R <- .checkrank(center, "center")[1, ]
    if (length(R) != M) stop("'center' must rank the same objects as 'X'")
  }
  d <- .dist2center(ag$U, R, dist)
  sw <- sum(ag$w)
  dbar <- sum(ag$w * d) / sw
  lam <- .solvelambda(dbar, M, dist, lambdamax)
  loglik <- sum(ag$w * (-lam * d)) - sw * .logpsi(lam, M, dist)
  v <- .vard(lam, M, dist)
  se <- if (v > 0) 1 / sqrt(sw * v) else NA_real_
  names(R) <- .objnames(X)
  df <- 1
  structure(list(center = R, ocenter = .toorderings(R, .objnames(X))[1, ],
                 lambda = lam, se = se, meandist = dbar,
                 loglik = loglik, df = df, BIC = 2 * loglik - df * log(sw),
                 n = sw, M = M, dist = dist, centerknown = !is.null(center),
                 call = cl), class = "mallows")
}
