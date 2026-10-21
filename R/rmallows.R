## -----------------------------------------------------------------------------
## [IT] rmallows (esportata): simula n ranking dal modello di Mallows con
## centro R e precisione lambda. Si simula attorno all'identita' e poi si
## rietichettano gli oggetti: se Y e' la permutazione relativa, il ranking e'
## x con x[order(R)] = Y (cosi' d(x, R) = d(Y, e)). Analoga a PerMallows::rmm.
## -----------------------------------------------------------------------------
#' Simulation and density of the Mallows model
#'
#' \code{rmallows} draws rankings from the Mallows (distance-based) model
#' \eqn{P(r) = C(\lambda) \exp\{-\lambda d(r, R)\}}{P(r) = C(lambda) exp(-lambda d(r, R))};
#' \code{dmallows} computes its probability (density) for given rankings.
#'
#' @param n Number of rankings to draw.
#' @param center The central (modal) ranking R (vector; rank of object j in
#'   position j). Its names, if any, become the column names.
#' @param lambda The precision parameter, \eqn{\lambda \ge 0}{lambda >= 0}.
#' @param dist The distance (see \code{\link{rankdist}}).
#' @param X A matrix (or vector) of complete rankings.
#' @param log If \code{TRUE} the log-probabilities.
#'
#' @details Sampling is exact for every distance: inversion vectors (Kendall),
#' the Chinese restaurant process for the Ewens distribution (Cayley), random
#' derangements (Hamming), random standard Young tableaux and inverse
#' Robinson-Schensted correspondence (Ulam, \eqn{M \le 50}), dynamic
#' programming (footrule and Spearman, \eqn{M \le 14}).
#'
#' @return \code{rmallows}: a n by M matrix of rankings. \code{dmallows}: a
#'   vector of probabilities.
#'
#' @examples
#' set.seed(1)
#' X <- rmallows(100, center = 1:6, lambda = 1, dist = "kendall")
#' head(X)
#' dmallows(X[1:5, ], center = 1:6, lambda = 1, dist = "kendall")
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mallows}}, \code{\link{rmixdbm}}
#' @export
rmallows <- function(n, center, lambda, dist = "kendall") {
  dist <- .matchdist(dist)
  R <- .checkrank(center, "center")[1, ]
  if (length(lambda) != 1 || lambda < 0) stop("'lambda' must be a single number >= 0")
  M <- length(R)
  X <- matrix(0L, n, M)
  X[, order(R)] <- .rdbmidentity(n, M, lambda, dist)
  colnames(X) <- names(center)
  X
}

## -----------------------------------------------------------------------------
## [IT] dmallows (esportata): probabilita' (o log-probabilita') dei ranking X
## sotto il modello di Mallows: log P = -lambda d(x, R) - log psi(lambda).
## Analoga a PerMallows::dmm.
## -----------------------------------------------------------------------------
#' @rdname rmallows
#' @export
dmallows <- function(X, center, lambda, dist = "kendall", log = FALSE) {
  dist <- .matchdist(dist)
  X <- .checkrank(X, "X")
  R <- .checkrank(center, "center")[1, ]
  if (length(R) != ncol(X)) stop("'X' and 'center' must rank the same objects")
  lf <- -lambda * .dist2center(X, R, dist) - .logpsi(lambda, ncol(X), dist)
  if (log) lf else exp(lf)
}
