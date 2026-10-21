## -----------------------------------------------------------------------------
## [IT] mallowsnorm, mallowsexpect, mallowslambda (esportate): versioni per
## l'utente di .logpsi, .expectd/.vard e .solvelambda, vettorizzate su lambda
## (o su dbar). Corrispondono a funzioni di PerMallows (expectation.mm e la
## stima di theta dentro lmm), estese a tutte le distanze.
## -----------------------------------------------------------------------------
#' Normalising constant, moments and precision of the Mallows model
#'
#' For the Mallows (distance-based) model
#' \eqn{P(r) = C(\lambda) \exp\{-\lambda d(r, R)\}}{P(r) = C(lambda) exp(-lambda d(r, R))},
#' \code{mallowsnorm} computes \eqn{\psi(\lambda) = 1/C(\lambda)}{psi(lambda) = 1/C(lambda)},
#' \code{mallowsexpect} the expected distance from the center (and its
#' variance) and \code{mallowslambda} the maximum likelihood estimate of
#' \eqn{\lambda}{lambda} given the mean distance of the data from the center.
#'
#' @param lambda Precision (dispersion) parameter(s), \eqn{\lambda \ge 0}{lambda >= 0}.
#' @param M Number of objects.
#' @param dist The distance (see \code{\link{rankdist}}).
#' @param log If \code{TRUE} (default) \eqn{\log\psi(\lambda)}{log psi(lambda)}.
#' @param var If \code{TRUE} the variance of the distance is returned too.
#' @param meandist Mean (possibly weighted) distance of the rankings from the
#'   center.
#' @param lambdamax Upper bound for \eqn{\lambda}{lambda}.
#'
#' @details By right-invariance the constant does not depend on the center.
#' For Kendall and Cayley closed forms are used (any M); for the other
#' distances the exact distribution of \code{\link{distfreq}}. The MLE of
#' \eqn{\lambda}{lambda} solves \eqn{E_\lambda[d] = \bar d}{E_lambda[d] = meandist}.
#'
#' @return \code{mallowsnorm}: \eqn{\log\psi}{log psi} or \eqn{\psi}{psi}.
#'   \code{mallowsexpect}: the expectation (a matrix with expectation and
#'   variance if \code{var = TRUE}). \code{mallowslambda}: the estimate(s).
#'
#' @examples
#' mallowsnorm(0.5, M = 10, dist = "kendall")
#' mallowsexpect(c(0, 0.5, 1), M = 10, dist = "cayley", var = TRUE)
#' mallowslambda(10, M = 10, dist = "kendall")
#'
#' @references Mallows, C.L. (1957). Non-null ranking models. I.
#'   Biometrika, 44, 114-130.
#' @references Fligner, M.A. and Verducci, J.S. (1986). Distance based ranking
#'   models. Journal of the Royal Statistical Society B, 48, 359-369.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mallows}}, \code{\link{distfreq}}
#' @export
mallowsnorm <- function(lambda, M, dist = "kendall", log = TRUE) {
  dist <- .matchdist(dist)
  lp <- vapply(lambda, .logpsi, 0, M = M, dist = dist)
  if (log) lp else exp(lp)
}

#' @rdname mallowsnorm
#' @export
mallowsexpect <- function(lambda, M, dist = "kendall", var = FALSE) {
  dist <- .matchdist(dist)
  e <- vapply(lambda, .expectd, 0, M = M, dist = dist)
  if (!var) return(e)
  cbind(lambda = lambda, expectation = e,
        variance = vapply(lambda, .vard, 0, M = M, dist = dist))
}

#' @rdname mallowsnorm
#' @export
mallowslambda <- function(meandist, M, dist = "kendall", lambdamax = 50) {
  dist <- .matchdist(dist)
  vapply(meandist, .solvelambda, 0, M = M, dist = dist, lambdamax = lambdamax)
}
