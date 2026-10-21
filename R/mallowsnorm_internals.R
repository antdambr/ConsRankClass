## -----------------------------------------------------------------------------
## [IT] .logpsi: log psi(lambda) = -log C(lambda), log della costante di
## normalizzazione. Formule chiuse, valide per qualunque M:
##   Kendall: psi = prod_{j=1}^M (1 - e^{-j lambda}) / (1 - e^{-lambda});
##   Cayley : psi = prod_{j=1}^{M-1} (1 + j e^{-lambda}).
## Per le altre distanze: log sum_d N(d) exp(-lambda d).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.logpsi <- function(lambda, M, dist) {
  switch(dist,
    kendall = {
      if (lambda == 0) return(lgamma(M + 1))
      j <- seq_len(M)
      sum(log(-expm1(-j * lambda))) - M * log(-expm1(-lambda))
    },
    cayley = {
      if (M < 2) return(0)
      sum(log1p(seq_len(M - 1) * exp(-lambda)))
    },
    {
      dd <- .distfreq(M, dist)
      .logsumexp(dd$logN - lambda * dd$d)
    })
}

## -----------------------------------------------------------------------------
## [IT] .expectd: E_lambda[d] = - d log psi / d lambda. Formule chiuse:
##   Kendall: sum_{j=1}^M [1/(e^lambda - 1) - j/(e^{j lambda} - 1)]
##            (per lambda -> 0 si usa lo sviluppo (j-1)/2 - lambda (j^2-1)/12);
##   Cayley : sum_{j=1}^{M-1} j / (e^lambda + j).
## Per le altre distanze: media della distribuzione esponenzialmente pesata.
## E_lambda[d] e' decrescente in lambda (da E_0 = media uniforme a 0).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.expectd <- function(lambda, M, dist) {
  switch(dist,
    kendall = {
      j <- seq_len(M)
      if (lambda < 1e-6) return(sum((j - 1) / 2 - lambda * (j^2 - 1) / 12))
      sum(1 / expm1(lambda) - j / expm1(j * lambda))
    },
    cayley = {
      if (M < 2) return(0)
      j <- seq_len(M - 1)
      sum(j / (exp(lambda) + j))
    },
    {
      dd <- .distfreq(M, dist)
      a <- dd$logN - lambda * dd$d
      p <- exp(a - max(a))
      sum(dd$d * p) / sum(p)
    })
}

## -----------------------------------------------------------------------------
## [IT] .vard: Var_lambda[d] = d^2 log psi / d lambda^2 (informazione di
## Fisher per lambda con il centro fissato; serve per l'errore standard).
##   Kendall: sum_j [1/(4 sinh^2(lambda/2)) - j^2/(4 sinh^2(j lambda/2))];
##   Cayley : sum_j j e^lambda / (e^lambda + j)^2.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.vard <- function(lambda, M, dist) {
  switch(dist,
    kendall = {
      j <- seq_len(M)
      if (lambda < 1e-6) return(M * (M - 1) * (2 * M + 5) / 72)
      sum(1 / (4 * sinh(lambda / 2)^2) - j^2 / (4 * sinh(j * lambda / 2)^2))
    },
    cayley = {
      if (M < 2) return(0)
      j <- seq_len(M - 1)
      sum(j * exp(lambda) / (exp(lambda) + j)^2)
    },
    {
      dd <- .distfreq(M, dist)
      a <- dd$logN - lambda * dd$d
      p <- exp(a - max(a)); p <- p / sum(p)
      m <- sum(dd$d * p)
      sum((dd$d - m)^2 * p)
    })
}

## -----------------------------------------------------------------------------
## [IT] .solvelambda: stima di massima verosimiglianza di lambda dato il
## centro. L'equazione di verosimiglianza e' E_lambda[d] = dbar (distanza
## media, eventualmente pesata), risolta con uniroot su [0, lambdamax].
## Se dbar >= media sotto l'uniforme, lambda = 0; se dbar e' troppo piccola
## (ad es. 0: tutti i dati uguali al centro) lambda = lambdamax.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.solvelambda <- function(dbar, M, dist, lambdamax = 50) {
  if (!is.finite(dbar)) return(NA_real_)
  if (dbar >= .expectd(0, M, dist)) return(0)
  if (dbar <= .expectd(lambdamax, M, dist)) return(lambdamax)
  stats::uniroot(function(l) .expectd(l, M, dist) - dbar, c(0, lambdamax),
                 tol = 1e-10)$root
}
