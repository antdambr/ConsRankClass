## -----------------------------------------------------------------------------
## [IT] distfreq (esportata): restituisce la distribuzione di frequenza della
## distanza sotto il modello uniforme, come data.frame (d, count, logcount).
## Analoga a PerMallows::count.perms e a MSmix::spear_dist_distr.
## -----------------------------------------------------------------------------
#' Distribution of a ranking distance under the uniform model
#'
#' Number of permutations of M objects at each possible distance from a
#' fixed ranking (by right-invariance it does not depend on that ranking).
#' It is the basis of the normalising constant of the Mallows model.
#'
#' @param M Number of objects.
#' @param dist The distance (see \code{\link{rankdist}}).
#'
#' @details Exact computation: Mahonian numbers (Kendall), Stirling numbers of
#' the first kind (Cayley), binomial coefficients times derangement numbers
#' (Hamming), sums of squared numbers of standard Young tableaux via the
#' Robinson-Schensted correspondence (Ulam, \eqn{M \le 50}), dynamic
#' programming over subsets (footrule and Spearman, \eqn{M \le 14}).
#'
#' @return A data frame with the distance values \code{d}, the number of
#'   permutations \code{count} and its logarithm \code{logcount}.
#'
#' @examples
#' distfreq(5, "kendall")
#' sum(distfreq(6, "ulam")$count) == factorial(6)
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mallowsnorm}}, \code{\link{rankdist}}
#' @export
distfreq <- function(M, dist = "kendall") {
  dist <- .matchdist(dist)
  dd <- .distfreq(M, dist)
  data.frame(d = dd$d, count = dd$N, logcount = dd$logN)
}
