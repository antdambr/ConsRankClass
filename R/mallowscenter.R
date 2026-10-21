## -----------------------------------------------------------------------------
## [IT] mallowscenter (esportata): mediana (eventualmente pesata) di un insieme
## di ranking rispetto a una distanza; e' la stima di massima verosimiglianza
## del centro del modello di Mallows. Con dist = "kendall" e' il consenso di
## Kemeny sulle sole permutazioni (come consrank(..., full = TRUE)).
## -----------------------------------------------------------------------------
#' Median ranking with respect to a distance
#'
#' Finds the ranking R minimising \eqn{\sum_i w_i d(r_i, R)}{sum_i w_i d(r_i, R)}
#' over all permutations: it is the maximum likelihood estimate of the
#' center of the Mallows model.
#'
#' @param X A n by M matrix of complete rankings.
#' @param dist The distance (see \code{\link{rankdist}}).
#' @param weights Optional non-negative weights (default: all equal to 1).
#' @param centermethod \code{"auto"} (default), \code{"exhaustive"},
#'   \code{"local"} or \code{"bb"} (see Details).
#'
#' @details Exact methods: weighted Borda count (\code{"spearmansq"}),
#' assignment problem solved with the Hungarian algorithm (\code{"footrule"},
#' \code{"hamming"}), Emond and Mason (2002) branch-and-bound of
#' \code{\link[ConsRank]{consrank}} (Kendall, used by \code{"auto"} when
#' M > 8), exhaustive search (\code{"auto"}: Kendall with M <= 8, the
#' other distances with M <= 7). Otherwise a local search (swaps and moves of
#' one object) from several starting points, which may stop at a local
#' optimum.
#'
#' @return The median ranking (vector named after the objects).
#'
#' @examples
#' data(Irish)
#' ## complete rankings only (rankings with ties are not supported yet)
#' X <- Irish$rankings[apply(Irish$rankings, 1, function(r) all(sort(r) == 1:5)), ]
#' mallowscenter(X, dist = "kendall")
#' mallowscenter(X, dist = "footrule")
#'
#' @references Emond, E.J. and Mason, D.W. (2002). A new rank correlation
#'   coefficient with application to the consensus ranking problem. Journal
#'   of Multi-Criteria Decision Analysis, 11, 17-28.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mallows}}, \code{\link[ConsRank]{consrank}}
#' @export
mallowscenter <- function(X, dist = "kendall", weights = NULL,
                          centermethod = c("auto", "exhaustive", "local", "bb")) {
  dist <- .matchdist(dist)
  centermethod <- match.arg(centermethod)
  X <- .checkrank(X, "X")
  ag <- .aggregaterank(X, weights)
  ctx <- .makectx(ag$U, ncol(X), dist, centermethod)
  R <- .centerestimate(ag$U, ag$w, dist, ctx)
  names(R) <- .objnames(X)
  R
}
