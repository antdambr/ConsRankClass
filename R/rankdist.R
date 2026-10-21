## -----------------------------------------------------------------------------
## [IT] rankdist: funzione generale per le distanze. Tre usi:
##   rankdist(X)       -> matrice n x n delle distanze tra le righe di X;
##   rankdist(X, y)    -> vettore delle distanze delle righe di X da y
##                        (un numero se anche X e' un vettore);
##   rankdist(X, Y)    -> matrice nrow(X) x nrow(Y) delle distanze incrociate.
## Le funzioni kendalld(), cayleyd(), ... sono scorciatoie di rankdist().
## -----------------------------------------------------------------------------
#' Distances between complete rankings
#'
#' Computes distances between complete rankings (no ties) for several
#' metrics on permutations. \code{rankdist} is the general function;
#' \code{kendalld}, \code{cayleyd}, \code{hammingd}, \code{ulamd},
#' \code{footruled} and \code{spearmand} are shortcuts for a given metric.
#'
#' @param X A n by M matrix of complete rankings (rank of object j in column
#'   j, 1 = most preferred), or a single ranking (vector).
#' @param Y Optional. \code{NULL} (default): the n by n matrix of distances
#'   between the rows of \code{X} is returned; a single ranking (vector): the
#'   distances of each row of \code{X} from \code{Y}; a matrix of rankings:
#'   the \code{nrow(X)} by \code{nrow(Y)} matrix of distances.
#' @param dist The distance: one of \code{"kendall"}, \code{"cayley"},
#'   \code{"hamming"}, \code{"ulam"}, \code{"footrule"}, \code{"spearman"}
#'   (square root of the sum of squared rank differences, as in Murphy and
#'   Martin, 2003) and \code{"spearmansq"} (sum of squared rank differences).
#' @param squared For \code{spearmand}: if \code{TRUE} the sum of squared
#'   rank differences, without square root.
#'
#' @details With \eqn{r} and \eqn{s} two rankings of M objects:
#' \itemize{
#'  \item Kendall: number of pairs of objects ranked in opposite order (on
#'   complete rankings it is one half of the Kemeny distance, see
#'   \code{kemenyd} in \pkg{ConsRank});
#'  \item Cayley: minimum number of transpositions (swaps of two objects)
#'   turning \eqn{r} into \eqn{s}, i.e. M minus the number of cycles;
#'  \item Hamming: number of objects with different rank;
#'  \item Ulam: minimum number of moves of one object to another position of
#'   the ordering, i.e. M minus the largest number of objects ranked in the
#'   same relative order by \eqn{r} and \eqn{s};
#'  \item Spearman footrule: \eqn{\sum_j |r_j - s_j|};
#'  \item Spearman: \eqn{\sqrt{\sum_j (r_j - s_j)^2}}; \code{"spearmansq"}:
#'   \eqn{\sum_j (r_j - s_j)^2}.
#' }
#' All these distances are right-invariant (they do not depend on the
#' labelling of the objects). Ties are not allowed.
#'
#' @return A matrix of distances, a vector, or a single number (see \code{Y}).
#'
#' @examples
#' r <- c(1, 2, 3, 4, 5)
#' s <- c(2, 1, 5, 3, 4)
#' kendalld(r, s)
#' cayleyd(r, s)
#' ulamd(r, s)
#' spearmand(r, s, squared = TRUE)
#' data(Irish)
#' ## complete rankings only (rankings with ties are not supported yet)
#' X <- Irish$rankings[apply(Irish$rankings, 1, function(r) all(sort(r) == 1:5)), ]
#' D <- rankdist(X[1:10, ], dist = "footrule")
#'
#' @references Diaconis, P. (1988). Group Representations in Probability and
#'   Statistics. IMS Lecture Notes - Monograph Series, 11.
#' @references Critchlow, D.E. (1985). Metric Methods for Analyzing
#'   Partially Ranked Data. Lecture Notes in Statistics 34, Springer.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mallows}}, \code{\link{mixdbm}}
#' @export
rankdist <- function(X, Y = NULL, dist = "kendall") {
  dist <- .matchdist(dist)
  xvec <- is.null(dim(X))
  X <- .checkrank(X, "X")
  M <- ncol(X)
  if (is.null(Y)) {                                  # n x n
    n <- nrow(X)
    D <- matrix(0, n, n, dimnames = list(rownames(X), rownames(X)))
    if (n > 1) for (i in 1:(n - 1)) {
      di <- .dist2center(X[(i + 1):n, , drop = FALSE], X[i, ], dist)
      D[(i + 1):n, i] <- di
      D[i, (i + 1):n] <- di
    }
    return(D)
  }
  yvec <- is.null(dim(Y))
  Y <- .checkrank(Y, "Y")
  if (ncol(Y) != M) stop("'X' and 'Y' must rank the same number of objects")
  if (yvec) {                                        # distanze da un ranking
    d <- .dist2center(X, Y[1, ], dist)
    names(d) <- if (!xvec) rownames(X)
    return(d)
  }
  D <- vapply(seq_len(nrow(Y)), function(k) .dist2center(X, Y[k, ], dist),
              numeric(nrow(X)))
  D <- matrix(D, nrow(X), nrow(Y), dimnames = list(rownames(X), rownames(Y)))
  D
}

## [IT] Scorciatoie: una funzione per ogni distanza (stessi argomenti X, Y di
## rankdist). spearmand ha in piu' l'argomento squared.

#' @rdname rankdist
#' @export
kendalld <- function(X, Y = NULL) rankdist(X, Y, "kendall")

#' @rdname rankdist
#' @export
cayleyd <- function(X, Y = NULL) rankdist(X, Y, "cayley")

#' @rdname rankdist
#' @export
hammingd <- function(X, Y = NULL) rankdist(X, Y, "hamming")

#' @rdname rankdist
#' @export
ulamd <- function(X, Y = NULL) rankdist(X, Y, "ulam")

#' @rdname rankdist
#' @export
footruled <- function(X, Y = NULL) rankdist(X, Y, "footrule")

#' @rdname rankdist
#' @export
spearmand <- function(X, Y = NULL, squared = FALSE)
  rankdist(X, Y, if (squared) "spearmansq" else "spearman")
