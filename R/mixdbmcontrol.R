## -----------------------------------------------------------------------------
## [IT] mixdbmcontrol (esportata): opzioni dell'algoritmo, come ccacontrol()
## per cca(). Vedi la documentazione per il significato di ciascuna.
## -----------------------------------------------------------------------------
#' Control parameters for mixdbm
#'
#' Options of the EM algorithm of \code{\link{mixdbm}} and
#' \code{\link{mixdbmselect}}.
#'
#' @param nstart Number of random starts (Murphy and Martin use 30).
#' @param init Starting values: \code{"mixed"} (default, alternates
#'   \code{"partition"} and \code{"centers"}), \code{"partition"} (random soft
#'   partition of the data followed by an M-step), \code{"centers"} (centers
#'   drawn from the observed rankings) or \code{"sem"} (random partition, then
#'   \code{semiter} iterations of the stochastic EM of Celeux and Diebolt,
#'   1985; the best visited parameters start the EM).
#' @param semiter Number of SEM iterations per start when \code{init = "sem"}.
#' @param maxiter Maximum number of EM iterations per start.
#' @param tol Relative tolerance on the log-likelihood for convergence.
#' @param centerdf Degrees of freedom charged to each central ranking in BIC
#'   and ICL (default 0: centers are discrete and not counted; the paper
#'   does not state it).
#' @param centermethod Method for the M-step of the centers (see
#'   \code{\link{mallowscenter}}).
#' @param lambdamax Upper bound for the precisions.
#' @param minsize A component is \emph{degenerate} when its precision reaches
#'   \code{lambdamax} and its expected size n p_g is below \code{minsize}
#'   (default 5): it is concentrated on a few copies of one ranking and is a
#'   spurious maximum at the boundary (lambda = Inf) of the parameter space.
#'   The best non-degenerate start is kept (see \code{\link{mixdbm}}).
#'   \code{minsize = 0} switches the check off.
#' @param verbose If \code{TRUE} prints the log-likelihood reached by each start.
#'
#' @return A list with the options.
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' @seealso \code{\link{mixdbm}}
#' @export
mixdbmcontrol <- function(nstart = 30,
                          init = c("mixed", "partition", "centers", "sem"),
                          semiter = 100, maxiter = 1000, tol = 1e-10,
                          centerdf = 0,
                          centermethod = c("auto", "exhaustive", "local", "bb"),
                          lambdamax = 50, minsize = 5, verbose = FALSE) {
  list(nstart = nstart, init = match.arg(init), semiter = semiter,
       maxiter = maxiter, tol = tol, centerdf = centerdf,
       centermethod = match.arg(centermethod), lambdamax = lambdamax,
       minsize = minsize, verbose = verbose)
}
