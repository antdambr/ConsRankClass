#' Options of the algorithm of mdscale
#'
#' Options that control the iterative algorithm of \code{mdscale}. Names and default
#' values are those of the options of the MATLAB function \code{mdscale}.
#'
#' @param maxiter the maximum number of iterations. Default 200.
#' @param tolfun the termination tolerance for the criterion, for its relative
#'   change and for the norm of its gradient. Default 1e-4.
#' @param tolx the termination tolerance for the norm of the step. Default 1e-4.
#' @param display the amount of output printed: \code{"off"} (default),
#'   \code{"final"} (iterations and final criterion of each replicate) or
#'   \code{"iter"} (one line per iteration).
#' @param window the number of iterations over which the relative change of the
#'   criterion is measured: the iterations stop when the criterion has decreased by
#'   less than \code{tolfun} (in relative terms) over the last \code{window}
#'   iterations. Default 10. With \code{window = 1} the rule is the one of MATLAB
#'   (relative change over one iteration); a longer window avoids stopping at an
#'   iteration that happens to make little progress far from the minimum.
#'
#' @return A list with the options.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}, Luca Coraggio \email{luca.coraggio@unina.it}
#'
#' @seealso \code{\link{mdscale}}
#'
#' @examples
#' ctl <- mdscalecontrol(maxiter = 500, tolfun = 1e-6)
#'
#' @export

mdscalecontrol <- function(maxiter = 200, tolfun = 1e-4, tolx = 1e-4,
                           display = c("off", "final", "iter"), window = 10) {
  display <- match.arg(display)                        # check the display option
  if (length(window) != 1 || window < 1 || window != round(window)) # a positive integer
    stop("window must be a positive integer.")
  if (length(maxiter) != 1 || maxiter < 0 || maxiter != round(maxiter)) # a non-negative integer
    stop("maxiter must be a non-negative integer.")
  if (length(tolfun) != 1 || !(tolfun >= 0))           # a non-negative number
    stop("tolfun must be a non-negative number.")
  if (length(tolx) != 1 || !(tolx >= 0))               # a non-negative number
    stop("tolx must be a non-negative number.")
  list(maxiter = maxiter, tolfun = tolfun, tolx = tolx, display = display, # the options
       window = as.integer(window))
}
