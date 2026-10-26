#' Plot a matrix of pairwise preferences
#'
#' Heatmap, in grey levels, of a matrix of pairwise preferences: the cell in
#' row a and column b is the proportion of judges preferring a to b (darker:
#' larger).
#'
#' @param x An object of class \code{pairmatrix}, as returned by
#'   \code{\link{pairmatrix}}.
#' @param type The proportions displayed: \code{"half"} (default; ties count one
#'   half, so the cells (a, b) and (b, a) sum to one), \code{"strict"} (strict
#'   preferences only) or \code{"ties"} (proportion of ties, a symmetric
#'   matrix).
#' @param order The order of the objects: \code{"score"} (default, from the
#'   largest normalised Borda score, so the darker cells are above the
#'   diagonal) or \code{"none"} (the order of the columns of the data).
#' @param values If \code{TRUE} (default), the proportions are written in the
#'   cells.
#' @param digits Number of decimal digits of the values. Default 2.
#' @param main Title of the plot. Default: none.
#' @param xlab,ylab Labels of the axes. Defaults: none.
#' @param \dots Other arguments passed to \code{\link[lattice]{levelplot}}.
#'
#' @return An object of class \code{trellis} (printed, so the plot is drawn,
#'   when not assigned).
#'
#' @seealso \code{\link{pairmatrix}}
#'
#' @examples
#' data(Univranks)
#' plot(pairmatrix(Univranks$rankings))
#'
#' @importFrom lattice levelplot panel.levelplot panel.text
#' @importFrom grDevices grey
#' @author Luca Coraggio \email{luca.coraggio@unina.it}
#'
#' @export

plot.pairmatrix <- function(x, type = c("half", "strict", "ties"),
                            order = c("score", "none"), values = TRUE,
                            digits = 2, main = NULL, xlab = NULL, ylab = NULL, ...) {

  type <- match.arg(type)                              # proportions displayed
  order <- match.arg(order)                            # order of the objects
  M <- switch(type, half = x$prop, strict = x$strict, ties = x$tieprop) # the matrix
  m <- length(x$items)                                 # number of objects
  ord <- if (order == "score") order(x$score, decreasing = TRUE) else seq_len(m) # best object first
  M <- M[ord, ord, drop = FALSE]                       # rows and columns in that order
  d <- data.frame(col = rep(seq_len(m), each = m),     # column of the cell
                  row = rep(m:1, times = m),           # row from the bottom: first object at the top
                  z = as.vector(M))                    # proportion of the cell
  labs <- x$items[ord]                                 # names of the objects in that order

  plt <- lattice::levelplot(z ~ col * row, data = d,
    at = seq(0, 1, length.out = 21),                   # 20 grey levels from 0 to 1
    col.regions = grDevices::grey(seq(1, 0.1, length.out = 20)), # white (0) to almost black (1)
    scales = list(x = list(at = seq_len(m), labels = labs, rot = 45), # objects on both axes
                  y = list(at = m:1, labels = labs)),
    main = main, xlab = xlab, ylab = ylab,
    aspect = 1,                                        # square cells
    panel = function(x, y, z, ...) {
      lattice::panel.levelplot(x, y, z, ...)           # the grey cells (diagonal left blank)
      if (values) {                                    # proportions written in the cells
        ok <- !is.na(z)                                # not on the diagonal
        lattice::panel.text(x[ok], y[ok], formatC(z[ok], format = "f", digits = digits),
                            cex = 0.6, col = ifelse(z[ok] > 0.55, "white", "black")) # readable on dark cells
      }
    }, ...)
  return(plt)
}
