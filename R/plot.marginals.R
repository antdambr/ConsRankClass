#' Plot a marginal matrix
#'
#' Bubble plot of the marginal matrix of a set of rankings: the objects are in
#' the rows and the positions in the columns, and the area of each circle is
#' proportional to the proportion of judges who put the object in that
#' position. The mean rank of each object can be marked on its row.
#'
#' @param x An object of class \code{marginals}, as returned by
#'   \code{\link{marginals}}.
#' @param order The order of the objects (rows, from the top): \code{"meanrank"}
#'   (default, from the smallest mean rank) or \code{"none"} (the order of the
#'   columns of the data).
#' @param meanrank If \code{TRUE} (default), the mean rank of each object is
#'   marked with a filled diamond on its row.
#' @param values If \code{TRUE}, the proportions (in percent) are written in the
#'   circles. Default \code{FALSE}.
#' @param maxradius Radius of the circle of the largest proportion, as a
#'   fraction of the distance between two positions. Default 0.45.
#' @param main Title of the plot. Default: none.
#' @param xlab,ylab Labels of the axes. Defaults: "Position" and none.
#' @param \dots Other arguments passed to \code{\link[lattice]{xyplot}}.
#'
#' @return An object of class \code{trellis} (printed, so the plot is drawn,
#'   when not assigned).
#'
#' @seealso \code{\link{marginals}}
#'
#' @examples
#' data(Univranks)
#' plot(marginals(Univranks$rankings), values = TRUE)
#'
#' @importFrom lattice xyplot panel.points panel.text panel.abline
#' @importFrom grid grid.circle unit gpar
#' @author Luca Coraggio \email{luca.coraggio@unina.it}
#'
#' @export

plot.marginals <- function(x, order = c("meanrank", "none"), meanrank = TRUE,
                           values = FALSE, maxradius = 0.45, main = NULL,
                           xlab = "Position", ylab = NULL, ...) {

  order <- match.arg(order)                            # order of the objects
  m <- length(x$items)                                 # number of objects
  ord <- if (order == "meanrank") order(x$meanrank) else seq_len(m) # objects from the top
  P <- x$prop[ord, , drop = FALSE]                     # proportions in that order
  d <- data.frame(pos = rep(seq_len(m), each = m),     # position (column of P)
                  obj = rep(m:1, times = m),           # row from the bottom: first object at the top
                  p = as.vector(P))                    # proportion of the cell
  mr <- x$meanrank[ord]                                # mean ranks in the same order
  pmax <- max(d$p)                                     # largest proportion

  plt <- lattice::xyplot(obj ~ pos, data = d,
    xlim = c(0.4, m + 0.6), ylim = c(0.4, m + 0.6),    # room for the circles
    scales = list(x = list(at = seq_len(m)),           # positions 1, ..., m
                  y = list(at = m:1, labels = x$items[ord])), # names of the objects
    main = main, xlab = xlab, ylab = ylab,
    aspect = 1,                                        # square cells
    panel = function(x, y, ...) {
      lattice::panel.abline(h = seq_len(m), v = seq_len(m), col = "grey90") # light grid
      ok <- d$p > 0                                    # cells with some judges
      r <- maxradius * sqrt(d$p[ok] / pmax)            # radius: area proportional to the proportion
      grid::grid.circle(x[ok], y[ok], r = grid::unit(r, "native"), default.units = "native",
                        gp = grid::gpar(col = "black", fill = NA)) # circles in data units
      if (values)                                      # proportions in percent
        lattice::panel.text(x[ok], y[ok], round(100 * d$p[ok]), cex = 0.6)
      if (meanrank)                                    # mean rank of each object on its row
        lattice::panel.points(mr, m:1, pch = 23, cex = 1, col = "black", fill = "black") # filled diamond
    }, ...)
  return(plt)
}
