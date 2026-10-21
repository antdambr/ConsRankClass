#' Plot a projection of rankings
#'
#' Plots the projection of a set of rankings on two axes of the permutation
#' polytope (a projected "polytope-o-gram", Marden, 1995, Section 2.2): each
#' distinct ranking is a symbol whose area is proportional to the number of
#' judges giving it, the poles of the objects are arrows from the center, and
#' the unit circle (where the permutations lie when the projection is exact) is
#' drawn as a reference.
#'
#' @param x An object of class \code{rankproj}, as returned by
#'   \code{\link{rankproj}}.
#' @param dims The two dimensions to plot. Default \code{c(1, 2)}.
#' @param groups Optional factor (or vector) with one value for each row of the
#'   data given to \code{rankproj} (e.g. the leaves of a tree, the clusters of a
#'   mixture, or a covariate): the judges are plotted with a different symbol
#'   for each group, and a legend is added.
#' @param labels If \code{TRUE}, each distinct ranking is labelled with its
#'   order of the objects (ties in brackets). Default \code{FALSE}; useful with
#'   few distinct rankings.
#' @param maxcex Size (\code{cex}) of the symbol of the most frequent ranking.
#'   Default 3.
#' @param main Title of the plot. Default: none.
#' @param xlab,ylab Labels of the axes. Default: the dimension and its
#'   proportion of the sum of the eigenvalues.
#' @param \dots Other arguments passed to \code{\link[lattice]{xyplot}}.
#'
#' @return An object of class \code{trellis} (printed, so the plot is drawn,
#'   when not assigned).
#'
#' @seealso \code{\link{rankproj}}
#'
#' @examples
#' data(sports, package = "ConsRank")
#' plot(rankproj(sports))
#'
#' @importFrom lattice xyplot panel.points panel.text panel.abline panel.arrows panel.lines
#' @export

plot.rankproj <- function(x, dims = c(1, 2), groups = NULL, labels = FALSE,
                          maxcex = 3, main = NULL, xlab = NULL, ylab = NULL, ...) {

  k <- ncol(x$scores)                                  # dimensions available
  if (length(dims) != 2 || any(!dims %in% seq_len(k)))
    stop("dims must be two of the dimensions computed by rankproj (1 to ", k, ").")
  if (is.null(xlab)) xlab <- sprintf("Dimension %d (%.1f%%)", dims[1], 100 * x$prop[dims[1]]) # default labels
  if (is.null(ylab)) ylab <- sprintf("Dimension %d (%.1f%%)", dims[2], 100 * x$prop[dims[2]])

  d <- rankgroups(x, groups)                           # symbols: distinct rankings (and groups)
  d$x <- x$scores[d$id, dims[1]]                       # coordinates of the symbols
  d$y <- x$scores[d$id, dims[2]]
  d$cex <- maxcex * sqrt(d$w / max(d$w))               # area proportional to the number of judges
  P <- x$poles[, dims, drop = FALSE]                   # projected poles
  ng <- nlevels(d$g)                                   # number of groups
  pchs <- rep(c(1, 2, 0, 5, 6, 3, 4, 8), length.out = ng) # black and white symbols, one per group
  rl <- if (labels) rankorders(x$rankings, x$items) else NULL # labels of the distinct rankings
  lx <- x$scores[, dims[1]]                            # coordinates of the labels of the rankings
  ly <- x$scores[, dims[2]]                            # (computed here: inside the panel x is its argument)

  key <- if (ng > 1) list(space = "right", points = list(pch = pchs, col = "grey40"),
                          text = list(levels(d$g))) else NULL # legend for the groups
  tt <- seq(0, 2 * pi, length.out = 200)               # unit circle
  plt <- lattice::xyplot(y ~ x, data = d,
    xlim = c(-1.15, if (labels) 1.5 else 1.15), ylim = c(-1.15, 1.15), # the unit ball, with room for the labels
    aspect = 1, main = main, xlab = xlab, ylab = ylab, key = key,
    panel = function(x, y, ...) {
      lattice::panel.abline(h = 0, v = 0, col = "grey85") # axes through the center
      lattice::panel.lines(cos(tt), sin(tt), col = "grey60", lty = 2) # unit circle
      lattice::panel.points(d$x, d$y, pch = pchs[as.integer(d$g)], cex = d$cex, col = "grey40") # rankings (below the poles)
      if (!is.null(rl))                                # labels of the distinct rankings
        lattice::panel.text(lx, ly, rl, cex = 0.5, pos = 4, col = "grey30")
      lattice::panel.arrows(0, 0, P[, 1], P[, 2], length = 0.08, col = "black", lwd = 1.5) # poles, on top
      for (j in seq_len(nrow(P)))                      # names of the objects at the tips of the poles
        lattice::panel.text(1.08 * P[j, 1], 1.08 * P[j, 2], rownames(P)[j], cex = 0.85, fontface = 2,
                            adj = c(if (P[j, 1] >= 0) 0 else 1, 0.5))
    }, ...)
  return(plt)
}
