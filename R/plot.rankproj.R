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
#' @param labpos Position of the names of the objects with respect to the tips
#'   of their poles, coded as in \code{\link[graphics]{text}}: 1 below, 2 left,
#'   3 above, 4 right. A single value for all the objects, one value for each
#'   object, or a named vector for some objects only (e.g.
#'   \code{c(Tennis = 1)}). Default \code{NULL}: left of the tip for the poles
#'   pointing left, right of it for the others.
#' @param laboffset Shifts of the names of the objects, to separate
#'   overlapping names: a named list of vectors \code{c(dx, dy)} in the units
#'   of the plot (the unit circle has radius 1), e.g.
#'   \code{list(Tennis = c(0, -0.08))}, or a matrix with two columns and one
#'   row for each object. Default \code{NULL} (no shift).
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
#' # names of the objects moved by hand
#' plot(rankproj(sports), labpos = c(Tennis = 2, Cycling = 3),
#'      laboffset = list(Swimming = c(0, -0.06)))
#'
#' @importFrom lattice xyplot panel.points panel.text panel.abline panel.arrows panel.lines
#' @author Luca Coraggio \email{luca.coraggio@unina.it}
#'
#' @export

plot.rankproj <- function(x, dims = c(1, 2), groups = NULL, labels = FALSE,
                          maxcex = 3, main = NULL, xlab = NULL, ylab = NULL,
                          labpos = NULL, laboffset = NULL, ...) {

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
  lp <- poslabels(labpos, rownames(P))                 # position of the names (NA = default)
  off <- offlabels(laboffset, rownames(P))             # shifts of the names (0 = none)
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
      for (j in seq_len(nrow(P))) {                    # names of the objects at the tips of the poles
        if (is.na(lp[j]))                              # default: left or right of the tip
          lattice::panel.text(1.08 * P[j, 1] + off[j, 1], 1.08 * P[j, 2] + off[j, 2], rownames(P)[j],
                              cex = 0.85, fontface = 2, adj = c(if (P[j, 1] >= 0) 0 else 1, 0.5))
        else                                           # position chosen by the user
          lattice::panel.text(1.08 * P[j, 1] + off[j, 1], 1.08 * P[j, 2] + off[j, 2], rownames(P)[j],
                              cex = 0.85, fontface = 2, pos = lp[j], offset = 0.3)
      }
    }, ...)
  return(plt)
}

# position of the names of the poles (internal): NA = default, else 1 to 4
poslabels <- function(labpos, items) {
  lp <- rep(NA_integer_, length(items))                # default for all the objects
  if (is.null(labpos)) return(lp)
  if (!is.numeric(labpos) || any(!labpos %in% 1:4))
    stop("labpos must contain the values 1 (below), 2 (left), 3 (above) or 4 (right).")
  if (!is.null(names(labpos))) {                       # named vector: only some objects
    bad <- setdiff(names(labpos), items)
    if (length(bad)) stop("labpos names objects that do not exist: ", paste(bad, collapse = ", "), ".")
    lp[match(names(labpos), items)] <- labpos
  } else if (length(labpos) == 1) {                    # one value for all the objects
    lp[] <- labpos
  } else if (length(labpos) == length(items)) {        # one value for each object
    lp <- as.integer(labpos)
  } else stop("labpos must have length 1, one value for each object, or names of objects.")
  lp
}

# shifts of the names of the poles (internal): matrix with one row per object
offlabels <- function(laboffset, items) {
  off <- matrix(0, length(items), 2, dimnames = list(items, NULL)) # no shift
  if (is.null(laboffset)) return(off)
  if (is.list(laboffset)) {                            # named list of c(dx, dy)
    bad <- setdiff(names(laboffset), items)
    if (is.null(names(laboffset)) || length(bad))
      stop("laboffset must be a list named with the names of the objects, e.g. list(A = c(0, -0.05)).")
    for (nm in names(laboffset)) {
      v <- laboffset[[nm]]
      if (!is.numeric(v) || length(v) != 2) stop("each element of laboffset must be a vector c(dx, dy).")
      off[nm, ] <- v
    }
  } else if (is.matrix(laboffset) && ncol(laboffset) == 2 && nrow(laboffset) == length(items)) {
    off[] <- laboffset                                 # one row for each object
  } else stop("laboffset must be a named list of c(dx, dy) or a matrix with two columns and one row for each object.")
  off
}
