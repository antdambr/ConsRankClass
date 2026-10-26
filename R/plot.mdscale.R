#' Plot method for mdscale
#'
#' Plots the configuration of an \code{mdscale} solution, or its Shepard diagram.
#'
#' @param x an object of the class \code{mdscale}, as returned by \code{mdscale}.
#' @param type \code{"configuration"} (default) for the points in two of the
#'   dimensions, or \code{"shepard"} for the Shepard diagram: the distances (points)
#'   and, for the non-metric criteria, the disparities (line) against the
#'   dissimilarities.
#' @param dims the two dimensions to plot. Default \code{c(1, 2)}.
#' @param labels the labels of the points. Default: the names of the rows of the
#'   configuration, or the numbers of the points if there are no names. Use
#'   \code{labels = FALSE} to plot the points without labels.
#' @param col the colour of the points (configuration) or of the distances
#'   (Shepard diagram). Default \code{"black"}.
#' @param cex the size of the labels and of the points. Default 0.8.
#' @param main the title of the plot. Default: none.
#' @param \dots other graphical parameters passed to \code{plot}.
#'
#' @details In the configuration plot the two axes have the same scale, so that the
#' distances in the plot are the distances of the solution.
#'
#' @return No return value, called for its side effect (a plot).
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}, Luca Coraggio \email{luca.coraggio@unina.it}
#'
#' @seealso \code{\link{mdscale}}
#'
#' @examples
#' data(Univranks)
#' out <- mdscale(ConsRank::kemenyd(Univranks$rankings))
#' plot(out)
#' plot(out, type = "shepard")
#'
#' @importFrom graphics plot text points lines abline
#' @export

plot.mdscale <- function(x, type = c("configuration", "shepard"), dims = c(1, 2),
                         labels = NULL, col = "black", cex = 0.8, main = NULL, ...) {
  type <- match.arg(type)                              # check the type of plot
  if (is.null(main)) main <- ""                        # no title by default
  if (type == "configuration") {                       # ---- configuration
    Y <- x$conf                                        # the configuration
    if (ncol(Y) < 2) {                                 # one dimension: points on a line
      Y <- cbind(Y, 0)                                 # second coordinate zero
      dims <- c(1, 2)                                  # the two columns
    }
    if (length(dims) != 2 || any(dims > ncol(Y)) || any(dims < 1)) # two valid dimensions
      stop("dims must be two of the dimensions of the solution.")
    if (is.null(labels)) {                             # default labels
      labels <- rownames(Y)                            # the names of the points
      if (is.null(labels)) labels <- seq_len(nrow(Y))  # or their numbers
    }
    lim <- range(Y[, dims])                            # same range on both axes
    op <- graphics::par(pty = "s")                     # square plotting region
    on.exit(graphics::par(op))                         # restore it at the end
    plot(Y[, dims[1]], Y[, dims[2]], type = if (isFALSE(labels)) "p" else "n", # points or room for labels
         xlim = lim, ylim = lim, asp = 1, col = col, cex = cex, pch = 20,
         xlab = paste("Dimension", dims[1]), ylab = paste("Dimension", dims[2]),
         main = main, ...)
    abline(h = 0, v = 0, col = "gray80", lty = 3)      # axes through the origin
    if (!isFALSE(labels))                              # labels at the points
      text(Y[, dims[1]], Y[, dims[2]], labels, col = col, cex = cex)
  } else {                                             # ---- Shepard diagram
    up <- upper.tri(x$dissimilarities)                 # each pair once
    D <- x$dissimilarities[up]                         # dissimilarities
    d <- x$distances[up]                               # distances
    dh <- x$disparities[up]                            # disparities
    ok <- !is.na(D) & !is.na(dh)                       # pairs used by the criterion
    if (!is.null(x$weights)) ok <- ok & x$weights[up] > 0 # positive weights only
    plot(D[ok], d[ok], col = col, cex = cex * 0.6, pch = 1, # distances against dissimilarities
         xlab = "Dissimilarities", ylab = "Distances / Disparities", main = main, ...)
    if (x$criterion %in% c("stress", "sstress")) {     # non-metric: the monotone regression
      o <- order(D[ok], dh[ok])                        # in order of dissimilarity
      lines(D[ok][o], dh[ok][o], col = "red", type = "s") # step function of the disparities
    } else {                                           # metric: the line d = D
      abline(0, 1, col = "red", lty = 2)
    }
  }
  invisible(NULL)                                      # nothing to return
}
