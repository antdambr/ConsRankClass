#' Multidimensional scaling of the judges
#'
#' Represents the judges in a low-dimensional space by multidimensional scaling
#' (\code{\link{mdscale}}) of the Kemeny distances between their rankings.
#' Judges with the same ranking share a point. Use \code{plot} to display the
#' configuration.
#'
#' @param X A n by m data matrix (or data frame), in which there are n judges
#'   and m objects to be judged. Each row is a ranking of the objects, which
#'   are represented by the columns; ties are allowed, and so are missing
#'   values (\code{NA} or \code{NaN}) for objects not ranked by a judge (e.g.
#'   top-k partial rankings), which are treated as tied in the last position.
#' @param wk Optional vector of n non-negative weights (e.g. the frequencies of
#'   the distinct rankings in the rows of \code{X}). Default \code{NULL}: every
#'   judge has weight 1.
#' @param p The number of dimensions. Default 2.
#' @param weighted If \code{TRUE} (default), the pair of distinct rankings i, j
#'   has weight \eqn{w_i w_j}{w_i w_j} in the criterion, where \eqn{w_i}{w_i} is
#'   the number (total weight) of judges giving ranking i, so the configuration
#'   represents the judges rather than the distinct rankings; if \code{FALSE}
#'   all the pairs of distinct rankings have the same weight.
#' @param \dots Other arguments passed to \code{\link{mdscale}} (e.g.
#'   \code{criterion}, \code{start}, \code{replicates}, \code{control}).
#'
#' @details The scaling is computed on the distinct rankings only, whose number
#'   is often much smaller than the number of judges, with the Kemeny distances
#'   between them (\code{ConsRank::kemenyd}). The default criterion is the one
#'   of \code{mdscale} (\code{"sstress"}). The output of \code{mdscale} is
#'   returned too, so that, e.g., its Shepard diagram can be plotted with
#'   \code{plot(x$mds, type = "shepard")}. To place the objects in the same
#'   space as the judges use \code{\link{kunfolding}}.
#'
#' @return An object of class \code{judgemds}, a list with:
#' \tabular{ll}{
#' conf \tab the configuration of the distinct rankings (one row each)\cr
#' rankings \tab the distinct rankings, in standard form (1, 2, ... for the buckets)\cr
#' weights \tab the number (total weight) of judges with each distinct ranking\cr
#' idx \tab the distinct ranking of each row of \code{X} (\code{NA} for rows with zero weight)\cr
#' wk \tab the weights of the rows of \code{X} (1 if not given)\cr
#' mds \tab the output of \code{\link{mdscale}}\cr
#' nj \tab the number of judges (the sum of the weights)\cr
#' items \tab the names of the objects
#' }
#'
#' @author Luca Coraggio \email{luca.coraggio@unina.it}
#'
#' @seealso \code{\link{plot.judgemds}}, \code{\link{mdscale}},
#'   \code{\link{rankproj}}, \code{\link{kunfolding}}
#'
#' @examples
#' data(Univranks)
#' jm <- judgemds(Univranks$rankings)
#' jm
#' plot(jm, groups = Univranks$predictors$Stud)
#'
#' @export

judgemds <- function(X, wk = NULL, p = 2, weighted = TRUE, ...) {

  R <- rankbuckets(X, wk, "judgemds")                  # distinct rankings, with weights
  nr <- nrow(R$Y)                                      # number of distinct rankings
  if (nr < 3)
    stop("judgemds() needs at least three distinct rankings.")
  D <- kemenyd(R$Y)                                    # Kemeny distances between the distinct rankings
  W <- if (weighted) outer(R$w, R$w) else NULL         # pair weights: products of the frequencies
  fit <- mdscale(D, p = p, weights = W, ...)           # multidimensional scaling
  conf <- fit$conf                                     # configuration of the distinct rankings
  dimnames(conf) <- list(NULL, paste0("Dim", seq_len(ncol(conf)))) # names of the axes

  out <- list(conf = conf, rankings = R$Y, weights = R$w, idx = R$idx, wk = R$wk,
              mds = fit, nj = sum(R$w), items = R$items) # the results
  class(out) <- c("judgemds", "ConsRankClass")         # class for print and plot
  return(out)
}

#-------------------------------------------------------------------------------
#' Print a multidimensional scaling of the judges
#'
#' Print method for objects of class \code{judgemds}.
#'
#' @param x An object of class \code{judgemds}.
#' @param \dots Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @export

print.judgemds <- function(x, ...) {

  cat("\nMultidimensional scaling of", format(x$nj), "judges (", nrow(x$rankings),
      "distinct rankings ) of", length(x$items), "objects\n") # header
  cat("Kemeny distances,", x$mds$criterion, "criterion =", format(x$mds$stress, digits = 4),
      "in", ncol(x$conf), "dimensions\n")              # fit of the scaling
  invisible(x)
}

#-------------------------------------------------------------------------------
#' Plot a multidimensional scaling of the judges
#'
#' Plots the configuration of the judges: each distinct ranking is a symbol
#' whose area is proportional to the number of judges giving it, optionally
#' with a different symbol for each group of judges.
#'
#' @param x An object of class \code{judgemds}, as returned by
#'   \code{\link{judgemds}}.
#' @param dims The two dimensions to plot. Default \code{c(1, 2)}.
#' @param groups Optional factor (or vector) with one value for each row of the
#'   data given to \code{judgemds} (e.g. the leaves of a tree, the clusters of a
#'   mixture, or a covariate): the judges are plotted with a different symbol
#'   for each group, and a legend is added.
#' @param labels If \code{TRUE}, each distinct ranking is labelled with its
#'   order of the objects (ties in brackets). Default \code{FALSE}.
#' @param maxcex Size (\code{cex}) of the symbol of the most frequent ranking.
#'   Default 3.
#' @param main Title of the plot. Default: none.
#' @param xlab,ylab Labels of the axes. Defaults: "Dimension 1" and
#'   "Dimension 2" (or the dimensions chosen).
#' @param \dots Other arguments passed to \code{\link[lattice]{xyplot}}.
#'
#' @return An object of class \code{trellis} (printed, so the plot is drawn,
#'   when not assigned).
#'
#' @seealso \code{\link{judgemds}}
#'
#' @examples
#' data(EVS)
#' plot(judgemds(EVS$rankings), labels = TRUE)
#'
#' @importFrom lattice xyplot panel.points panel.text panel.abline
#' @export

plot.judgemds <- function(x, dims = c(1, 2), groups = NULL, labels = FALSE,
                          maxcex = 3, main = NULL, xlab = NULL, ylab = NULL, ...) {

  p <- ncol(x$conf)                                    # dimensions available
  if (length(dims) != 2 || any(!dims %in% seq_len(p)))
    stop("dims must be two of the dimensions computed by judgemds (1 to ", p, ").")
  if (is.null(xlab)) xlab <- paste("Dimension", dims[1]) # default labels
  if (is.null(ylab)) ylab <- paste("Dimension", dims[2])
  d <- rankgroups(x, groups)                           # symbols: distinct rankings (and groups)
  d$x <- x$conf[d$id, dims[1]]                         # coordinates of the symbols
  d$y <- x$conf[d$id, dims[2]]
  d$cex <- maxcex * sqrt(d$w / max(d$w))               # area proportional to the number of judges
  ng <- nlevels(d$g)                                   # number of groups
  pchs <- rep(c(1, 2, 0, 5, 6, 3, 4, 8), length.out = ng) # black and white symbols, one per group
  rl <- if (labels) rankorders(x$rankings, x$items) else NULL # labels of the distinct rankings
  lx <- x$conf[, dims[1]]                              # coordinates of the labels
  ly <- x$conf[, dims[2]]
  key <- if (ng > 1) list(space = "right", points = list(pch = pchs, col = "grey40"),
                          text = list(levels(d$g))) else NULL # legend for the groups

  rx <- range(lx); ry <- range(ly)                     # ranges of the coordinates
  xl <- rx + c(-0.08, if (labels) 0.35 else 0.08) * diff(rx) # room on the right for the labels
  yl <- ry + c(-0.08, 0.08) * diff(ry)                 # some room above and below
  plt <- lattice::xyplot(y ~ x, data = d, aspect = "iso", # same scale on the two axes
    xlim = xl, ylim = yl,
    main = main, xlab = xlab, ylab = ylab, key = key,
    panel = function(x, y, ...) {
      lattice::panel.abline(h = 0, v = 0, col = "grey85") # axes through the center
      lattice::panel.points(d$x, d$y, pch = pchs[as.integer(d$g)], cex = d$cex, col = "grey40") # rankings
      if (!is.null(rl))                                # labels of the distinct rankings
        lattice::panel.text(lx, ly, rl, cex = 0.5, pos = 4, col = "grey20")
    }, ...)
  return(plt)
}

#-------------------------------------------------------------------------------
## Symbols of the plots of rankproj and judgemds: one per distinct ranking, or
## one per distinct ranking and group when groups (one value per row of the
## data) are given. Returns a data frame with the distinct ranking (id), the
## number (total weight) of judges (w) and the group (g, a factor).

rankgroups <- function(x, groups) {

  if (is.null(groups))                                 # one symbol per distinct ranking
    return(data.frame(id = seq_len(nrow(x$rankings)), w = x$weights, g = factor(1)))
  if (length(groups) != length(x$idx))
    stop("groups must have one value for each row of the data (", length(x$idx), ").")
  ok <- !is.na(x$idx) & !is.na(groups)                 # rows used
  g <- factor(groups)                                  # groups as a factor
  agg <- stats::aggregate(list(w = x$wk[ok]), by = list(id = x$idx[ok], g = g[ok]), FUN = sum) # weight per ranking and group
  agg <- agg[agg$w > 0, , drop = FALSE]                # combinations with some judges
  data.frame(id = agg$id, w = agg$w, g = factor(agg$g, levels = levels(g)))
}

## Orders of the objects of rankings in standard form, as labels: objects from
## the best, ties in brackets, e.g. "A>(B C)>D".

rankorders <- function(Y, items) {
  apply(Y, 1, function(r) {
    b <- split(items, r)                               # objects in each bucket, best first
    paste(sapply(b, function(o) if (length(o) > 1) paste0("(", paste(o, collapse = " "), ")") else o),
          collapse = ">")
  })
}
