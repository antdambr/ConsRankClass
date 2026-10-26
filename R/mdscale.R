#' Metric and non-metric multidimensional scaling
#'
#' Metric and non-metric multidimensional scaling of a dissimilarity matrix, with the
#' same goodness-of-fit criteria, options and output as the function \code{mdscale}
#' of MATLAB (Statistics and Machine Learning Toolbox). It can be used on its own
#' (for instance on a matrix of Kemeny distances between rankings) and as the MDS
#' engine of \code{kunfolding}.
#'
#' @param D the dissimilarities: a symmetric N by N matrix with zero diagonal and
#'   non-negative entries, an object of class \code{dist}, or a vector with the
#'   N(N-1)/2 elements of the upper triangle in the order of \code{dist} (the same
#'   order as the MATLAB function \code{pdist}). \code{NA} entries are missing
#'   values and are ignored. A symmetric matrix with ones along the diagonal and all
#'   other elements smaller than one is treated as a similarity matrix S and
#'   transformed into the dissimilarities \code{sqrt(1 - S)}.
#' @param p the number of dimensions of the solution. Default \code{p = 2}. It is
#'   ignored when \code{start} is a matrix or an array.
#' @param criterion the goodness-of-fit criterion to minimise, which also sets the
#'   kind of scaling. Non-metric scaling: \code{"sstress"} (default; squared stress,
#'   normalised by the sum of fourth powers of the distances) or \code{"stress"}
#'   (Kruskal's stress-1, normalised by the sum of squared distances; the default of
#'   MATLAB, and the kind of criterion minimised by \code{smacof}). Metric scaling:
#'   \code{"metricstress"} (stress normalised by the sum of squared dissimilarities),
#'   \code{"metricsstress"} (squared stress normalised by the sum of fourth powers of
#'   the dissimilarities), \code{"sammon"} (Sammon's nonlinear mapping; the
#'   dissimilarities must be strictly positive) or \code{"strain"} (the criterion of
#'   classical scaling). See Details.
#' @param weights optional non-negative weights for the dissimilarities, in the same
#'   form as \code{D}. Pairs with zero weight are ignored. The diagonal is ignored.
#' @param start the starting configuration: \code{"cmdscale"} (default, the classical
#'   scaling solution; not allowed with missing dissimilarities or zero weights),
#'   \code{"random"} (independent normal coordinates, scaled so that the expected
#'   squared distance equals the mean squared dissimilarity), an N by p matrix, or an
#'   N by p by R array of R starting configurations.
#' @param replicates the number of times the scaling is repeated, each time with a
#'   new starting configuration; the solution with the lowest criterion is returned.
#'   Default 1. With \code{start = "cmdscale"} all the replicates start from the same
#'   configuration; with an array it is the third dimension of the array.
#' @param ties the treatment of tied dissimilarities in non-metric scaling:
#'   \code{"primary"} (default, as in MATLAB: tied dissimilarities may get different
#'   disparities) or \code{"secondary"} (tied dissimilarities get the same
#'   disparity). Ignored by the metric criteria.
#' @param control a list of options of the algorithm, as returned by
#'   \code{mdscalecontrol}: maximum number of iterations, tolerances and display.
#' @param \dots arguments passed to \code{mdscalecontrol}.
#'
#' @details Let d_ij be the Euclidean distances between the points of the
#' configuration, D_ij the dissimilarities, w_ij the weights, and dh_ij the
#' disparities, i.e. the weighted least squares monotone (non-decreasing) regression
#' of the distances on the order of the dissimilarities. All sums are over the pairs
#' i < j with non-missing dissimilarity and positive weight. The criteria are
#' \tabular{ll}{
#' stress \tab sqrt( sum w (dh - d)^2 / sum w d^2 )\cr
#' sstress \tab sqrt( sum w (dh^2 - d^2)^2 / sum w d^4 )\cr
#' metricstress \tab sqrt( sum w (D - d)^2 / sum w D^2 )\cr
#' metricsstress \tab sqrt( sum w (D^2 - d^2)^2 / sum w D^4 )\cr
#' sammon \tab sum w (D - d)^2 / D / sum w D\cr
#' strain \tab sqrt( sum over all i, j of v_ij (b_ij - y_i'y_j)^2 )
#' }
#' where, for \code{"strain"}, b_ij = -D_ij^2/2, v_ij = w_ij with v_ii = 1, and the
#' points y_i are centred at their mean weighted by the row sums of the weights.
#' Without weights the minimum of the strain is the classical scaling solution
#' (\code{cmdscale}). For \code{"sstress"} the disparities are, as in MATLAB, the
#' monotone regression of the distances (not of the squared distances).
#'
#' The criterion is minimised as a function of the configuration only (for the
#' non-metric criteria the disparities are recomputed at every evaluation) by a
#' limited-memory BFGS method with a line search for the strong Wolfe conditions,
#' written in C++. The iterations stop when the criterion is smaller than
#' \code{tolfun}, or its relative change over the last \code{window} iterations is
#' smaller than \code{tolfun}, or the norm of the gradient is smaller than
#' \code{tolfun}, or the norm of the step is smaller than \code{tolx}, or after
#' \code{maxiter} iterations (with a warning). These are the stopping rules and the
#' default values of MATLAB, except that MATLAB measures the relative change over one
#' iteration (\code{window = 1}); the default \code{window = 10} avoids stopping at an
#' iteration that happens to make little progress. The algorithm is not the one of
#' MATLAB, so the two functions can stop at different local minima; from the same
#' starting configuration and with small tolerances they usually reach the same
#' solution. Since the criteria can have several local minima, \code{replicates}
#' with \code{start = "random"} is recommended.
#'
#' The returned configuration is centred at the origin and rotated to its principal
#' axes, with the sign of each axis chosen so that its largest absolute coordinate is
#' positive. For the non-metric criteria, which do not depend on the scale of the
#' configuration, the configuration is scaled so that the weighted sum of squared
#' distances equals the weighted sum of squared dissimilarities.
#'
#' @return An object of the class \code{mdscale}, a list with:
#' \tabular{ll}{
#' conf \tab the N by p configuration\cr
#' stress \tab the value of the criterion at the configuration\cr
#' disparities \tab the N by N matrix of disparities (non-metric criteria) or the
#'   dissimilarities (metric criteria); \code{NA} for the ignored pairs\cr
#' dissimilarities \tab the N by N matrix of dissimilarities\cr
#' distances \tab the N by N matrix of distances between the points of \code{conf}\cr
#' weights \tab the N by N matrix of weights (\code{NULL} if not given)\cr
#' criterion \tab the criterion\cr
#' ties \tab the treatment of ties (non-metric criteria)\cr
#' iterations \tab the number of iterations of the best replicate\cr
#' converged \tab \code{TRUE} if a stopping rule other than the iteration limit was met\cr
#' message \tab the stopping rule that ended the iterations\cr
#' allstress \tab the value of the criterion for each replicate\cr
#' call \tab the matched call
#' }
#'
#' @references Kruskal, J.B. (1964). Nonmetric multidimensional scaling: a numerical
#' method. Psychometrika, 29, 115-129.
#'
#' Takane, Y., Young, F.W., and de Leeuw, J. (1977). Nonmetric individual differences
#' multidimensional scaling: an alternating least squares method with optimal scaling
#' features. Psychometrika, 42, 7-67.
#'
#' Sammon, J.W. (1969). A nonlinear mapping for data structure analysis. IEEE
#' Transactions on Computers, 18, 401-409.
#'
#' Torgerson, W.S. (1958). Theory and Methods of Scaling. Wiley, New York.
#'
#' Borg, I., and Groenen, P.J.F. (2005). Modern Multidimensional Scaling: Theory and
#' Applications (2nd ed.). Springer, New York.
#'
#' Nocedal, J., and Wright, S.J. (2006). Numerical Optimization (2nd ed.). Springer,
#' New York.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}, Luca Coraggio \email{luca.coraggio@unina.it}
#'
#' @seealso \code{\link{mdscalecontrol}}, \code{\link{plot.mdscale}},
#'   \code{\link{kunfolding}}, \code{\link[stats]{cmdscale}}
#'
#' @examples
#' data(Univranks)
#' # Kemeny distances between the rankings of the judges
#' DK <- ConsRank::kemenyd(Univranks$rankings)
#' # non-metric scaling
#' out <- mdscale(DK, p = 2)
#' out
#' plot(out)
#' plot(out, type = "shepard")
#' # metric scaling
#' out2 <- mdscale(DK, criterion = "metricstress")
#' out2$stress
#'
#' @keywords Multidimensional scaling
#'
#' @useDynLib ConsRankClass, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats cmdscale rnorm
#' @export

mdscale <- function(D, p = 2,
                    criterion = c("sstress", "stress", "metricstress",
                                  "metricsstress", "sammon", "strain"),
                    weights = NULL, start = "cmdscale", replicates = 1,
                    ties = c("primary", "secondary"),
                    control = mdscalecontrol(...), ...) {

  cl <- match.call()                                   # the call, returned in the output
  criterion <- match.arg(criterion)                    # check the criterion
  ties <- match.arg(ties)                              # check the treatment of ties
  allcrit <- c("stress", "sstress", "metricstress",    # the criteria, in the order
               "metricsstress", "sammon", "strain")    # used by the C++ engine
  crit <- match(criterion, allcrit)                    # code of the criterion (1 to 6)
  nonmetric <- crit <= 2                               # stress and sstress are non-metric

  # ---- dissimilarities: matrix, dist object or vector of the upper triangle
  D <- mdsasmatrix(D, "D")                             # always an N by N matrix
  labs <- rownames(D)                                  # labels of the points, if any
  if (is.null(labs)) labs <- colnames(D)               # otherwise the column names
  N <- nrow(D)                                         # number of points
  if (N < 3)                                           # at least three points are needed
    stop("D must have at least 3 rows and columns: give the dissimilarities between at least 3 objects.")
  if (any(is.infinite(D)))                             # Inf is not accepted
    stop("D contains infinite values: give finite dissimilarities, and use NA for the missing ones.")
  offd <- row(D) != col(D)                             # off-diagonal positions
  if (all(!is.na(diag(D))) && all(diag(D) == 1) &&     # ones on the diagonal and
      all(D[offd] < 1, na.rm = TRUE)) {                # all other entries below one:
    D <- sqrt(1 - D)                                   # similarities, transformed as in MATLAB
  }
  if (any(D[offd] < 0, na.rm = TRUE))                  # dissimilarities must be non-negative
    stop("D contains negative values: dissimilarities must be non-negative.")
  if (!isTRUE(all.equal(D, t(D), check.attributes = FALSE))) # symmetry (NA in the same places)
    stop("D must be symmetric: give a symmetric matrix of dissimilarities.")
  diag(D) <- 0                                         # the diagonal is not used

  # ---- weights
  if (is.null(weights)) {                              # no weights:
    W <- matrix(1, N, N)                               # every pair has weight one
  } else {
    W <- mdsasmatrix(weights, "weights")               # same forms as D
    if (!all(dim(W) == c(N, N)))                       # same size as D
      stop("weights must have the same size as D.")
    if (anyNA(W) || any(W < 0))                        # non-negative, no missing values
      stop("weights must be non-negative numbers without missing values.")
    if (!isTRUE(all.equal(W, t(W), check.attributes = FALSE))) # symmetry
      stop("weights must be symmetric.")
  }
  diag(W) <- 0                                         # the diagonal is not used

  # ---- active pairs: upper triangle, non-missing dissimilarity, positive weight
  up <- which(upper.tri(D), arr.ind = TRUE)            # all pairs i < j
  dv <- D[up]                                          # their dissimilarities
  wv <- W[up]                                          # their weights
  act <- !is.na(dv) & wv > 0                           # pairs used by the criterion
  if (sum(act) < 1)                                    # nothing to fit
    stop("There are no pairs with a non-missing dissimilarity and a positive weight.")
  if (crit == 5 && any(dv[act] <= 0))                  # Sammon divides by the dissimilarities
    stop("Sammon's criterion requires strictly positive dissimilarities: use another criterion or remove the zero dissimilarities.")

  # ---- starting configurations
  if (is.character(start)) {                           # "cmdscale" or "random"
    start <- match.arg(start, c("cmdscale", "random")) # check the name
    if (length(p) != 1 || p < 1 || p != round(p) || p >= N) # p must be a valid dimension
      stop("p must be a positive integer smaller than the number of rows of D.")
    if (length(replicates) != 1 || replicates < 1 || replicates != round(replicates))
      stop("replicates must be a positive integer.")
    if (start == "cmdscale") {                         # classical scaling start
      if (any(!act))                                   # it needs all the dissimilarities
        stop("start = \"cmdscale\" cannot be used with missing dissimilarities or zero weights: use start = \"random\" or give a starting configuration.")
      Y0 <- suppressWarnings(cmdscale(D, k = p))       # classical scaling (it may return fewer columns)
      Y0 <- as.matrix(Y0)                              # make sure it is a matrix
      if (ncol(Y0) < p)                                # fewer positive eigenvalues than p:
        Y0 <- cbind(Y0, matrix(0, N, p - ncol(Y0)))    # complete with zero columns
      starts <- rep(list(Y0), replicates)              # the same start for every replicate
    } else {                                           # random start
      sc <- sqrt(mean(dv[act]^2) / (2 * p))            # E(d^2) = 2 p sc^2 = mean squared dissimilarity
      starts <- lapply(seq_len(replicates),            # one random configuration per replicate
                       function(r) matrix(rnorm(N * p), N, p) * sc)
    }
  } else if (length(dim(start)) == 3) {                # array of starting configurations
    if (dim(start)[1] != N)                            # one row per point
      stop("start must have as many rows as D.")
    p <- dim(start)[2]                                 # dimensions from the array
    starts <- lapply(seq_len(dim(start)[3]), function(r) start[, , r, drop = TRUE]) # one per replicate
    starts <- lapply(starts, function(s) matrix(s, N, p)) # keep them as matrices
  } else {                                             # a single starting matrix
    start <- as.matrix(start)                          # make sure it is a matrix
    if (nrow(start) != N)                              # one row per point
      stop("start must have as many rows as D.")
    p <- ncol(start)                                   # dimensions from the matrix
    starts <- rep(list(start), replicates)             # the same start for every replicate
  }
  if (any(sapply(starts, function(s) any(!is.finite(s))))) # starting values must be finite
    stop("start contains missing or infinite values.")

  # ---- matrices of the strain criterion
  if (crit == 6) {                                     # strain only
    B2 <- -0.5 * D^2                                   # b_ij = -D_ij^2 / 2
    B2[is.na(B2)] <- 0                                 # missing: zero (its weight is zero)
    W1 <- W                                            # weights of the strain
    W1[is.na(D)] <- 0                                  # missing dissimilarities get weight zero
    rw <- rowSums(W1)                                  # row sums, used to centre the configuration
    diag(W1) <- 1                                      # the diagonal terms have weight one
  } else {                                             # other criteria: not used
    B2 <- matrix(0, 0, 0); W1 <- matrix(0, 0, 0); rw <- numeric(0)
  }

  # ---- minimisation, one replicate at a time
  disp <- match(control$display, c("off", "final", "iter")) - 1 # display code for the engine
  best <- NULL                                         # best replicate so far
  allstress <- numeric(length(starts))                 # criterion of each replicate
  for (r in seq_along(starts)) {                       # for each replicate
    out <- mdscale_engine(dv[act], wv[act],            # active dissimilarities and weights
                          up[act, 1] - 1L, up[act, 2] - 1L, # their points (0-based)
                          N, starts[[r]], crit,        # sizes, start, criterion
                          match(ties, c("primary", "secondary")), # ties
                          B2, W1, rw,                  # strain matrices
                          control$maxiter, control$tolfun, control$tolx, disp,
                          10L, control$window)          # L-BFGS memory, window of the relative change
    allstress[r] <- out$stress                         # criterion of this replicate
    if (disp >= 1)                                     # final display, as in MATLAB
      cat(out$iter, " iterations, final ", criterion, " criterion = ",
          format(out$stress, digits = 6), "\n", sep = "")
    if (is.null(best) || out$stress < best$stress) best <- out # keep the best one
  }
  if (best$reason == 0)                                # iteration limit
    warning("Iteration limit exceeded: the minimisation of the criterion did not converge. Increase maxiter in mdscalecontrol().")

  # ---- final configuration: centred, rotated to principal axes, signs fixed
  Y <- sweep(best$conf, 2, colMeans(best$conf))        # centre at the origin
  Y <- Y %*% svd(Y)$v                                  # rotate to the principal axes
  for (a in seq_len(ncol(Y)))                          # for each axis
    if (Y[which.max(abs(Y[, a])), a] < 0) Y[, a] <- -Y[, a] # largest absolute coordinate positive
  dh <- best$dhat                                      # disparities of the active pairs
  if (nonmetric) {                                     # scale-invariant criteria:
    dY <- as.matrix(stats::dist(Y))[up[act, , drop = FALSE]]  # distances of the active pairs
    sc <- sqrt(sum(wv[act] * dv[act]^2) / sum(wv[act] * dY^2)) # sum w d^2 = sum w D^2
    if (is.finite(sc) && sc > 0) { Y <- Y * sc; dh <- dh * sc } # rescale configuration and disparities
  }
  dimnames(Y) <- list(labs, paste0("D", seq_len(ncol(Y)))) # names of points and dimensions

  # ---- disparities as a matrix
  if (nonmetric) {                                     # non-metric: the fitted disparities
    DH <- matrix(NA_real_, N, N)                       # ignored pairs stay NA
    diag(DH) <- 0                                      # zero diagonal
    DH[up[act, , drop = FALSE]] <- dh                  # upper triangle
    DH[up[act, 2:1, drop = FALSE]] <- dh               # lower triangle
  } else {                                             # metric: the dissimilarities themselves,
    DH <- D                                            # as in MATLAB
  }
  dimnames(DH) <- list(labs, labs)                     # names
  dimnames(D) <- list(labs, labs)                      # names

  msgs <- c("iteration limit reached",                 # meaning of the stopping codes
            "criterion less than tolfun",
            "relative change in the criterion less than tolfun",
            "norm of the gradient less than tolfun",
            "norm of the step less than tolx")
  res <- list(conf = Y,                                # configuration
              stress = best$stress,                    # value of the criterion
              disparities = DH,                        # disparities
              dissimilarities = D,                     # dissimilarities
              distances = as.matrix(stats::dist(Y)),          # distances in the configuration
              weights = if (is.null(weights)) NULL else W, # weights
              criterion = criterion,                   # criterion
              ties = ties,                             # treatment of ties
              iterations = best$iter,                  # iterations of the best replicate
              converged = best$reason != 0,            # convergence
              message = msgs[best$reason + 1],         # stopping rule
              allstress = allstress,                   # criterion of each replicate
              call = cl)                               # call
  class(res) <- c("mdscale", "ConsRankClass")          # class of the output
  res                                                  # return it
}

# Internal: turns a dist object, a vector of the upper triangle or a matrix into an
# N by N matrix (what is the name of the argument, for the error messages).
mdsasmatrix <- function(x, what) {
  if (inherits(x, "dist")) {                           # dist object:
    M <- as.matrix(x)                                  # full matrix with labels
  } else if (is.null(dim(x))) {                        # vector of the upper triangle
    K <- length(x)                                     # its length, N(N-1)/2
    N <- (1 + sqrt(1 + 8 * K)) / 2                     # the corresponding N
    if (K < 1 || N != round(N))                        # not a valid length
      stop(what, " must be a square matrix, a dist object or a vector of length N(N-1)/2.")
    M <- matrix(0, N, N)                               # empty matrix
    M[lower.tri(M)] <- x                               # same order as dist (and pdist)
    M <- M + t(M)                                      # symmetric
  } else {                                             # matrix or data frame
    M <- as.matrix(x)                                  # as a matrix
    if (nrow(M) != ncol(M))                            # must be square
      stop(what, " must be a square matrix, a dist object or a vector of length N(N-1)/2.")
  }
  if (!is.numeric(M))                                  # must be numeric
    stop(what, " must be numeric.")
  storage.mode(M) <- "double"                          # double precision
  M                                                    # return it
}
