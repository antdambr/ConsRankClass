#' Print and summary methods for mdscale
#'
#' Print and summary methods for objects of the class \code{mdscale}.
#'
#' @param x,object an object of the class \code{mdscale}, as returned by \code{mdscale}.
#' @param \dots not used.
#'
#' @details \code{print} shows the criterion, its value at the solution and the
#' convergence; \code{summary} also shows the value of the criterion for each
#' replicate and the first rows of the configuration.
#'
#' @return The object, invisibly.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#'
#' @seealso \code{\link{mdscale}}, \code{\link{plot.mdscale}}
#'
#' @examples
#' data(Univranks)
#' out <- mdscale(ConsRank::kemenyd(Univranks$rankings))
#' print(out)
#' summary(out)
#'
#' @export

print.mdscale <- function(x, ...) {
  kind <- if (x$criterion %in% c("stress", "sstress")) "Non-metric" else "Metric" # kind of scaling
  cat("\n", kind, " multidimensional scaling\n\n", sep = "") # title
  cat("Criterion:", x$criterion)                       # criterion
  if (kind == "Non-metric") cat(" (", x$ties, " approach to ties)", sep = "") # treatment of ties
  cat("\n")
  cat("Number of points:", nrow(x$conf), "\n")         # points
  cat("Dimensions:", ncol(x$conf), "\n")               # dimensions
  cat("Value of the criterion:", format(x$stress, digits = 6), "\n") # fit
  cat("Iterations:", x$iterations, "\n")               # iterations of the best replicate
  cat("Convergence:", x$message, "\n")                 # stopping rule
  if (length(x$allstress) > 1)                         # with replicates
    cat("Replicates:", length(x$allstress), "\n")      # how many
  cat("\n")
  invisible(x)                                         # return the object
}

#' @rdname print.mdscale
#' @export

summary.mdscale <- function(object, ...) {
  print(object)                                        # the short description
  if (length(object$allstress) > 1) {                  # with replicates
    cat("Criterion of each replicate:\n")              # the value of each one
    print(signif(object$allstress, 6))
    cat("\n")
  }
  cat("Configuration (first rows):\n")                 # the first rows of the configuration
  print(utils::head(object$conf, 10))
  cat("\n")
  invisible(object)                                    # return the object
}
