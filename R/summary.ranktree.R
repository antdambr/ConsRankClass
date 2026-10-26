#'S3 methods for ranktree
#'
#'Summary methods for objects of class \code{ranktree}
#'
#'@param object An object of the class "ranktree"
#'@param \dots not used
#'
#'@return it shows the summary of the prediction tree. For a tree built with \code{model = "mallows"}
#'(see \code{\link{ranktree}}) it shows, for each leaf, the modal ranking and the precision
#'\eqn{\lambda}{lambda} of the Mallows model, \eqn{\phi = e^{-\lambda}}{phi = exp(-lambda)}, the expected
#'Kendall distance of a complete ranking from the modal ranking and the log-likelihood of the rankings in
#'the leaf, and it returns (invisibly) these values as a data frame. For the other trees it returns
#'(invisibly) the paths of the leaves, as given by \code{\link{treepaths}}
#'
#'@examples
#'
#' data("Univranks")
#'tree <- ranktree(Univranks$rankings,Univranks$predictors,num=50)
#'summary(tree)
#'
#'@export
#'

summary.ranktree <- function(object,...){
  #object is of the class "ranktree"
  Tree <- object
  isleaf <- unlist(Tree$node$terminal)
  forsum <- treepaths(Tree)
  ll <- length(forsum$paths)

  if (identical(Tree$control$model,"mallows")){ #distance-based tree (Lee and Yu, 2010)
    mallowsheader(Tree)
    cat("\n")
    for (j in seq_len(ll)){
      cat(forsum$paths[[j]], "\n\n")
    }
    cat("Mallows models of the terminal nodes:\n")
    lv <- mallowsleaves(Tree)
    lv[,c("lambda","phi","exp.distance")] <- round(lv[,c("lambda","phi","exp.distance")],4)
    lv$loglik <- round(lv$loglik,2)
    print(lv, row.names=FALSE)
    return(invisible(lv))
  }

  obj <- "\nRecursive partitioning method for preference rankings: \n"
  cat(obj)
  cat("Number of terminal nodes: ",length(forsum$leaves),"\n")
  cat("Averaged tau_X rank correlation coefficient: ", Tree$goodness$tau, "\n")
  cat("Error (normalized Kemeny distance): ", Tree$goodness$error, "\n")
  cat("Variables generating splits: ", paste(unique(unlist(Tree$node$varsplit[!isleaf])),collapse=', '), "\n\n")

  for (j in seq_len(ll)){
    cat(forsum$paths[[j]], "\n\n")
  }

  #return the paths of the leaves invisibly (before, the value NULL of the loop was printed)
  return(invisible(forsum))

}
