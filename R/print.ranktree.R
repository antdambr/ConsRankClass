#'S3 methods for ranktree
#'
#'Print methods for objects of class \code{ranktree}
#'
#'@param x An object of the class "ranktree"
#'@param \dots not used
#'
#'@return print a brief summary of the prediction tree. For a tree built with \code{model = "mallows"}
#'(see \code{\link{ranktree}}) the summary concerns the Mallows models of the leaves: log-likelihood,
#'deviance, number of parameters, AIC and BIC
#'
#'@examples
#'
#'data("Univranks")
#'tree <- ranktree(Univranks$rankings,Univranks$predictors,num=50)
#'tree
#'
#'
#'
#'@export


print.ranktree <- function(x,...){

  Tree <- x

  isleaf <- unlist(Tree$node$terminal)

  if (identical(Tree$control$model,"mallows")){ #distance-based tree (Lee and Yu, 2010)
    mallowsheader(Tree)
    cat("Leaves number:",paste(unlist(Tree$node$number[isleaf]), collapse=', '))
    cat("\n")
    return(invisible(x))
  }

  obj <- "\nRecursive partitioning method for preference rankings: \n"
  cat(obj)
  cat("Number of items in the rank matrix: ",ncol(Tree$Y), "\n")
  cat("Number of predictors: ", ncol(Tree$X), "\n")
  cat("Algorithm used to compute median rankings: ", Tree$control$algorithm, "\n")
  if (identical(Tree$control$impurity, "median")) #printed only for the new impurity, so the output is unchanged otherwise
    cat("Impurity: sum of the Kemeny distances from the median ranking of the node", "\n")
  if (!is.null(Tree$control$nbuckets)) #constrained consensus (printed only when there is a constraint)
    cat("Consensus of the nodes: median bucket order with", Tree$control$nbuckets, "buckets", "\n")
  if (!is.null(Tree$control$bucketsizes))
    cat("Consensus of the nodes: median bucket order with bucket sizes",
        paste(ifelse(is.na(Tree$control$bucketsizes), "free", Tree$control$bucketsizes), collapse=", "), "\n")
  if (!is.null(Tree$control$topk))
    cat("Consensus of the nodes: median top-", Tree$control$topk, " list", "\n", sep="")
  cat("\n")

  cat("Number of terminal nodes: ",(Tree$numnodes+1)/2,"\n")
  cat("Averaged tau_X rank correlation coefficient: ", Tree$goodness$tau, "\n")
  cat("Error (normalized Kemeny distance): ", Tree$goodness$error, "\n")
  cat("Variables generating splits: ", paste(unique(unlist(Tree$node$varsplit[!isleaf])),collapse=', '), "\n")
  cat("Leaves number:",paste(unlist(Tree$node$number[isleaf]), collapse=', '))
  cat("\n")

}
