#'S3 methods for kunfolding
#'
#'Print methods for objects of class \code{kunfolding}
#'
#'@param x An object of the class "kunfolding"
#'@param \dots not used
#'
#'@return print a brief summary of the Kemeny equivalent augmented unfolding
#'
#'
#'@export


print.kunfolding <- function(x,...){
  
  out <- x
  
  obj <- "\nKemeny-equivalent augmented unfolding: \n"
  cat(obj)
  if(identical(out$engine, "mdscale")){ #mdscale engine
    cat("mdscale engine minimizing", out$criterion, "\n\n")
  } else if(out$transf=="primary" || out$transf=="secondary" || out$transf=="tertiary"){
    cat("smacof engine minimizing raw stress: ", out$model, "mds with", out$transf, "approach to ties transformation.", "\n\n")
  } else {
    cat("smacof engine minimizing raw stress: ", out$model, "mds with", out$transf, "transformation.", "\n\n")
  }
  cat("Dimensions:", ncol(out$rowcoord), sep=" ", "\n")
  cat("Number of individuals:", nrow(out$rowcoord),sep=" ", "\n")
  cat("Number of itrems:", nrow(out$colcoord), sep=" ", "\n\n")
  
  cat("Augmented dissimilarity matrix measures:", "\n")
  print(out$resumeaug)
  cat("\n")
  
  cat("Unfolding measures:", "\n")
  cat("Raw Stress",out$resume[1,], "\n", sep=" ")
  for (j in 2:nrow(out$resume)){
    if(j==2){
      cat("Normalized Raw Stress",out$resume[j,], "\n", sep=" ")
    } else {
    cat(rownames(out$resume)[j],out$resume[j,], "\n", sep=" ")
    }
  }
  cat("\n")
  cat("Recovery measures","\n")
  cat(rownames(out$resumerec)[1],out$resumerec[1,], "\n", sep=" ")
  for (j in 2:nrow(out$resumerec)){
    cat(rownames(out$resumerec)[j],out$resumerec[j,], "\n", sep=" ")
  }
  

}
