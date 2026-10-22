#'Validation of the tree for preference rankings
#'
#'Validation of the tree either with a test set procedure or with v-fold cross validation
#'
#' @param Tree An object of the class "ranktree" coming form the function \code{ranktree}
#' @param testX The data frame containing the test set (predictors)
#' @param testY The matrix obtaining the test set (response)
#' @param method One between "test" (default) or "cv"
#' @param V The cross-validation parameter. Default V=5
#' @param plotting With the default option plotting=TRUE, the pruning sequence plot is visualized
#' 
#' 
#' @return A list containing:
#' \tabular{lll}{
#' tau \tab \tab the Tau_x rank correlation coefficient of the sequence of the trees\cr
#' error \tab  \tab the error of the sequence of the trees\cr
#' termnodes\tab \tab the number of terminal nodes of the sequence of the trees\cr
#' best_tau \tab \tab the best tree in terms of Tau_x rank correlation coefficient\cr
#' best_error \tab \tab the best tree in terms of error (it is the same)\cr
#' validation \tab \tab information about the validation procedure}
#' If the tree was built with \code{model = "mallows"} (see \code{\link{ranktree}}) the subtrees are
#' evaluated by the deviance of the Mallows models of their leaves, and the list contains:
#' \tabular{lll}{
#' termnodes\tab \tab the number of terminal nodes of the sequence of the trees\cr
#' deviance \tab \tab the mean deviance per ranking, \eqn{-2} times the mean log-likelihood of the test
#' (or cross-validated) rankings under the Mallows models of the leaves, of the sequence of the trees\cr
#' devse \tab \tab the standard error of the mean deviance\cr
#' best_deviance \tab \tab the best tree in terms of deviance\cr
#' best_deviance1se \tab \tab the smallest tree whose deviance is within one standard error of the
#' minimum (1-SE rule, as in Lee and Yu, 2010)\cr
#' validation \tab \tab information about the validation procedure}
#'
#' @details With \code{model = "mallows"} the test rankings can be top-k partial rankings as for
#' \code{\link{ranktree}} (unranked objects coded as \code{NA} or tied in the last position).
#' With \code{model = "kemeny"} the objects not ranked (\code{NA}) in the test rankings are treated as in the tree
#' (option \code{unranked} of \code{\link{ranktreecontrol}}): with \code{"last"} they are tied after the last ranked
#' object; with \code{"unknown"} tau and error are computed on the observed pairs only.
#'
#' @references Lee, P.H., and Yu, P.L.H. (2010). Distance-based tree models for ranking data.
#' Computational Statistics and Data Analysis, 54, 1672-1682.
#'
#' 
#' @examples
#' \donttest{
#' data(EVS)
#' EVS$rankings[is.na(EVS$rankings)] <- 3
#' set.seed(654)
#' training=sample(1911,1434)
#' tree <- ranktree(EVS$rankings[training,],EVS$predictors[training,],decrmin=0.001,num=50)
#' #test set validation
#' vtreetest <- validatetree(tree,testX=EVS$predictors[-training,],EVS$rankings[-training,]) 
#' #cross-validation
#' vtreecv <- validatetree(tree,method="cv",V=10) 
#' }
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @importFrom janitor compare_df_cols
#' 
#' @export



validatetree <- function(Tree,testX=NULL,testY=NULL,method="test",V=5,plotting=TRUE){
  
  #Tree: the tree as returned by ranktree
  #testX: test sample 
  #method: one among "testset" or "cv"
  
  if (length(Tree$subtrees) < 2){stop("The tree has only the root node: there is nothing to validate")}
  nitems <- ncol(Tree$Y)
  #rule "mallows": deviance of the Mallows models of the leaves
  mallowsrule <- identical(Tree$control$model,"mallows")
  #objects not ranked in unknown order: tau and error on the observed pairs only
  unknown <- !mallowsrule && identical(Tree$control$unranked,"unknown")
  errnorm <- nitems*(nitems-1) #normalization of the error (already done with unknown)
  if (unknown) {errnorm <- 1}
  if (method=="test"){#test set 
    
    if (is(testX,"NULL")){stop("The test set must be included")}
    if (is(testY,"NULL")){stop("The test set must be included")}
    checktest <- compare_df_cols(Tree$X,testX)
    if (isFALSE(identical(checktest[[2]],checktest[[3]]))){stop("The test set in not conform to the training set")}
    if (isFALSE(identical(colnames(testY),colnames(Tree$Y)))){stop("The items in the ranking matrices must be equal and in the same order")}
    
    treeseq <- length(Tree$subtrees)
    taos <- errors <- rep(0,treeseq)  
    
    if (mallowsrule){
      #model "mallows": deviance of the test rankings under the Mallows models of the leaves
      TK <- topkrankings(as.matrix(testY))
      llik <- matrix(0,nrow(TK$Y),treeseq)
      for (j in 1:treeseq){
        llik[,j] <- treeloglik(Tree$subtrees[[j]],testX,TK$Y,TK$k)
      }
      validate <- c(list(termnodes=Tree$pruneinfo$termnodes),
                    devsummary(llik,Tree$pruneinfo$termnodes),
                    list(validation="test set"))
      if (isTRUE(plotting)) {plot.ranktree(Tree, plot.type="pruningseq",valtree=validate)}
      return(validate)
    }

    if (!unknown && anyNA(testY)) {testY <- fillunranked(as.matrix(testY))} #unranked = "last"
    for (j in 1:length(taos)){
      yfit <- predict.ranktree(Tree$subtrees[[j]],testX)$rankings
      if (unknown) {
        uf <- unknownfit(as.matrix(testY),yfit) #observed pairs only
        taos[j] <- uf$tau
        errors[j] <- uf$error
        next
      }
      taos[j] <- mean(diag(tau_x(testY,yfit)))



      errors[j] <- mean(diag(kemenyd(testY,yfit)))


    }
    
    
    
    best <- which.max(taos)
    beste <- which.min(errors)
    validation="test set"
    
    validate=list()
    validate$tau <- taos
    validate$error <- errors/errnorm #normalized error
    validate$termnodes <- Tree$pruneinfo$termnodes
    validate$best_tau <- Tree$pruneinfo$termnodes[best]
    validate$best_error <- Tree$pruneinfo$termnodes[beste]
    validate$validation <- validation



  }else if (method=="cv"){ #end test set, begin CV
    
    treeseq <- length(Tree$subtrees)
    taos <- errors <- rep(0,treeseq) 
    
    #get alpha and geometric means of the alpha boundary points
    alpha <- Tree$alpha
    avalpha <- c(sqrt(alpha[1:(length(alpha)-1)]*alpha[2:length(alpha)]),Inf)
    
    N <- nrow(Tree$X)
    
    cvid <- 1+(seq(1:N)%%V)
    
    cvid <- cvid[sample(N)]
    
    fitcv <- list( matrix(0,nrow=nrow(Tree$X),ncol=ncol(Tree$Y)) ) 
    
    fitcv <- rep(fitcv,length(alpha))
    
    if (mallowsrule){
      kcv <- topkrankings(Tree$Y)$k
      llik <- matrix(0,N,length(alpha))
    }
    
    for (j in 1:V){
      
      idtest <- which(cvid==j)
      idtrain <- which(cvid!=j)
      controlcv <- Tree$control
      #controlcv$printinfo <- FALSE
      cvtree <- ranktree(Tree$Y[idtrain,],Tree$X[idtrain,],control=controlcv)
      
      nsubtree <- findcvtree(cvtree,avalpha)
      
      level <- sort(cvtree$pruneinfo$termnodes[nsubtree])
      
      #sapply(cvtree$subtrees,function(x) unlist(x$numnodes))
      
      
      for (p in 1:length(nsubtree)){

        if (mallowsrule){
          llik[idtest,p] <- treeloglik(cvtree$subtrees[[ nsubtree[p] ]],Tree$X[idtest,],
                                       Tree$Y[idtest, , drop=FALSE],kcv[idtest])
          next
        }

        yfit <- predict.ranktree( cvtree$subtrees[[ nsubtree[p] ]], Tree$X[idtest,])

        fitcv[[p]][idtest,] <- yfit$rankings

      }


    }

    if (mallowsrule){
      #model "mallows": cross-validated deviance under the Mallows models of the leaves
      validate <- c(list(termnodes=Tree$pruneinfo$termnodes),
                    devsummary(llik,Tree$pruneinfo$termnodes),
                    list(validation="cross validation"))
      if (isTRUE(plotting)) {plot.ranktree(Tree, plot.type="pruningseq",valtree=validate)}
      return(validate)
    }

    
    if (unknown) {
      #observed pairs only
      taos <- sapply(fitcv,function(x) unknownfit(Tree$Y,x)$tau)
      errors <- sapply(fitcv,function(x) unknownfit(Tree$Y,x)$error)
    } else {
    taos <- sapply(fitcv,function(x) mean(diag(tau_x(x,Tree$Y))))

    errors <- sapply(fitcv,function(x) mean(diag(kemenyd(x,Tree$Y))))
    }
    
    best <- which.max(taos)
    beste <- which.min(errors)
    
    
    validation="cross validation"
    
    validate=list()
    validate$tau <- taos
    validate$error <- errors/errnorm #normalized error
    validate$termnodes <- Tree$pruneinfo$termnodes
    validate$best_tau <- Tree$pruneinfo$termnodes[best]
    validate$best_error <- Tree$pruneinfo$termnodes[beste]

    validate$validation <- validation
    
    
    
    
    
  }#end cv
  
  if( isTRUE(plotting)) {
    
    plot.ranktree(Tree, plot.type="pruningseq",valtree=validate)
    
  }
  
  
  return(validate)
  
}#end function

#-----------------

findcvtree <-function(Tree,alpha0){
  
  adj <- 1+100*( 2^(-52) )
  alphas <- Tree$alpha
  k <- rep(0,length(alpha0))
  
  for (j in 1:length(alpha0)){
    
    k[j] <- sum(alphas <= alpha0[j]*adj)
    
  }
  
  return(k)
  
}

#-----------------

#log-likelihood of each new judge under the Mallows model of its leaf (rule "mallows")

treeloglik <- function(Tree,newx,Ys,k){
  
  leaf <- downtree(Tree,newx)$leafid
  centers <- do.call(rbind,Tree$node$nclass[leaf])
  lambda <- unlist(Tree$node$lambda[leaf])
  
  return(mallowsloglik(Ys,k,centers,lambda))
  
}

#-----------------

#mean deviance per judge of each subtree, its standard error, the best subtree
#and the smallest subtree within one standard error of the best (1-SE rule)

devsummary <- function(llik,termnodes){
  
  dv <- -2*llik
  dev <- colMeans(dv)
  devse <- apply(dv,2,stats::sd)/sqrt(nrow(dv))
  best <- which.min(dev)
  ok <- which(dev <= dev[best] + devse[best])
  best1se <- ok[which.min(termnodes[ok])]
  
  return(list(deviance=dev, devse=devse, best_deviance=termnodes[best],
              best_deviance1se=termnodes[best1se]))
  
}
