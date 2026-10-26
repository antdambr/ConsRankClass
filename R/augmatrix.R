#'Kemeny-equivalent augmented dissimilarity matrix
#'
#'Kemeny-equivalent augmented dissimilarity matrix  
#'
#' @param X A n by m data matrix, in which there are n judges and m objects to be judged. Each row is a ranking of the objects which are represented by the columns. 
#'  
#' @return A list containing the dissimilarity matrix and othe information about the augmented matrix. See details for detailed information.
#' 
#' @details First the matrix is transformed with the tau_X rank correlation coeficient, then it is normalized. The output contains:
#' \tabular{llll}{
#' Delta \tab \tab \tab the augmented dissimilarity matrix\cr
#' Interaction \tab \tab  \tab the submatrix containnig the interactions individuals-items\cr
#' Objects \tab  \tab \tab the submatrix containing the within-items proximities\cr
#' Indiv \tab \tab \tab the submatrix containing the within-individuals proximities \cr
#' beta \tab \tab \tab the beta parameter \cr
#' alpha \tab \tab \tab the alpha parameter \cr
#' csi \tab \tab \tab the csi parameter \cr
#' res \tab \tab \tab the resume of th eaugmentation in terms of:\cr
#' \tab TauX \tab \tab tau_x rank correlation coefficient\cr
#' \tab Kendall \tab \tab kendall rank correlation coefficient\cr
#' \tab Spearman \tab \tab Spearman correlation coefficient
#'} 
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @references D'Ambrosio, A., Vera, J. F., & Heiser, W. J. (2022). Avoiding degeneracies in ordinal unfolding using Kemeny-equivalent dissimilarities for two-way two-mode preference rank data. Multivariate Behavioral Research, 57(4), 679-699.
#' 
#'
#' @keywords Unfolding 
#' @keywords Kemeny-equivalent dissimilarity
#' 
#' 
#' @export 


augmatrix <- function(X){
  # Process a rectangular (individuals x items) matrix to be analyzed by
  # standard MDS programs to perform multidimensional unfolding.
  # First the matrix is transformed with the tau_X rank correlation
  # coeficient, then it is normalized.
  
  if(is(X,"data.frame")){
    X <- as.matrix(X)
  }
  
  nobj <- ncol(X)
  nind <- nrow(X)
  maxd <- nobj*(nobj-1);
  
  
  
  Outputs <- PreprocessingTauX(X)
  Interaction <- Outputs$Interaction
  Objects <- Outputs$Objects
  Indiv <- Outputs$Indiv
  
  Outputs2 <- keqtransf(Objects,Indiv,Interaction)
  Delta <- Outputs2$Delta
  beta <- Outputs2$beta
  alpha <- Outputs2$alpha
  csi <- Outputs2$csi
  
  
  TauX <- diag(tau_x(X,Interaction))
  CorK <- diag(cor((X),(Interaction),method="kendall",use="pairwise"))
  CorS <- diag(cor((X),(Interaction),method="spearman",use="pairwise"))
  
  MTau <- mean(TauX)
  MKendall <- mean(CorK)
  MSpearman <- mean(CorS)
  
  nr <- c("TauX","Kendall","Spearman")
  nc <- c("Average","Min","Max")
  
  res <- data.frame(cbind( rbind(MTau,MKendall,MSpearman), rbind(min(TauX),min(CorK),min(CorS)),
                              rbind(max(TauX),max(CorK),max(CorS)) ))
  colnames(res) <- nc
  row.names(res) <- nr
  
  
  return(list (Delta=Delta, Interaction=Interaction, Objects=Objects,
               Indiv=Indiv, beta=beta, alpha=alpha, csi=csi, res=res) )
  
}
#-------------------------------
PreprocessingTauX <- function(X){
  
  # X is a n times c rectangular matrix,
  # with n individuals and c items
  
  c <- ncol(X)
  
  # generate matrix of centers of gravity of the objects
  obj <- matrix(2,c,c)
  
  for (k in 1:c){
    
    obj[k,k] <- 1
    
  }
  
  Interaction <- 1-tau_x(X,obj)
  Interaction <- as.matrix(Interaction)
  
  Objects <- 1-TauXDistObj(X) 
  Objects <- as.matrix(Objects)
  
  Indiv = 1-tau_x(X)
  Indiv=as.matrix(Indiv)
  
  Ds <- rbind(  cbind(Indiv,Interaction), cbind(t(Interaction),Objects) )
  
  return(list(Interaction=Interaction,Objects=Objects,Indiv=Indiv,Ds=Ds))
  
}

#-------------------------------------------------------------

keqtransf <- function(Objects,Indiv,Interaction){
  
  #Kemeny-equivalent transformation
  
  Objects <- as.matrix(Objects)
  Indiv <- as.matrix(Indiv)
  Interaction <- as.matrix(Interaction)
  
  
  nobj <- nrow(Objects)
  nind <- nrow(Indiv)
  
  beta <- sqrt(nobj^2/sum(sum(Objects^2)))
  Objects1 <- beta*Objects
  alpha <- sqrt(nind^2/sum(sum(Indiv^2)))
  Indiv1 <- alpha*Indiv
  csi <- sqrt((nind*nobj)/sum(sum(Interaction^2)))
  Interaction1 <- csi*Interaction
  
  Delta1 <- rbind(Objects1,Interaction1) #Matriz bloque columna: n+R x n
  Delta2 <- rbind(t(Interaction1),Indiv1) #Matriz bloque columna: R x n+R
  
  #individuals on the top    
  
  Delta <- rbind( cbind(Indiv1,Interaction1), cbind(t(Interaction1), Objects1) )
  
  return( list (Delta=Delta,beta=beta, alpha=alpha, csi=csi))
  
}

#------------------------------------