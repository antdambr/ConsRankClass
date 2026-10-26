#'Kemeny-equivalent augmented unfolding
#'
#'Kemeny-equivalent augmented unfolding.   
#'
#' @param X A n by m data matrix, in which there are n judges and m objects to be judged. Each row is a ranking of the objects which are represented by the columns. 
#' @param p the dimensionality of the solution. Default p=2
#' @param control a list of options that control details of the \code{kunfolding} model governed by the function
#' \code{mdscontrol}. The options govern the MDS model ("ordinal" or "metric"), 
#' the initial configuration ("torgerson", "random" or "user"), the transformation ("primary", "secondary", "tertiary","spline","ratio",
#' "interval","none"), the set of weights, and other related options
#' @param \dots arguments passed bypassing \code{mdscontrol}
#' 
#' 
#' @details The MDS engine is \code{smacofSym} from \code{smacof} (default) or
#' \code{\link{mdscale}} from this package, which minimises the criteria of the
#' MATLAB function \code{mdscale} (choose it with \code{engine = "mdscale"} in
#' \code{\link{mdscontrol}}).
#' The output consists in a object of the class "kunfolding". It contains:
#' \tabular{llll}{
#' rawstress \tab \tab \tab raw stress\cr
#' nrawstress\tab \tab  \tab normalized raw stress\cr
#' stress1 \tab  \tab \tab Stress-1\cr
#' rowcoord \tab \tab \tab row (individuals) coordinates\cr
#' colcord \tab \tab \tab column (items) coordinates\cr
#' dhat \tab \tab \tab dhat\cr
#' dij\tab \tab \tab configuration distance\cr
#' shepardD \tab \tab \tab DeSarbo I Index \cr
#' kendallfit \tab  \tab \tab Kendall tau_b between transformed and fitted proximities\cr
#' tauxfit \tab \tab \tab Tau_X between transformed and fitted proximities\cr
#' avgrecov \tab \tab \tab Averaged recovery measure between raw preference data and fitted proximities\cr
#' avgedpearson \tab \tab \tab Averaged Pearson correlation between raw preference data and fitted proximities\cr
#' avgspearman \tab \tab \tab Averaged Spearman rho between raw preference data and fitted proximities\cr 
#' avgkendall \tab \tab \tab Averaged Kendall taub between raw preference data and fitted proximities\cr
#' avgtaux \tab \tab \tab Averaged Tau_X between raw preference data and fitted proximities\cr
#' resume \tab \tab \tab Resume meausures\cr
#' resumerec \tab tab \tab Resume of recovery measures\cr
#' resumeaug \tab \tab \tab Resume of augmentation matrix\cr
#' kDelta \tab \tab \tab Kemeny equivalent dissimilarity matrix\cr
#' beta \tab \tab \tab beta parameter\cr
#' alpha \tab \tab \tab alpha parameter \cr
#' interactions \tab \tab \tab n x m interaction submatrix \cr
#' csi \tab \tab \tab csi parameter \cr
#' mdssol \tab \tab \tab mds solution as returned by \code{smacof} package \cr
#' n_i \tab \tab \tab number of individuals\cr
#' n_c \tab \tab \tab number of items \cr
#' tots \tab \tab \tab total \cr
#' model \tab \tab \tab mds model \cr
#' transf \tab \tab \tab transformation used\cr
#' engine \tab \tab \tab MDS engine ("smacof" or "mdscale")\cr
#' criterion \tab \tab \tab criterion minimised by \code{mdscale} (NULL with smacof)
#'} 
#'  
#' @return An object of the class kunfolding. See details for detailed information.
#' 
#' @seealso \code{augmatrix}
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @references D'Ambrosio, A., Vera, J. F., & Heiser, W. J. (2022). Avoiding degeneracies in ordinal unfolding using Kemeny-equivalent dissimilarities for two-way two-mode preference rank data. Multivariate Behavioral Research, 57(4), 679-699.
#' 
#' 
#'
#' @keywords Unfolding 
#' @keywords Kemeny-equivalent augmented unfolding
#' 
#' @importFrom smacof smacofSym
#' @importFrom stats cor
#' 
#' @examples
#' data("breakfast", package="smacof")
#' unfout <- kunfolding(breakfast)
#' itemsl <- colnames(breakfast)
#' plot(unfout,labs=itemsl)
#' 
#' @export 

kunfolding=function(X,p=2,control=mdscontrol(...),...){
  
  if (isTRUE(control$printscr)) cat("Kemeny equivalent augmented unfolding", "\n") #message only if the summary is displayed

  #Kemeny-equivalent augmented unfolding
  
  Delta <- augmatrix(X) #build Delta matrix using Kemeny-equivalent dissimilarity
  
  n_ind <- nrow(X) # number of individuals
  n_col <- ncol(X) #number of items
  tots <- n_ind+n_col #row_points + col_points
  
  
  model <- control$model
  init <- control$init
  transf <- control$transf
  w <- control$w
  minstress <- control$minstress
  itermax <- control$itermax
  nrep <- control$nrep
  minstress = control$minstress
  itermax = control$itermax
  printscr = control$printscr
  spline.degree = control$spline.degree 
  spline.intKnots = control$spline.intKnots
  relax = control$relax
  modulus = control$modulus
  
  model <- match.arg(model, c("ordinal", "metric"), several.ok = FALSE)
  init <- match.arg(init, c("torgerson", "random","user"), several.ok = FALSE)
  transf <- match.arg(transf, c("primary", "secondary","tertiary","spline","ratio","interval","none"), several.ok = FALSE)
  
  engine <- control$engine #MDS engine: "smacof" or "mdscale"
  if (is.null(engine)) engine <- "smacof" #control lists not made by mdscontrol
  engine <- match.arg(engine, c("smacof", "mdscale"), several.ok = FALSE)
  criterion <- NULL #criterion of mdscale (NULL with smacof)

  if (engine == "smacof") {
  #smacof, minimizing normalized raw stress
   control$engine <- "smacof"
   if (transf=="none" || transf=="ratio"){
      stype <- "ratio"
      sties <- "primary"
      model <- "metric"
    }
    if(model=="ordinal"){
      stype <- "ordinal"
      sties <- "primary"
    }
    if(transf=="interval"){
      stype="interval"
      sties <- "primary"
      model <- "metric"
    }
    if(transf=="spline"){
      stype="mspline"
      sties <- "primary"
    }
    if(transf=="primary" || transf=="secondary" || transf=="tertiary"){
      sties=transf #approach to ties passed to smacofSym (tertiary too)
      stype <- "ordinal" #the transformation wins over model = "metric"
      model <- "ordinal"
    }
    sinit <- init
    if(init=="user"){ #user starting configuration
      if (is.null(control$userinit)) stop("init = \"user\" requires a starting configuration: give it with userinit in mdscontrol().")
      sinit <- control$userinit
    }
    if(init=="random"){
      sinit <- matrix(runif(tots*p,min=-2,max=2),nrow=tots,ncol=p)
      
      if(is.null(w)){
        
        cs <- mean(Delta$Delta^2)/(2*p)
        
      } else {
        
        cs <- mean((Delta$Delta[w>0])^2) / (2*p)
        
      }
      #Scale random starting configuration to have an average squared distance
      #equal to the average squared dissimilarity.
      
      sinit <- sinit*sqrt(cs)
      
    } #end if random
    
    if (printscr) { #messages only if the summary is displayed
    if(transf=="primary" || transf=="secondary" || transf=="tertiary"){

      cat("smacof engine minimizing raw stress: ", model, "mds with", transf, "approach to  ties transformation", "\n")

    } else {

      cat("smacof engine minimizing raw stress: ", model, "mds with", transf, "transformation", "\n")

    }
    }
    
    
    mds_sol <- smacofSym(Delta$Delta, ndim=p, type=stype, ties=sties,init=sinit,
                   itmax=itermax, eps=minstress, weightmat=w,
                   spline.degree=spline.degree, spline.intKnots=spline.intKnots,
                   verbose = FALSE, relax = relax, modulus = modulus)

  } else {
  #mdscale engine (criteria of the MATLAB function mdscale)
    criterion <- control$criterion #criterion chosen by the user, if any
    if (transf == "tertiary") #not available in mdscale
      stop("The tertiary approach to ties is not available with engine = \"mdscale\": use transf = \"primary\" or \"secondary\", or engine = \"smacof\".")
    if (is.null(criterion)) { #otherwise implied by model and transformation
      if (transf == "ratio" || transf == "none") {
        criterion <- "metricsstress" #metric: squared stress normalized by the dissimilarities
      } else if (transf == "primary" || transf == "secondary") {
        criterion <- "sstress" #non-metric: squared stress, the default of mdscale (the transformation wins over model, as with smacof)
      } else {
        stop("With engine = \"mdscale\" use transf = \"primary\" or \"secondary\" (ordinal), or transf = \"ratio\" or \"none\" (metric), or choose the criterion in mdscontrol().")
      }
    }
    mties <- if (transf == "secondary") "secondary" else "primary" #approach to ties (non-metric criteria)
    if (init == "torgerson") { #classical scaling start
      mstart <- "cmdscale"
    } else if (init == "random") { #random start
      mstart <- "random"
    } else { #user start
      if (is.null(control$userinit)) stop("init = \"user\" requires a starting configuration: give it with userinit in mdscontrol().")
      mstart <- control$userinit
    }
    mrep <- if (is.null(nrep)) 1 else nrep #replicates
    if (printscr) cat("mdscale engine minimizing", criterion, "\n") #message only if the summary is displayed
    mds_sol <- mdscale(Delta$Delta, p = p, criterion = criterion, weights = w,
                       start = mstart, replicates = mrep, ties = mties,
                       control = mdscalecontrol(maxiter = itermax, tolfun = minstress, tolx = minstress))
  }

  #get necessary info
  if (engine == "mdscale") {
    dhat <- mds_sol$disparities #disparities (the dissimilarities for metric criteria)
    dimnames(dhat) <- NULL
  } else if(is(mds_sol$dhat,"dist")){
    dhat <- as.matrix(mds_sol$dhat) #dhat
  } else {
    dhat <- vec2square(c(mds_sol$dhat))
  }
  
  fitd=as.matrix(dist(mds_sol$conf)) #fitted distances
  
  #fitd_norm=as.matrix(mds_sol$confdist) #configuration distance
  
  fitR=mds_sol$conf[1:n_ind,] #coordinates of individuals
  
  fitC=mds_sol$conf[(n_ind+1):tots,] #coordinates of items
  
  dhat_unf=dhat[1:n_ind,(n_ind+1):tots] #dhat interactions
  
  fitd_unf=fitd[1:n_ind,(n_ind+1):tots] #fitted distances unfolding
  
  
  
  #-----Stress measures for unfolding------------
  
  #stress measures for unfolding
  RawStress_unf=sum((dhat_unf-fitd_unf)^2) #Raw Stress unfolding
  
  NormRawStress_unf=RawStress_unf/sum(dhat_unf^2) #Normalized raw Stress unfolding
  
  Stress_1unf=sqrt(NormRawStress_unf) #Stress-1 Unfolding
  
  #other measures for unfolding
  if(is.null(w)){w <- matrix(1,nrow=n_ind,ncol=n_col)}
  
  
  #----------------------------------------------
  #DAF, Tucker, VAF
  
  DAF_unf <- 1-NormRawStress_unf #Dispersion Accounted For solution unfolding
  Tucker_unf <- sqrt(1-NormRawStress_unf) #Tucker coefficient unfolding
  VAF <- cor(c(dhat_unf),c(fitd_unf))^2  #Variance Accounted For
  #------------------------------------------------------------------------
  
  #loop for computing Shepard index
  h=0
  xy <- as.vector(fitd_unf)
  for (j in 1:length(xy)){
    h <- h+sum((1* ((xy[j]-xy[setdiff(1:length(xy),j)])/(xy[j]+xy[setdiff(1:length(xy),j)]) > 0.1)))
  }
  ShepardDindex <- 2*h/(length(xy)*(length(xy)-1)) #Shepard index
  #--------------------------------------------------------------------
  #Preparing De Sarbo's index
  
  dxy <- mean(xy)
  dx <- mean(vec2square(fitd[1:n_ind,1:n_ind]))
  dy <- mean(vec2square(fitd[(n_ind+1):tots,(n_ind+1):tots]))
  
  I1 <- log((dx/dxy))
  I2 <- log((dy/dxy))
  I3 <- log((dx/dy))
  DeSarboIndex <- sum(c(I1^2, I2^2, I3^2)) #De Sarbo's Index
  
  #-------------------------------------------------------------------------
  #goodness of fit measures between transformed and fitted proximities
  
  dij=fitd_unf
  trdij=dhat_unf
  
  fitTauX=vecTaux(as.vector(dij),as.vector(trdij)) 
  fitSpearman=cor(as.vector(dij),as.vector(trdij),method="spearman",use="pairwise") 
  fitKendall=cor(as.vector(dij),as.vector(trdij),method="kendall",use="pairwise") 
  
  #----------------------------------------------------------------------------------------
  #first resume measures
  
  nc_unf=c("RawStress","NormalizedRawStress","Stress-1","DAF","Tucker","VAF","Spearman",
           "Kendall","TauX","Shepard","DeSarbo")
  resume_unf = data.frame(rbind(RawStress_unf,NormRawStress_unf,Stress_1unf, DAF_unf, Tucker_unf,VAF,fitSpearman,fitKendall,
                                fitTauX,ShepardDindex,DeSarboIndex))
  row.names(resume_unf)=nc_unf
  colnames(resume_unf)="Unfolding Measures"
  #---------------------------------------------------------------------------------
  #Recovery measures (only for unfolding)
  
  AV_recov_pref <- matrix(1,n_ind,(n_col*(n_col-1)/2))*999
  AV_pearson <- matrix(0,n_ind,1)
  AV_Spearman <- AV_pearson
  AV_Kendall <- AV_pearson
  AV_Tau_X <- AV_pearson
  #dij=fitd_unf
  x <- as.matrix(X)
  #judges who tie all the items: Pearson, Spearman and Kendall correlations are not
  #defined for them (their ranking has zero variance)
  alltie <- apply(x, 1, function(r) length(unique(r)) == 1)
  for (j in 1:n_ind){
    XX <- t(X[j,])
    DD <- t(dij[j,])
    ps <- 1
    for (k in 1:(n_col-1)){
      AV_recov <- ( (sign(XX[k]-XX[setdiff(k:length(XX),k)])>=0)*2-1 ) ==  ( (sign(DD[k]-DD[setdiff(k:length(DD),k)])>=0)*2-1 )
      ending <- (ps-1)+length(AV_recov)
      AV_recov_pref[j,ps:ending] <- t(AV_recov) #averaged recovered preferences
      ps <- ending+1
    }
    #    print(j)
    #    print(class(t(X[j,])))
    #    print(dim(t(X[j,])))
    #    print(class(t(dij[j,])))
    #    print(dim(t(dij[j,])))
    if (alltie[j]){ #correlations not defined: NA (no call to cor, which would warn)
      AV_pearson[j,1] <- AV_Spearman[j,1] <- AV_Kendall[j,1] <- NA
    } else {
    AV_pearson[j,1] <- cor(as.matrix(x[j,]),as.matrix(dij[j,])) #averaged Pearson correlation
    AV_Spearman[j,1] <- cor(as.matrix(x[j,]),as.matrix(dij[j,]),method='spearman',use="pairwise") #averaged Spearman correlation
    AV_Kendall[j,1] <- cor(as.matrix(x[j,]),as.matrix(dij[j,]),method='kendall',use="pairwise") #averaged Kendall tau correlation
    }
    AV_Tau_X[j,1] <- tau_x(t(x[j,]),t(dij[j,])) #averaged Tau_X correlation (defined also with all ties)
  }
  #the averaged correlations are computed on the judges for which they are defined
  if (any(alltie)){
    warning(sum(alltie), " judge(s) tie all the items: they are excluded from the averaged ",
            "Pearson, Spearman and Kendall correlations (they are kept in the solution ",
            "and in the averaged Tau_X).", call. = FALSE)
  }
  avpearson <- mean(AV_pearson, na.rm = TRUE)   #averaged Pearson correlation
  avspearman <- mean(AV_Spearman, na.rm = TRUE) #averaged Spearman correlation
  avkendall <- mean(AV_Kendall, na.rm = TRUE)   #averaged Kendall correlation
  #
  ncrecov <- c("Av recovered preferences","AV Pearson correlation","AV Spearman correlation",
               "AV Kendall correlation","AV Tau_X correlation")


  resume_recov = data.frame(rbind(mean(AV_recov_pref),avpearson,avspearman
                                  ,avkendall,mean(AV_Tau_X)))
  row.names(resume_recov)=ncrecov
  colnames(resume_recov)="Recovery Measures"
  
  
  
  if (printscr){
  }
  
  interactions <- Delta$Interaction
  colnames(interactions) <- colnames(X)
  
  out <- list(rawstress=RawStress_unf, #raw stress
       nrawstress=NormRawStress_unf, #normalized raw stress
       stress1=Stress_1unf, #Stress1
       rowcoord = fitR, #row coordinates
       colcoord = fitC, #items coordinates
       dhat = dhat_unf, #dhat
       dij = dij, #configuration distance
       desarboI = DeSarboIndex, #DeSarbo I Index 
       shepardD = ShepardDindex, #Shepard I index
       spearmanfit = fitSpearman, #Spearman rho between transformed and fitted proximities
       kendallfit = fitKendall, #Kendall tau_b between transformed and fitted proximities
       tauxfit = fitTauX, #Tau_X between transformed and fitted proximities
       avgrecov=mean(AV_recov_pref), #Averaged recovery measure between raw preference data and fitted proximities
       avgedpearson=avpearson, #Averaged Pearson correlation between raw preference data and fitted proximities
       avgspearman=avspearman, #Averaged Spearman rho between raw preference data and fitted proximities
       avgkendall=avkendall, #Averaged Kendall taub between raw preference data and fitted proximities
       avgtaux=mean(AV_Tau_X), #Averaged Tau_X between raw preference data and fitted proximities
       resume=resume_unf, #Resume meausures
       resumerec=resume_recov, #Resume of recovery measures
       resumeaug=Delta$res, #Resume of augmentation matrix
       kDelta=Delta$Delta, #Kemeny equivalent dissimilarity matrix
       beta = Delta$beta, #beta parameter
       alpha = Delta$alpha, #alpha parameter
       interactions = interactions, #Interactions
       csi = Delta$csi, #csi parameter
       mdssol=mds_sol, #mds solution
       n_i=n_ind, #number of individuals
       n_c=n_col, #number of items
       tots=tots, # total objects
       model=model,
       transf=transf,
       engine=engine, #MDS engine
       criterion=criterion #criterion of mdscale (NULL with smacof)
  )

  
  class(out) <- c("kunfolding","ConsRankClass")
  
  if (printscr){print(out)}
  
  
  return(out)
  
}


#-------------------------------------
TauXDistObj <- function(X){
  
  
  c <- ncol(X)
  
  pw <- gtools::combinations(c,2)
  
  D <- matrix(0,c,c)
  
  for (j in 1:nrow(pw)){
    
    D[pw[j,1],pw[j,2]] <- vecTaux(X[,pw[j,1]],X[,pw[j,2]])
    
  }
  
  D <- D+t(D)+diag(1,c);
  
  return(D)
  
}
#------------------------------------

vec2square <- function(x) {
  
  stopifnot(is(x,"numeric") || is(x,"matrix"))
  
  if (is(x,"vector")) {
    n <- length(x)
    m <- (1+sqrt(1+8*n))/2
    if(!isnatural(m) ){
      stop("The size of the vector x is not correct to be a valid distance squared symmetric matrix")
    }
    inds <- c()
    for (i in 1:(m-1)){
      inds <- c(inds, (1+i+(i-1)*m):(i*m))
    }
    y <- numeric(m*m)
    y[inds] <- x
    y <- matrix(y, m, m) + t(matrix(y, m, m))
    
  } else if (is(x,"matrix")) {
    m <- nrow(x); n <- ncol(x)
    if (m != n)
      stop("Argument 'x' must be a vector or a square matrix.")
    if (any(diag(x) != 0))
      stop("Argument 'x' can only have 0s on the diagonal.")
    if (n == 1) return(c())
    inds <- c()
    for (i in 1:(n-1))
      inds <- c(inds, (1+i+(i-1)*n):(i*n))
    y <- x[inds]
    
  } else
    stop("Argument 'x' must be a numeric vector or matrix.")
  
  return(y)
}
#-----------------
isnatural <- function(x, tol = .Machine$double.eps^0.5) {x > tol & abs(x - round(x)) < tol}







