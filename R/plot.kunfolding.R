#' Plot Kemeny equivalent augmented unfolding solution
#'
#' Plot the Kemeny equivalent augmented unfolding coming from \code{kunfolding}
#'
#' @param x An object of the class "kunfolding"
#' @param labs The labels of the items. Defalut is NULL. If not provided, a sequence "o1,...,on" is printed, with n=number of items
#' @param labsrow The labels of the individuals. Defalut is NULL. If not provided, a sequence "1,...,m" is printed, with m=number of individuals
#' @param main Main title of the plot. Default NULL 
#' @param cols Color of the individuals. It must be numeric. Default is NULL (dark gray).
#' @param cexind cex of the individuals. Default 1
#' @param cexitems cex of the items. defaul 1
#' @param pchcol pch parameter for items points. Default 15
#' @param \dots System reserved (No specific usage)
#'
#' @return the plot of unfolding solution
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @seealso \code{\link{kunfolding}},  \code{\link{augmatrix}}
#' 
#' @examples
#' data("breakfast", package="smacof")
#' unfout <- kunfolding(breakfast)
#' itemsl <- colnames(breakfast)
#' plot(unfout,labs=itemsl,cexitems=0.8)
#' 
#' @keywords Unfolding
#' 
#' @importFrom graphics plot
#' @importFrom graphics par
#' 
#' @export

#---------

plot.kunfolding <- function(x,labs=NULL,labsrow=NULL,main=NULL,cols=NULL,cexind=1,cexitems=1,pchcol=15,...){
  unfout <- x
  
  if (!is(main,"NULL")){
    title <- main
  } else {
    title <- " "
  }
  
  if(is(labs,"NULL")){
    labs <- paste("o",seq(1:nrow(unfout$colcoord)),sep="")
  }
  
  mds_sol <- unfout$mdssol
  fitR <- unfout$rowcoord
  fitC <- unfout$colcoord
  nr <- nrow(fitR)

  mins=min(mds_sol$conf)
  maxs=max(mds_sol$conf)
  lims=max(abs(c(mins,maxs)))
  if(is(cols,"NULL")){
    cols="darkgray"
  } else {
    cols=1+cols
  }
  par(pty="s")
  plot(fitR[,1],fitR[,2],type = "p", pch = 20, cex = .5, 
       xlim=c(-lims,lims), ylim=c(-lims,lims),col=cols,main=title,
       xlab="Dimension 1", ylab="Dimension 2",xaxt='n',yaxt='n', ann="FALSE")
  if(is.null(labsrow)){
    text(fitR[,1],fitR[,2],c(1:nr),col=cols,cex=cexind)
  } else {
    text(fitR[,1],fitR[,2],labsrow,col=cols,cex=cexind)
  }
  points(fitC[,1],fitC[,2],col=2, pch=pchcol)
  text(fitC[,1],fitC[,2],labs,cex=cexitems)
  par(pty="m")
  
}

