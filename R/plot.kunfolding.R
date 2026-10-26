#' Plot Kemeny equivalent augmented unfolding solution
#'
#' Plot the Kemeny equivalent augmented unfolding coming from \code{kunfolding}
#'
#' @param x An object of the class "kunfolding"
#' @param labs The labels of the items. Defalut is NULL. If not provided, a sequence "o1,...,on" is printed, with n=number of items
#' @param labsrow The labels of the individuals. Defalut is NULL. If not provided, a sequence "1,...,m" is printed, with m=number of individuals
#' @param main Main title of the plot. Default NULL 
#' @param cols Groups of the individuals, drawn with different colours (with \code{bw = TRUE}, different symbols): numeric codes, or a factor. Default is NULL (dark gray).
#' @param cexind cex of the individuals. Default 1
#' @param cexitems cex of the items. defaul 1
#' @param pchcol pch parameter for items points. Default 15
#' @param bw Logical. If \code{TRUE}, a black and white plot: the individuals are grey points (without
#' numbers, unless \code{labsrow} is given), the items black symbols with their labels above them. With
#' \code{cols} (numeric groups of individuals) the groups are drawn with different symbols, with a legend.
#' Default \code{FALSE} (the plot with colours).
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
#' plot(unfout,labs=itemsl,cexitems=0.8,bw=TRUE)
#'
#' @keywords Unfolding
#' 
#' @importFrom graphics plot
#' @importFrom graphics par
#' 
#' @export

#---------

plot.kunfolding <- function(x,labs=NULL,labsrow=NULL,main=NULL,cols=NULL,cexind=1,cexitems=1,pchcol=15,bw=FALSE,...){
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
  if (bw){ #black and white version (internal function below)
    plotkunfbw(fitR, fitC, lims, labs, labsrow, title, cols, cexind, cexitems, pchcol)
    return(invisible(NULL))
  }
  if(is(cols,"NULL")){
    cols="darkgray"
  } else {
    #groups given as a factor or as text: their codes 1, 2, ... (before: warning
    #"'+' not meaningful for factors" and colours NA)
    if (!is.numeric(cols)) {cols <- as.integer(factor(cols))}
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

#---------
#black and white plot of an unfolding solution (internal)

#' @importFrom graphics legend points text
plotkunfbw <- function(fitR, fitC, lims, labs, labsrow, title, cols, cexind, cexitems, pchcol){
  op <- par(pty = "s")                       #square plotting region
  on.exit(par(op))                           #restore the graphical parameters
  if (is.null(cols)){                        #no groups: one symbol for all the individuals
    pchind <- rep(20, nrow(fitR))
  } else {                                   #groups: one symbol for each group
    grp <- factor(cols)                      #groups of the individuals
    symb <- c(1, 2, 0, 5, 6, 3, 4, 8)        #open symbols, distinguishable in black and white
    if (nlevels(grp) > length(symb)){
      stop("With bw = TRUE at most 8 groups of individuals can be drawn: give fewer groups in cols.")
    }
    pchind <- symb[as.integer(grp)]          #symbol of each individual
  }
  plot(fitR[,1], fitR[,2], type = "p", pch = pchind, cex = if (is.null(cols)) .6 else .8,
       col = "gray45", xlim = c(-lims, lims), ylim = c(-lims, lims), main = title,
       xlab = "Dimension 1", ylab = "Dimension 2", xaxt = 'n', yaxt = 'n', ann = FALSE)
  if (!is.null(labsrow)){                    #labels of the individuals only if given
    text(fitR[,1], fitR[,2], labsrow, col = "gray30", cex = cexind, pos = 3, offset = 0.3)
  }
  points(fitC[,1], fitC[,2], col = "black", pch = pchcol)            #the items
  text(fitC[,1], fitC[,2], labs, cex = cexitems, pos = 3, offset = 0.5) #their labels, above
  if (!is.null(cols)){                       #legend of the groups
    legend("topright", legend = levels(grp), pch = symb[seq_len(nlevels(grp))],
           col = "gray45", bty = "n", cex = 0.8)
  }
  invisible(NULL)
}

