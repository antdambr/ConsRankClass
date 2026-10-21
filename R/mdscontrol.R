#' Utility function
#'
#' Utility function to use to set the control arguments of \code{kunfolding}
#'
#' @param model Specifies MDS model. One among "ordnal (default)" or "metric"
#' @param init Initial configuration. One among "torgerson" (degault), "random" or "user"
#' @param transf The transformation. One among "primary" (default), "secondary","tertiary","spline","ratio","interval","none"
#' @param userinit The user initial configuration if "init" has been set as "user"
#' @param w The set of weigths. Default: NULL
#' @param minstress the minimum stress (for stress method). Default 1e-5
#' @param itermax Maximum number if iterations. Default 500
#' @param printscr Display the messages and the summary of the model. Default TRUE
#' @param spline.degree Degree of spline transformation. Default 2
#' @param spline.intKnots Interior knots. default 2
#' @param relax Relax the solution. Default FALSE
#' @param modulus Modulus. Default 1
#' @param engine The MDS engine: "smacof" (default, the function \code{smacofSym} of
#'   the package \code{smacof}) or "mdscale" (the function \code{\link{mdscale}} of
#'   this package, which minimises the criteria of the MATLAB function
#'   \code{mdscale}). With both engines \code{transf = "primary"},
#'   \code{"secondary"} or \code{"tertiary"} gives an ordinal MDS whatever
#'   \code{model} is. With "mdscale", \code{transf = "primary"} or
#'   \code{"secondary"} minimises the squared stress
#'   (\code{criterion = "sstress"}, the default of \code{\link{mdscale}}) with that
#'   approach to ties, and
#'   \code{transf = "ratio"} or \code{"none"} minimises the metric squared stress
#'   (\code{criterion = "metricsstress"}); the other transformations are not
#'   available. \code{minstress} and \code{itermax} are passed to
#'   \code{\link{mdscalecontrol}} as \code{tolfun} (and \code{tolx}) and
#'   \code{maxiter}.
#' @param criterion Only for \code{engine = "mdscale"}: the criterion of
#'   \code{\link{mdscale}} ("stress", "sstress", "metricstress", "metricsstress",
#'   "sammon" or "strain"). Default NULL: the criterion implied by \code{model} and
#'   \code{transf}, as described above.
#' @param nrep Only for \code{engine = "mdscale"}: the number of replicates (see the
#'   argument \code{replicates} of \code{\link{mdscale}}). Default 1.
#'
#' @return A list containing all the control parameters
#'
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#'
#' @seealso  \code{\link{kunfolding}}, \code{\link{mdscale}}
#'
#' @export


mdscontrol <- function(model = "ordinal",
                       init = "torgerson",
                       transf="primary",
                       userinit = NULL,
                       w = NULL,
                       minstress = 1e-5,
                       itermax = 500,
                       printscr = TRUE,
                       spline.degree = 2,
                       spline.intKnots = 2,
                       relax = FALSE,
                       modulus = 1,
                       engine = "smacof",
                       criterion = NULL,
                       nrep = 1)
{

  engine <- match.arg(engine, c("smacof", "mdscale")) # check the engine

  return(list(model = model, init = init, transf = transf, w = w,
              minstress = minstress, itermax = itermax,
              printscr=printscr, spline.degree=spline.degree, spline.intKnots=spline.intKnots,
              relax=relax, modulus=modulus,
              userinit = userinit, engine = engine, criterion = criterion, nrep = nrep))

}
