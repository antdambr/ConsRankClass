#' Utility function
#'
#' Utility function to use to set the control arguments of \code{ranktree}
#'
#' @param num The maximum number of observations in a node to be split: default, 10\% of the sample size
#' @param decrmin Minimum decrease in impurity (\code{model = "kemeny"}) or in deviance divided by the
#' sample size (\code{model = "mallows"}) to split a node
#' @param algorithm (\code{model = "kemeny"} only) The algorithm used to compute the median ranking. One among "BB", "quick" (default),
#' "fast", "decor" and "BOLT" (the name is matched ignoring the case). "BOLT" is the deterministic heuristic of \code{\link[ConsRank]{bolt}},
#' fast also with many objects
#' @param otp (\code{model = "kemeny"}, for "BOLT" only) If TRUE (default), the optimal tie partitioning step of BOLT is used
#' (see \code{\link[ConsRank]{bolt}})
#' @param maxjudges (\code{model = "kemeny"}, for "BOLT" only) Maximum number of distinct rankings of the judges of a node used
#' by BOLT as starting solutions. Default Inf (all of them). Not used with \code{impurity = "median"}, where the median
#' rankings are computed from the combined input matrices of the nodes
#' @param full (\code{model = "kemeny"} only) Specifies if the median ranking must be searched in the universe of rankings including all the possible ties. Default: FALSE
#' @param itermax (\code{model = "kemeny"} only) Number of iterations for "fast" and "decor" algorithms. itermax=10 is the default option.
#' @param np (\code{model = "kemeny"}, for "decor" only) the number of population individuals. np=15 is the default option.
#' @param gl (\code{model = "kemeny"}, for"decor" only) generations limit, maximum number of consecutive generations without improvement. gl=100 is the default option.
#' @param ff (\code{model = "kemeny"}, for"decor" only) the scaling rate for mutation. Must be in [0,1]. ff=0.4 is the default option.
#' @param cr (\code{model = "kemeny"}, for"decor" only) the crossover range. Must be in [0,1]. cr=0.9 is the default option.
#' @param proc (\code{model = "kemeny"}, for "BB" only) proc=TRUE allows the branch and bound algorithm to work in difficult cases, i.e. when the number of objects is larger than 15 or 25. proc=FALSE is the default option
#' @param ps (\code{model = "kemeny"} only) If PS=TRUE, on the screen some information about how many branches are processed are displayed. Default value: FALSE
#' @param prunecost (\code{model = "kemeny"} only) The cost used by the cost-complexity pruning: "impurity" (default, the within-node
#' impurity used to grow the tree, as in D'Ambrosio and Heiser, 2016) or "error" (the Kemeny distance
#' between the rankings and the median ranking of each node, i.e. the same measure used by
#' \code{validatetree} to evaluate the subtrees). With \code{model = "mallows"} the pruning always uses the
#' deviance of the leaves
#' @param model The approach used to build the tree: \code{"kemeny"} (default; the impurity of a node is
#' based on Kemeny distances, D'Ambrosio and Heiser, 2016) or \code{"mallows"} (a Mallows model is
#' fitted within each node and the impurity is its deviance, Lee and Yu, 2010). See Details of
#' \code{\link{ranktree}}
#' @param impurity (\code{model = "kemeny"} only) The impurity of a node: \code{"pairwise"} (default, the
#' sum of the Kemeny distances between all the pairs of rankings within the node, divided by its size, as in
#' D'Ambrosio and Heiser, 2016) or \code{"median"} (the sum of the Kemeny distances of the rankings of the
#' node from its median ranking, computed with \code{algorithm}). With \code{"pairwise"} the splits do not
#' depend on the median rankings; with \code{"median"} every candidate split requires the median rankings of
#' the two child nodes. See Details of \code{\link{ranktree}}
#' @param maxcuts (\code{model = "kemeny"} with \code{impurity = "median"} only) The maximum number of cutting points
#' considered for an ordered predictor (numeric or ordered factor) in a node. Default \code{NULL}: all the cutting points
#' (the exact search). With an integer K, when a predictor has more than K cutting points in a node only K of them are
#' considered, those that leave in the left child the proportions 1/(K+1), 2/(K+1), ..., K/(K+1) of the judges of the node
#' (or the nearest ones): the search is much faster with continuous predictors, at the price of possibly missing the
#' best split. Nominal predictors are not affected
#' @param nbuckets (\code{model = "kemeny"} only) If given, the consensus of each node is the median bucket order with
#' exactly \code{nbuckets} buckets (computed by \code{\link[ConsRank]{mcbo}} with \code{algorithm}, which must be "BB",
#' "quick", "decor" or "BOLT"). Default \code{NULL}: no constraint
#' @param bucketsizes (\code{model = "kemeny"} only) If given, the consensus of each node is the median bucket order with
#' this pattern of bucket sizes, from the first (top) bucket to the last: each element is a prescribed size or \code{NA}
#' (free size), as in \code{\link[ConsRank]{mcbo}}. Requires \code{algorithm = "BOLT"}. Default \code{NULL}: no constraint
#' @param topk (\code{model = "kemeny"} only) If given, the consensus of each node is the median top-k list of length
#' \code{topk} (the top-k \emph{consensus}, computed by \code{\link[ConsRank]{consranktopk}}: the \code{topk} objects in
#' positions 1, ..., k and the others tied in position k + 1), with the exact algorithm of \code{consranktopk}, or with
#' its BOLT heuristic when \code{algorithm = "BOLT"}. Not to be confused with top-k \emph{data} (see \code{model =
#' "mallows"}). Default \code{NULL}: no constraint. At most one of \code{nbuckets}, \code{bucketsizes} and \code{topk}
#' can be given, and not with \code{full = TRUE}
#' @param dist (\code{model = "mallows"} only) The distance of the Mallows model. Only \code{"kendall"}
#' is available at present
#' @param minprop (\code{model = "mallows"} only) Minimum size of each child node, as a proportion of
#' the size of the node to be split. Default (\code{NULL}): 0.1, as in Lee and Yu (2010)
#' @param lambdamax (\code{model = "mallows"} only) Upper bound for the precision parameter
#' \eqn{\lambda}{lambda} of the Mallows model (reached when all the rankings in a node are equal to its
#' modal ranking). Default: 50
#'
#' @details The options that do not concern the chosen \code{model} are ignored, with a warning if they
#' are set explicitly.
#'
#' @return A list containing all the control parameters
#'
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#'
#' @references D'Ambrosio, A., and Heiser W.J. (2016). A recursive partitioning method for the prediction
#' of preference rankings based upon Kemeny distances. Psychometrika, vol. 81 (3), pp.774-94.
#' @references Lee, P.H., and Yu, P.L.H. (2010). Distance-based tree models for ranking data.
#' Computational Statistics and Data Analysis, 54, 1672-1682.
#'
#' @seealso  \code{\link{ranktree}}
#'
#' @export


ranktreecontrol <- function(num=NULL, decrmin=0.01, algorithm="quick",full=FALSE,
                             itermax=10,np=15,gl=100,ff=0.4,cr=0.9,proc=FALSE,
                             ps=FALSE, prunecost=c("impurity","error"),
                             model=c("kemeny","mallows"), impurity="pairwise",
                             dist="kendall", minprop=NULL, lambdamax=50,
                             otp=TRUE, maxjudges=Inf, maxcuts=NULL,
                             nbuckets=NULL, bucketsizes=NULL, topk=NULL)

{

  setbyuser <- names(match.call())[-1]
  #algorithm of the median ranking, matched ignoring the case as in consrank
  algs <- c("BB","quick","fast","decor","BOLT")
  ialg <- match(tolower(algorithm), tolower(algs))
  if (length(algorithm) != 1 || is.na(ialg)){
    stop("algorithm must be one of \"BB\", \"quick\", \"fast\", \"decor\" or \"BOLT\".")
  }
  algorithm <- algs[ialg] #name as written in consrank
  if (!is.logical(otp) || length(otp) != 1 || is.na(otp)){
    stop("otp must be TRUE or FALSE.")
  }
  if (!is.numeric(maxjudges) || length(maxjudges) != 1 || is.na(maxjudges) || maxjudges < 1){
    stop("maxjudges must be a number not smaller than 1 (Inf to use all the distinct rankings).")
  }
  prunecost <- match.arg(prunecost)
  model <- match.arg(model)
  impurity <- match.arg(impurity, c("pairwise","median"))
  dist <- .matchdist(dist)
  if (model=="mallows" && dist!="kendall"){
    stop("model = \"mallows\" is available for the Kendall distance only")
  }
  if (!is(minprop,"NULL") && (!is.numeric(minprop) || minprop < 0 || minprop >= 0.5)){
    stop("minprop must be a number in [0, 0.5)")
  }
  if (!is(maxcuts,"NULL") && (!is.numeric(maxcuts) || length(maxcuts) != 1 || is.na(maxcuts) ||
                               maxcuts < 1 || maxcuts != round(maxcuts))){
    stop("maxcuts must be NULL (all the cutting points) or a positive integer.")
  }

  #constraints on the consensus of the nodes: at most one of them
  isint <- function(v) is.numeric(v) && all(is.finite(v)) && all(v == round(v)) #whole numbers
  nconstr <- (!is(nbuckets,"NULL")) + (!is(bucketsizes,"NULL")) + (!is(topk,"NULL"))
  if (nconstr > 1){
    stop("give at most one of nbuckets, bucketsizes and topk: they are alternative constraints on the consensus of the nodes.")
  }
  if (!is(nbuckets,"NULL") && (length(nbuckets) != 1 || !isint(nbuckets) || nbuckets < 2)){
    stop("nbuckets must be a whole number not smaller than 2: the number of buckets of the consensus of the nodes.")
  }
  if (!is(bucketsizes,"NULL")){
    fixed <- bucketsizes[!is.na(bucketsizes)] #prescribed sizes (NA: free size)
    if (!is.numeric(bucketsizes) || length(bucketsizes) < 2 || (length(fixed) > 0 && (!isint(fixed) || any(fixed < 1)))){
      stop("bucketsizes must be a vector of at least two elements, each a positive whole number (prescribed size) or NA (free size), as in ConsRank::mcbo.")
    }
    if (algorithm != "BOLT"){
      stop("bucketsizes is available with algorithm = \"BOLT\" only (as in ConsRank::mcbo).")
    }
  }
  if (!is(topk,"NULL") && (length(topk) != 1 || !isint(topk) || topk < 1)){
    stop("topk must be a positive whole number: the length of the top-k consensus of the nodes.")
  }
  if (nconstr == 1 && isTRUE(full)){
    stop("full = TRUE cannot be used with nbuckets, bucketsizes or topk, which constrain the consensus to have ties.")
  }
  if ((!is(nbuckets,"NULL")) && algorithm == "fast"){
    stop("with nbuckets the algorithm must be \"BB\", \"quick\", \"decor\" or \"BOLT\" (the algorithms of ConsRank::mcbo).")
  }

  #options that do not concern the chosen model
  kemenyonly <- c("algorithm","full","itermax","np","gl","ff","cr","proc","ps","prunecost","impurity",
                  "otp","maxjudges","maxcuts","nbuckets","bucketsizes","topk")
  mallowsonly <- c("dist","minprop","lambdamax")
  ignored <- intersect(setbyuser, if (model=="mallows") kemenyonly else mallowsonly)
  if (model=="kemeny" && impurity=="pairwise" && "maxcuts" %in% setbyuser){
    ignored <- c(ignored, "maxcuts") #used only by the median impurity
  }
  if (length(ignored) > 0){
    warning("option(s) ignored with model = \"", model, "\"",
            if (model=="kemeny" && "maxcuts" %in% ignored) " and impurity = \"pairwise\"" else "",
            ": ", paste(ignored, collapse=", "), call.=FALSE)
  }

  return(list(num=num,decrmin=decrmin, algorithm=algorithm, full=full, itermax=itermax, np=np, gl=gl, ff=ff,
              cr=cr, proc=proc, ps=ps, prunecost=prunecost, model=model, impurity=impurity, dist=dist,
              minprop=minprop, lambdamax=lambdamax, otp=otp, maxjudges=maxjudges, maxcuts=maxcuts,
              nbuckets=nbuckets, bucketsizes=bucketsizes, topk=topk))
}
