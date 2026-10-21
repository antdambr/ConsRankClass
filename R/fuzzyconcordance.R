#'Normalized Degree of Concordance (NDC) and Adjusted Concordance Index (ACI)
#'
#'Given two fuzzy (Ruspini) partitions, it compute the NDC and the ACI. NDC is the fuzzy version
#'of the Rand Index, as well as ACI is the fuzzy version of the Adjusted Rand Index
#'
#' @param P A fuzzy partition. It has to be a matrix with n rows and k columns. 
#' Each column is expression of the degree of membership of the i-th row over the k partitions (see details).
#' @param Q A fuzzy partition. It has to be a matrix with n rows and h columns. 
#' Each column is expression of the degree of membership of the i-th row over the h partitions (see details).  
#' @param nperms \code{NULL} (default): the expected NDC under random permutations, needed by ACI,
#' is computed exactly. A number: it is estimated with \code{nperms} random permutations,
#' as in versions <= 1.1.0 (slower, kept for reproducibility).
#'
#' @details Both P and Q, or only one of those, can be crisp (or hard) partitions. In this case, each row must contain
#' either 0 or 1, and the sum of the i-th row must be 1. In other words, either P or Q (or both) are
#' expressed in terms of dummy coding. If both partitions are crisp, then NDC is equal to Rand Index 
#' and ACI is equal to Adjusted Rand Index.
#' This function can be used to externally validate the output of any fuzzy clustering method
#'
#' ACI = (NDC - E[NDC]) / (1 - E[NDC]), where E[NDC] is the expected NDC when the rows of Q
#' are randomly permuted. Under a random permutation each pair of units of P is matched with
#' a pair of units of Q drawn uniformly, so E[NDC] = 1 - mean |e_P(p) - e_Q(q)| over all
#' couples of pairs (p, q). This mean is computed exactly in O(N log N), N = n(n-1)/2, by
#' sorting the two sets of pairwise concordances: it is the limit of the permutation
#' estimate when \code{nperms} goes to infinity.
#'
#' @return A list containing:
#' \tabular{lll}{
#' ACI \tab \tab the Adjusted Concordance Index\cr
#' NDC \tab \tab the Normalized Degree of Concordance
#' }
#' 
#' @examples
#' #two random fuzzy partitions
#' P = rbind(c(0.5259,  0.1656,    0.3085),
#' c(0.5623,    0.1036,    0.3341),
#' c(0.2508,    0.1849,    0.5643),
#' c(0.5654,    0.1934,    0.2413),
#' c(0.4529,    0.1679,    0.3792),
#' c(0.2390,    0.1758,    0.5852),
#' c(0.3114,    0.1743,    0.5143),
#' c(0.4188,    0.1392,    0.4420),
#' c(0.5830,    0.1655,    0.2514),
#' c(0.5860,    0.1171,    0.2969),
#' c(0.2630,    0.1706,    0.5664),
#' c(0.5882,    0.1032,    0.3086),
#' c(0.5829,    0.1277,    0.2894),
#' c(0.3942,    0.1046,    0.5012),
#' c(0.5201,    0.1097,    0.3702),
#' c(0.2568,    0.1823,    0.5609),
#' c(0.3687,    0.1695,    0.4618),
#' c(0.5663,    0.1317,    0.3020),
#' c(0.5169,    0.1950,    0.2881),
#' c(0.5838,    0.1034,    0.3128))
#'
#' Q = rbind(c(0.4494,    0.3755,    0.1751),
#' c(0.5219,    0.3526,    0.1255),
#' c(0.3432,    0.5062,    0.1506),
#' c(0.3120,    0.5181,    0.1699),
#' c(0.5362,    0.2747,    0.1891),
#' c(0.4082,    0.3959,    0.1959),
#' c(0.4670,    0.3782,    0.1547),
#' c(0.4276,    0.4585,    0.1139),
#' c(0.4013,    0.4837,    0.1149),
#' c(0.3724,    0.5019,    0.1258),
#' c(0.5055,    0.3104,    0.1841),
#' c(0.4027,    0.4719,    0.1254),
#' c(0.3565,    0.4620,    0.1814),
#' c(0.6106,    0.2650,    0.1244),
#' c(0.5595,    0.2476,    0.1929),
#' c(0.4657,    0.3993,    0.1350),
#' c(0.2964,    0.5839,    0.1197),
#' c(0.5387,    0.3362,    0.1251),
#' c(0.4043,    0.4341,    0.1616),
#' c(0.5631,    0.2895,    0.1473))
#'          
#' ci <- fuzzyconcordance(P,Q) 
#' 
#' #generate a random fuzzy partition with two components (clusters)
#' Q2 <- matrix(runif(20),ncol=1)
#' Q2 <- cbind(Q2,1-Q2)
#' 
#' ci2 <- fuzzyconcordance(P,Q2)
#' 
#' #generate a random crisp partition
#' P2 <- t(rmultinom(20,1,c(0.3,0.3,0.4)))
#' 
#' ci3 <- fuzzyconcordance(P2,Q)
#' #--------------------
#' \dontrun{
#' # install.packages("Rankcluster")
#' library("Rankcluster") # model-based clustering algorithm for
#'                      #  ranking data by Biernacki and Jacques (2013)
#'                      #  <doi:10.1016/j.csda.2012.08.008>
#' data(APA)
#' set.seed(136) #for reproducibility
#' rcres <- rankclust(APA$data,K=3)  # solution with 3 centers, it takes about 75 seconds
#' ##
#' ccares <- cca(APA$data,k=3) #solution with 3 components, it takes about 7 seconds
#' ##
#' ci <- fuzzyconcordance(rcres[3]@tik,ccares$pk)
#' ci$ACI  # 0.0226 means that the two partitions are similar (see NDC below), 
#'         # but their similarity is mainly due to chance
#' ci$NDC
#' }
#'  
#'   
#'     
#' @references 
#' D’Ambrosio, A., Amodio, S., Iorio, C., Pandolfo, G. and Siciliano, R. (2021).
#'  Adjusted Concordance Index: an Extension of the Adjusted Rand Index to Fuzzy Partitions. 
#'  Journal of Classification vol. 38(1), pp. 112–128 (2021). DOI: 10.1007/s00357-020-09367-0
#' 
#' Hullermeier, E., Rifqi, M., Henzgen, S., and Senge, R. (2012). Comparing fuzzy partitions: a generalization of the 
#' Rand index and related measures. IEEE Transactions on Fuzzy Systems, 20(3), 546–556. DOI: 10.1109/TFUZZ.2011.2179303
#' 
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' 
#' @seealso \code{\link{cca}}
#' 
#' @importFrom stats rmultinom runif
#' @importFrom methods as is
#' 
#' @export
#' 


fuzzyconcordance <- function(P, Q, nperms = NULL){

  if(is(P,"data.frame")){P <- as.matrix(P)}
  if(is(Q,"data.frame")){Q <- as.matrix(Q)}
  m <- nrow(P)

  if (m != nrow(Q)){
    stop("The partitions must have the same number of rows. The partitions are referred to the same statistical units")
  }

  ## concordances of all the N = m(m-1)/2 pairs of units, in 'dist' order
  ep <- e_matrix(P)
  eq <- e_matrix(Q)
  N <- m * (m - 1) / 2

  NDC <- 1 - sum(abs(ep - eq)) / N  #Hullermeier NDC

  if (is.null(nperms)) {
    ## exact expected NDC under random permutations of the units of Q:
    ## sum over all couples of pairs (p, q) of |ep[p] - eq[q]|, via sorting
    mNDC <- 1 - sumabsdiff(ep, eq) / N^2
  } else {
    ## permutation estimate, as in versions <= 1.1.0 (full m x m matrices)
    pm <- matrix(runif(m * nperms), nrow = m)
    allp <- apply(pm, 2, order)
    EP <- as.matrix(1 - 0.5 * dist(P, method = "manhattan"))
    EQ <- as.matrix(1 - 0.5 * dist(Q, method = "manhattan"))
    NDCs <- vapply(seq_len(nperms), function(col) {
      d <- allp[, col]
      1 - sum(abs(EP - EQ[d, d])) / (m * (m - 1))
    }, numeric(1))
    mNDC <- mean(NDCs)
  }
  ACI <- (NDC - mNDC) / (1 - mNDC)

  return(list(ACI = ACI, NDC = NDC))

}

#----------------------

#' @importFrom proxy dist

e_matrix <- function(x){

  # concordance of each pair of units: 1 - half the L1 distance of the memberships
  as.vector(1 - (0.5 * dist(x, method = "manhattan")))

}

#----------------------

## sum of |a[p] - b[q]| over all the length(a) * length(b) couples, O(N log N):
## with b sorted and its cumulative sums, for each a[p] the k values b <= a[p]
## give k * a[p] - sum(b[1:k]), the others sum(b[-(1:k)]) - (nb - k) * a[p].
## a is sorted too, so that findInterval() scans b in order (much faster).

sumabsdiff <- function(a, b){

  a <- sort(a, method = "radix")
  b <- sort(b, method = "radix")
  nb <- length(b)
  cb <- c(0, cumsum(b))
  k <- findInterval(a, b)
  sum(a * k - cb[k + 1]) + sum((cb[nb + 1] - cb[k + 1]) - a * (nb - k))

}
