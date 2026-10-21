## Splitting rule "mallows" of ranktree (Lee and Yu, 2010): a Mallows model
## (Kendall distance) is fitted within each node and the impurity is its
## deviance -2 loglik.
##
## Complete rankings and top-k partial rankings are allowed. A top-k ranking
## is built stage by stage: at stage j one of the M - j + 1 objects not yet
## placed is chosen, and V_j counts the objects among them that the center
## prefers to the chosen one. Under the Mallows model with Kendall distance the
## V_j are independent, P(V_j = v) = exp(-lambda v) / psi_{M-j+1}(lambda), with
## psi_m(lambda) = (1 - exp(-m lambda)) / (1 - exp(-lambda)) (Fligner and
## Verducci, 1986). The probability of a top-k ranking, i.e. the sum over all
## its completions, is then the product over the first k stages only:
##   log P = -lambda sum_{j <= k} V_j - sum_{j <= k} log psi_{M-j+1}(lambda),
## where sum_j V_j is the number of pairs ordered by the top-k ranking (a ranked
## object against an object ranked after it or unranked) on which the ranking
## disagrees with the center. For a set of judges with pair counts n_ab (judges
## preferring a to b) and c_j (judges with at least j ranked positions):
##   - the MLE of the center minimises T = sum of the disagreements, whatever
##     lambda > 0 (Kemeny consensus restricted to permutations);
##   - the MLE of lambda solves sum_j c_j E_lambda[V_j] = T;
##   - the deviance is D = 2 [lambda T + sum_j c_j log psi_{M-j+1}(lambda)].
## Everything depends on counts only, so that all the cutting points of a
## predictor are evaluated at once from cumulative sums.

#--------------------------------------

## Checks the rankings for the "mallows" rule and puts them in a standard form:
## ranked objects 1, ..., k, unranked objects (NA, or tied in the last
## position) all tied in position k + 1. Ties among ranked objects are not
## allowed. Returns the standardised rankings and k (number of stages).

topkrankings <- function(Y){

  M <- ncol(Y)
  out <- matrix(0, nrow(Y), M, dimnames=dimnames(Y))
  k <- integer(nrow(Y))
  bad <- logical(nrow(Y))

  for (i in seq_len(nrow(Y))){
    r <- Y[i,]
    unr <- is.na(r)
    if (!any(unr)){
      last <- r == max(r)
      if (sum(last) > 1) {unr <- last}
    }
    v <- r[!unr]
    if (length(v) == 0 || anyDuplicated(v) > 0) {bad[i] <- TRUE; next}
    out[i, !unr] <- rank(v)
    out[i, unr] <- length(v) + 1
    k[i] <- min(length(v), M - 1)
  }

  if (any(bad)){
    stop("model = \"mallows\" requires complete rankings or top-k partial rankings ",
         "('pick k out of M': unranked objects coded as NA or tied in the last position): ",
         sum(bad), " ranking(s) with ties among the ranked objects (first one: row ",
         which(bad)[1], ")")
  }

  return(list(Y=out, k=k))

} #end function

#--------------------------------------

## Statistics of each judge: pair indicators (as in pairsigns) and stage
## indicators K[i, j] = 1 if judge i ranked at least j positions (j = 1..M-1)

mallowsstats <- function(Ys, k){

  M <- ncol(Ys)
  PS <- pairsigns(Ys)
  K <- 1*outer(k, seq_len(M-1), ">=")

  return(list(P=PS$P, N=PS$N, K=K))

} #end function

#--------------------------------------

## log psi_m(lambda) and E_lambda[V] for the stage with m objects left,
## for a vector of lambdas (rows) and the stages m = M, ..., 2 (columns)

stagelogpsi <- function(lambda, M){

  m <- M:2
  out <- outer(lambda, m, function(l, m) log(-expm1(-m*l)) - log(-expm1(-l)))
  z <- lambda == 0
  if (any(z)) {out[z,] <- matrix(log(m), sum(z), length(m), byrow=TRUE)}

  return(out)

} #end function

stageexpect <- function(lambda, M){

  m <- M:2
  out <- outer(lambda, m, function(l, m) 1/expm1(l) - m/expm1(m*l))
  small <- lambda < 1e-6
  if (any(small)){
    out[small,] <- outer(lambda[small], m, function(l, m) (m-1)/2 - l*(m^2-1)/12)
  }

  return(out)

} #end function

#--------------------------------------

## Mallows fit (Kendall) for many sets of judges at once. P, N: pair counts
## (one row per set); C: stage counts; ctx: context of the center search.
## Returns the minimal disagreement T, the center (index in ctx$perms, or the
## ranking itself when M > 8), lambda and the deviance of each set.

mallowsfit <- function(P, N, C, ctx, lambdamax=50){

  P <- matrix(P, ncol=ncol(ctx$pairs0))
  N <- matrix(N, ncol=ncol(ctx$pairs0))
  C <- matrix(C, ncol=ctx$M-1)
  ns <- nrow(P)

  #center: minimal number of disagreements over the permutations
  if (!is.null(ctx$Kall)){
    #cost of pair q (a < b): N[q] if a is put before b, P[q] otherwise
    Tall <- ctx$Kall %*% t(P - N) + matrix(rowSums(N), nrow(ctx$Kall), ns, byrow=TRUE)
    best <- apply(Tall, 2, which.min)
    Tmin <- Tall[cbind(best, seq_len(ns))]
    centers <- ctx$perms[best, , drop=FALSE]
  } else {
    centers <- matrix(0, ns, ctx$M)
    Tmin <- numeric(ns)
    for (s in seq_len(ns)){
      lc <- kendalllocal(P[s,], N[s,], ctx)
      centers[s,] <- lc$center
      Tmin[s] <- lc$T
    }
    best <- NULL
  }

  #lambda: sum_j c_j E_lambda[V_j] = T, decreasing in lambda (bisection)
  e0 <- rowSums(C * stageexpect(rep(0, ns), ctx$M))
  emax <- rowSums(C * stageexpect(rep(lambdamax, ns), ctx$M))
  lo <- rep(0, ns)
  hi <- rep(lambdamax, ns)
  for (it in 1:60){
    mid <- (lo + hi)/2
    up <- rowSums(C * stageexpect(mid, ctx$M)) > Tmin
    lo[up] <- mid[up]
    hi[!up] <- mid[!up]
  }
  lambda <- (lo + hi)/2
  lambda[Tmin >= e0] <- 0
  lambda[Tmin <= emax] <- lambdamax

  dev <- 2*(lambda*Tmin + rowSums(C * stagelogpsi(lambda, ctx$M)))

  return(list(T=Tmin, best=best, centers=centers, lambda=lambda, deviance=dev))

} #end function

#--------------------------------------

## Context of the center search: exhaustive for M <= 8 (all permutations,
## Kall[k, q] = 1 if permutation k puts object pairs[q, 1] after pairs[q, 2]),
## local search for M > 8 in the split search (the center of the nodes of the
## tree is then computed exactly with the branch-and-bound of ConsRank)

mallowsctx <- function(M){

  pairs <- which(upper.tri(diag(M)), arr.ind=TRUE)
  ctx <- list(M=M, pairs=pairs, pairs0=matrix(0, 1, nrow(pairs)))
  if (M <= 8){
    ctx$perms <- .allperms(M)
    ctx$Kall <- 1*(ctx$perms[, pairs[,1], drop=FALSE] > ctx$perms[, pairs[,2], drop=FALSE])
  }

  return(ctx)

} #end function

#--------------------------------------

## Local search of the center from the pair counts (M > 8): start from the
## Copeland-like score, then move single objects to the best position until no
## move decreases the number of disagreements

kendalllocal <- function(Pv, Nv, ctx){

  M <- ctx$M
  W <- matrix(0, M, M) #W[a, b]: judges preferring a to b
  W[ctx$pairs] <- Pv
  W[ctx$pairs[,2:1]] <- Nv
  D <- W - t(W)

  o <- order(-rowSums(D))
  repeat{
    improved <- FALSE
    for (i in 1:M){
      x <- o[i]
      rest <- o[-i]
      #cost change of inserting x in each position of rest (relative to first)
      dd <- c(0, cumsum(D[x, rest]))
      j <- which.min(dd)
      if (dd[j] < dd[i] - 1e-9){
        o <- append(rest, x, after=j-1)
        improved <- TRUE
      }
    }
    if (!improved) {break}
  }

  center <- order(o)
  Tval <- sum(W[cbind(o[col(diag(M))[upper.tri(diag(M))]], o[row(diag(M))[upper.tri(diag(M))]])])

  return(list(center=center, T=Tval))

} #end function

#--------------------------------------

## Exact Mallows fit of a node (rows id of the standardised rankings Ys):
## center (all the optimal centers when M <= 8), lambda, deviance, and the
## tau_X and Kemeny error of the center as for the other rules

mallowsnode <- function(Ys, ST, id, ctx, lambdamax){

  M <- ncol(Ys)
  P <- colSums(ST$P[id, , drop=FALSE])
  N <- colSums(ST$N[id, , drop=FALSE])
  C <- colSums(ST$K[id, , drop=FALSE])

  fit <- mallowsfit(P, N, C, ctx, lambdamax)

  if (!is.null(ctx$Kall)){
    Tall <- drop(ctx$Kall %*% (P - N)) + sum(N)
    mc <- ctx$perms[which(Tall <= fit$T + 1e-9), , drop=FALSE]
  } else {
    #exact center: Kemeny consensus restricted to permutations (pairs tied in
    #the data cost the same for every permutation)
    CR <- consrank(Ys[id, , drop=FALSE], algorithm="BB", full=TRUE, ps=FALSE)
    mc <- matrix(CR$Consensus, ncol=M)
    if (!is.null(mc) && nrow(mc) > 0){
      fit <- mallowsfit(P, N, C, ctx, lambdamax) #T of the local search
      W <- matrix(0, M, M)
      W[ctx$pairs] <- P
      W[ctx$pairs[,2:1]] <- N
      Texact <- sum(W[outer(mc[1,], mc[1,], ">")])
      if (Texact <= fit$T){
        fit <- mallowsfitT(Texact, C, M, lambdamax)
      } else {
        mc <- fit$centers
      }
    } else {
      mc <- fit$centers
    }
  }
  colnames(mc) <- colnames(Ys)
  center <- mc[1,]

  Yn <- Ys[id, , drop=FALSE]
  tau <- mean(tau_x(Yn, center))
  error <- sum(kemenyd(Yn, center))/(length(id)*M*(M-1))

  return(list(center=center, mcenter=mc, lambda=fit$lambda, deviance=fit$deviance,
              tau=tau, error=error))

} #end function

## lambda and deviance for a given T

mallowsfitT <- function(Tval, C, M, lambdamax){

  C <- matrix(C, nrow=1)
  f <- function(l) sum(C * stageexpect(l, M)) - Tval
  if (Tval >= sum(C * stageexpect(0, M))) {
    lambda <- 0
  } else if (Tval <= sum(C * stageexpect(lambdamax, M))) {
    lambda <- lambdamax
  } else {
    lambda <- stats::uniroot(f, c(0, lambdamax), tol=1e-12)$root
  }
  dev <- 2*(lambda*Tval + sum(C * stagelogpsi(lambda, M)))

  return(list(T=Tval, lambda=lambda, deviance=dev))

} #end function

#--------------------------------------

## Best split of a node for the "mallows" rule. Same output as ranksplit; the
## node impurity is the mean deviance D/m, and the decrease is
## (D_father - D_left - D_right)/N. Splits leaving fewer than minprop*m judges
## in a child node are not allowed.

ranksplitmallows <- function(Y,X,N,nomin,namespreds,xID,ST,ctx,control){

  npred <- ncol(X)
  m <- length(xID)
  minsize <- max(1, control$minprop*m)
  lmax <- control$lambdamax

  #statistics of the judges of the node, all in one matrix
  np <- ncol(ST$P)
  S <- cbind(ST$P[xID, , drop=FALSE], ST$N[xID, , drop=FALSE], ST$K[xID, , drop=FALSE])
  DevFather <- mallowsfit(colSums(S[, 1:np, drop=FALSE]), colSums(S[, np+(1:np), drop=FALSE]),
                          colSums(S[, -(1:(2*np)), drop=FALSE]), ctx, lmax)$deviance

  #deviance decrease of many splits (sums of S over the left/right judges)
  decrease <- function(SL, SR, mL, mR){
    fL <- mallowsfit(SL[, 1:np, drop=FALSE], SL[, np+(1:np), drop=FALSE], SL[, -(1:(2*np)), drop=FALSE], ctx, lmax)
    fR <- mallowsfit(SR[, 1:np, drop=FALSE], SR[, np+(1:np), drop=FALSE], SR[, -(1:(2*np)), drop=FALSE], ctx, lmax)
    dec <- (DevFather - fL$deviance - fR$deviance)/N
    dec[mL < minsize | mR < minsize] <- -Inf
    return(dec)
  }

  best <- list(dec=-Inf)

  for (i in 1:npred){

    x <- X[xID, i]

    if (nomin[i] == 0){ #ordinal or numerical

      h <- sort(unique(x))
      if (length(h) < 2) {next}
      g <- match(x, h)
      ok <- which(!is.na(g))
      nh <- length(h) - 1
      SL <- colcumsum(rowsum(S[ok, , drop=FALSE], g[ok], reorder=TRUE))
      mL <- cumsum(tabulate(g[ok], nbins=length(h)))[1:nh]
      mR <- length(ok) - mL
      SR <- matrix(SL[nh+1,], nh, ncol(SL), byrow=TRUE) - SL[1:nh, , drop=FALSE]
      dec <- decrease(SL[1:nh, , drop=FALSE], SR, mL, mR)
      md <- which.max(dec)
      if (dec[md] > best$dec){
        best <- list(dec=dec[md], var=i, cut=h[md],
                     left=xID[which(x <= h[md])], right=xID[which(x > h[md])])
      }

    } else { #nominal: the candidate subsets of nominalsplit

      cand <- splittingvariables(x)
      if (is(cand, "NULL")) {next}
      if (length(cand) == 1){
        leftid <- list(which(x == cand))
        rightid <- list(which(x != cand))
        vsplit <- list(cand)
      } else {
        lev <- if (is(x, "factor")) function(j) levels(x)[cand[j,]] else function(j) cand[j,]
        leftid <- lapply(1:nrow(cand), function(j) which(!is.na(match(x, lev(j)))))
        rightid <- lapply(1:nrow(cand), function(j) which(is.na(match(x, lev(j)))))
        vsplit <- lapply(1:nrow(cand), function(j) {v <- lev(j); v[!is.na(v)]})
      }
      SL <- do.call(rbind, lapply(leftid, function(l) colSums(S[l, , drop=FALSE])))
      SR <- do.call(rbind, lapply(rightid, function(r) colSums(S[r, , drop=FALSE])))
      dec <- decrease(SL, SR, lengths(leftid), lengths(rightid))
      md <- which.max(dec)
      if (dec[md] > best$dec){
        best <- list(dec=dec[md], var=i, cut=vsplit[[md]],
                     left=xID[leftid[[md]]], right=xID[rightid[[md]]])
      }

    }

  } #end for

  if (!is.finite(best$dec)){
    return(list(DecImp=0))
  }

  return(list(DecImp=best$dec, indexing=list(left=best$left, right=best$right),
              custSplit=best$cut, NodeImpurity=list(left=NA, right=NA),
              SplitVar=namespreds[best$var], idSplitVar=best$var))

} #end function

#--------------------------------------

## Log-likelihood of each judge (standardised rankings Ys, stages k) under the
## Mallows model of the leaf it falls in (centers: rankings, one row per
## judge; lambda: one value per judge)

mallowsloglik <- function(Ys, k, centers, lambda){

  M <- ncol(Ys)
  pairs <- which(upper.tri(diag(M)), arr.ind=TRUE)
  PS <- pairsigns(Ys)
  before <- centers[, pairs[,1], drop=FALSE] < centers[, pairs[,2], drop=FALSE]
  V <- rowSums(ifelse(before, PS$N, PS$P))
  K <- 1*outer(k, seq_len(M-1), ">=")
  lp <- matrix(0, length(k), M-1)
  for (l in unique(lambda)){
    w <- which(lambda == l)
    lp[w,] <- matrix(stagelogpsi(l, M), length(w), M-1, byrow=TRUE)
  }

  return(-lambda*V - rowSums(K*lp))

} #end function

#--------------------------------------

## Goodness of fit of a tree with model "mallows": deviance and log-likelihood
## (sum over the leaves), number of parameters and AIC, BIC. As in Lee and Yu
## (2010, Table 6) the parameters are one precision per leaf plus one for each
## split; the modal rankings are not counted.

mallowsgoodness <- function(Tree){

  idt <- Tree$tsynth$children$measures$idt
  dev <- sum(unlist(Tree$node$deviance[idt]))
  nleaves <- length(idt)
  npar <- 2*nleaves - 1
  n <- Tree$node$size[[1]]

  return(list(deviance=dev, loglik=-dev/2, npar=npar, AIC=dev + 2*npar, BIC=dev + npar*log(n)))

} #end function

#--------------------------------------

## Mallows models of the leaves of a tree (model "mallows"): size, modal
## ranking, lambda, phi = exp(-lambda), expected Kendall distance of a complete
## ranking from the modal ranking, log-likelihood of the rankings in the leaf

mallowsleaves <- function(Tree){

  idt <- Tree$tsynth$children$measures$idt
  M <- ncol(Tree$Y)
  lambda <- unlist(Tree$node$lambda[idt])

  return(data.frame(node=unlist(Tree$node$number[idt]), size=unlist(Tree$node$size[idt]),
                    modal.ranking=sapply(Tree$node$class[idt], paste, collapse=" > "),
                    lambda=lambda, phi=exp(-lambda),
                    exp.distance=sapply(lambda, .expectd, M=M, dist="kendall"),
                    loglik=-unlist(Tree$node$deviance[idt])/2,
                    stringsAsFactors=FALSE))

} #end function

#--------------------------------------

## Header of print and summary for a tree with model "mallows"

mallowsheader <- function(Tree){

  g <- Tree$goodness
  M <- ncol(Tree$Y)
  nleaves <- sum(unlist(Tree$node$terminal))
  npart <- sum(Tree$ktop < M - 1)

  f2 <- function(v) formatC(v, format="f", digits=2)
  distname <- switch(Tree$control$dist, kendall="Kendall", Tree$control$dist)

  cat("\nDistance-based recursive partitioning for preference rankings: \n")
  cat("Model within nodes: Mallows model with", distname, "distance", "\n")
  cat("Number of items in the rank matrix: ", M, "\n")
  cat("Number of predictors: ", ncol(Tree$X), "\n")
  cat("Number of rankings: ", nrow(Tree$Y), " (complete: ", nrow(Tree$Y) - npart,
      ", top-k partial: ", npart, ")", "\n", sep="")
  cat("\n")
  cat("Number of terminal nodes: ", nleaves, "\n")
  cat("Log-likelihood: ", f2(g$loglik), "   Deviance: ", f2(g$deviance), "\n")
  cat("Number of parameters: ", g$npar, " (precisions: ", nleaves, ", splits: ", nleaves - 1, ")", "\n", sep="")
  cat("AIC: ", f2(g$AIC), "   BIC: ", f2(g$BIC), "\n")
  cat("Variables generating splits: ",
      paste(unique(unlist(Tree$node$varsplit[!unlist(Tree$node$terminal)])), collapse=', '), "\n")

} #end function
