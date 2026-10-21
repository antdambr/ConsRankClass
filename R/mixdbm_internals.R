## -----------------------------------------------------------------------------
## [IT] .distmatrix: matrice nu x G delle distanze dei ranking distinti dai
## centri delle componenti.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.distmatrix <- function(U, centers, dist) {
  D <- matrix(0, nrow(U), nrow(centers))
  for (g in seq_len(nrow(centers))) D[, g] <- .dist2center(U, centers[g, ], dist)
  D
}

## -----------------------------------------------------------------------------
## [IT] .estep: passo E. Per ogni ranking distinto calcola log p_g + log
## f_g(r) (la componente uniforme, se c'e', ha densita' 1/M!), le probabilita'
## a posteriori z (nu x Gt) e la log-verosimiglianza osservata (pesata con le
## frequenze w).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.estep <- function(U, w, par, ctx) {
  nu <- nrow(U); G <- length(par$lambda); Gt <- length(par$p)
  logf <- matrix(-lgamma(ctx$M + 1), nu, Gt)
  for (g in seq_len(G))
    logf[, g] <- -par$lambda[g] * par$D[, g] - .logpsi(par$lambda[g], ctx$M, ctx$dist)
  lp <- logf + matrix(log(par$p), nu, Gt, byrow = TRUE)
  lse <- .logsumexprows(lp)
  list(z = exp(lp - lse), loglik = sum(w * lse))
}

## -----------------------------------------------------------------------------
## [IT] .mstep: passo M dati i pesi wz[u, g] (numero atteso, o simulato nel
## SEM, di copie del ranking distinto u assegnate alla componente g):
## proporzioni = medie dei pesi; centri con .centerupdate; precisioni
## risolvendo E_lambda[d] = distanza media pesata (vincolo "free": una per
## componente; "equal": una sola, con la distanza media di tutte le
## componenti non uniformi).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.mstep <- function(U, wz, par, lambdamodel, ctx) {
  G <- length(par$lambda)
  par$p <- colSums(wz) / sum(wz)
  for (g in seq_len(G)) {
    if (sum(wz[, g]) < 1e-10) next
    new <- .centerupdate(U, wz[, g], ctx$dist, par$centers[g, ], ctx)
    if (any(new != par$centers[g, ])) {
      par$centers[g, ] <- new
      par$D[, g] <- .dist2center(U, new, ctx$dist)
    }
  }
  if (G > 0) {
    num <- colSums(wz[, seq_len(G), drop = FALSE] * par$D)
    den <- colSums(wz[, seq_len(G), drop = FALSE])
    if (lambdamodel == "free") {
      for (g in seq_len(G)) if (den[g] > 1e-10)
        par$lambda[g] <- .solvelambda(num[g] / den[g], ctx$M, ctx$dist, ctx$lambdamax)
    } else {
      par$lambda[] <- .solvelambda(sum(num) / sum(den), ctx$M, ctx$dist, ctx$lambdamax)
    }
  }
  par
}

## -----------------------------------------------------------------------------
## [IT] .emone: algoritmo EM da un singolo punto di partenza. Si ferma quando
## l'aumento relativo della log-verosimiglianza e' sotto tol o dopo maxiter
## iterazioni. Restituisce parametri, posteriori, log-verosimiglianza e la
## sua traiettoria (che deve essere non decrescente).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.emone <- function(U, w, lambdamodel, ctx, init, maxiter, tol) {
  par <- init
  par$D <- .distmatrix(U, par$centers, ctx$dist)
  llold <- -Inf; lltrace <- numeric(0)
  for (it in seq_len(maxiter)) {
    e <- .estep(U, w, par, ctx)
    lltrace <- c(lltrace, e$loglik)
    if (abs(e$loglik - llold) < tol * (1 + abs(e$loglik))) break
    llold <- e$loglik
    par <- .mstep(U, w * e$z, par, lambdamodel, ctx)
  }
  list(centers = par$centers, lambda = par$lambda, p = par$p, z = e$z,
       loglik = e$loglik, iter = it, lltrace = lltrace,
       converged = it < maxiter)
}

## -----------------------------------------------------------------------------
## [IT] .semone: EM stocastico (Celeux & Diebolt, 1985), usato solo per
## generare buone partenze. Dopo il passo E ogni osservazione viene assegnata
## a una componente estratta dalle sue probabilita' a posteriori (passo S) e
## il passo M usa questa partizione. La catena puo' uscire dai massimi
## locali; si restituiscono i parametri visitati con la log-verosimiglianza
## piu' alta (da raffinare poi con l'EM).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.semone <- function(U, w, G, noise, lambdamodel, ctx, init, niter) {
  par <- init
  par$D <- .distmatrix(U, par$centers, ctx$dist)
  Gt <- G + noise
  best <- NULL; bestll <- -Inf
  for (it in seq_len(niter)) {
    e <- .estep(U, w, par, ctx)
    if (e$loglik > bestll) { bestll <- e$loglik; best <- par }
    ## passo S: divisione multinomiale delle w_u copie di ogni ranking
    ## distinto, come binomiali condizionate in sequenza
    cnt <- matrix(0, nrow(U), Gt)
    rem <- w; prem <- rep(1, nrow(U))
    for (g in seq_len(Gt - 1)) {
      pr <- pmin(1, pmax(0, e$z[, g] / pmax(prem, 1e-300)))
      cnt[, g] <- stats::rbinom(nrow(U), rem, pr)
      rem <- rem - cnt[, g]; prem <- prem - e$z[, g]
    }
    cnt[, Gt] <- rem
    ## componenti vuote: si tengono in vita con una piccola pseudo-osservazione
    for (g in which(colSums(cnt) == 0)) {
      u <- sample.int(nrow(U), 1, prob = w)
      cnt[u, g] <- cnt[u, g] + 1e-3
    }
    par <- .mstep(U, cnt, par, lambdamodel, ctx)
  }
  best
}

## -----------------------------------------------------------------------------
## [IT] .randominit: un punto di partenza casuale.
##   "partition": partizione casuale "soft" dei dati (pesi esponenziali
##                normalizzati) seguita da un passo M (sui dati APA arriva
##                all'ottimo molto piu' spesso);
##   "centers"  : centri estratti tra i ranking osservati (prob. proporzionale
##                alla frequenza, distinti se possibile), lambda casuale
##                comune, proporzioni uguali.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.randominit <- function(U, w, G, noise, lambdamodel, ctx,
                        init = c("partition", "centers")) {
  init <- match.arg(init)
  M <- ncol(U); Gt <- G + noise
  if (init == "partition") {
    z <- matrix(stats::rexp(nrow(U) * Gt), nrow(U))
    wz <- w * z / rowSums(z)
    centers <- matrix(0L, G, M); lambda <- num <- den <- numeric(G)
    for (g in seq_len(G)) {
      centers[g, ] <- .centerupdate(U, wz[, g], ctx$dist, sample.int(M), ctx)
      num[g] <- sum(wz[, g] * .dist2center(U, centers[g, ], ctx$dist))
      den[g] <- sum(wz[, g])
      lambda[g] <- .solvelambda(num[g] / den[g], M, ctx$dist, ctx$lambdamax)
    }
    if (G > 0 && lambdamodel == "equal")
      lambda[] <- .solvelambda(sum(num) / sum(den), M, ctx$dist, ctx$lambdamax)
    return(list(centers = centers, lambda = lambda, p = colSums(wz) / sum(w)))
  }
  centers <- matrix(0L, G, M)
  if (G > 0) {
    nd <- min(G, nrow(U))
    id <- sample.int(nrow(U), nd, prob = w)
    centers[seq_len(nd), ] <- U[id, ]
    if (G > nd) for (g in (nd + 1):G) centers[g, ] <- sample.int(M)
  }
  list(centers = centers, lambda = rep(stats::runif(1, 0.3, 1.5), G),
       p = rep(1 / Gt, Gt))
}

## -----------------------------------------------------------------------------
## [IT] .nestedinits: partenze per G componenti costruite da una soluzione con
## G - 1 componenti: si tengono i suoi parametri e si aggiunge una componente
## centrata su uno dei ranking osservati piu' "sottostimati" dal modello
## (frequenza osservata - probabilita' stimata massima). Usata da
## mixdbmselect(nested = TRUE) e dall'argomento startfrom di mixdbm().
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.nestedinits <- function(fit, U, w, G, noise, lambdamodel, nnew = 3) {
  fhat <- exp(.logsumexprows(.mixlogdens(fit, U)))
  excess <- w / sum(w) - fhat
  cand <- utils::head(order(-excess), nnew)
  lapply(cand, function(u) {
    pg <- fit$p[seq_len(G - 1)]
    lamnew <- if (lambdamodel == "equal" && G > 1) fit$lambda[1] else 1
    p <- c(0.9 * pg, 0.1, if (noise) 0.9 * fit$p[length(fit$p)])
    list(centers = rbind(unname(fit$centers), U[u, ]),
         lambda = c(unname(fit$lambda), lamnew), p = unname(p / sum(p)))
  })
}

## -----------------------------------------------------------------------------
## [IT] .mixlogdens: matrice n x (G + noise) di log p_g + log f(r | R_g,
## lambda_g) per un modello stimato (serve a predict e a .nestedinits).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.mixlogdens <- function(object, X) {
  Gt <- object$G + object$noise
  lp <- matrix(-lgamma(object$M + 1), nrow(X), Gt)
  for (g in seq_len(object$G))
    lp[, g] <- -object$lambda[g] * .dist2center(X, object$centers[g, ], object$dist) -
      .logpsi(object$lambda[g], object$M, object$dist)
  lp + matrix(log(object$p), nrow(X), Gt, byrow = TRUE)
}

## -----------------------------------------------------------------------------
## [IT] .modelname: etichetta del modello come nel paper, ad es.
## "Cayley-free, 5" o "Kendall-equal, 3+N".
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.modelname <- function(dist, lambdamodel, G, noise) {
  sprintf("%s-%s, %d%s", .distlabel(dist), lambdamodel, G, if (noise) "+N" else "")
}

## -----------------------------------------------------------------------------
## [IT] .degenerate: indici delle componenti degeneri di una soluzione:
## precisione al limite lambdamax e numerosita' attesa n p_g < minsize. Sono
## componenti concentrate su poche copie di un solo ranking, con massimo sul
## bordo (lambda = Inf) dello spazio dei parametri; con molte copie dello
## stesso ranking (n p_g >= minsize) la componente descrive invece un gruppo
## vero e non e' considerata degenere.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.degenerate <- function(fit, n, ctx, control) {
  G <- length(fit$lambda)
  minsize <- if (is.null(control$minsize)) 5 else control$minsize
  which(fit$lambda >= (1 - 1e-6) * ctx$lambdamax & n * fit$p[seq_len(G)] < minsize)
}

## -----------------------------------------------------------------------------
## [IT] .mixdbmfit: il lavoro vero di mixdbm() su dati gia' controllati e
## aggregati (cosi' mixdbmselect() puo' riusare aggregazione e contesto):
## genera le partenze, esegue l'EM da ciascuna, tiene la migliore non
## degenere (o la migliore in assoluto se sono tutte degeneri), calcola
## gradi di liberta', BIC e ICL come nel paper e riordina le componenti per
## proporzione decrescente (la componente uniforme resta ultima).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.mixdbmfit <- function(ag, n, labs, G, lambdamodel, noise, ctx, control,
                       startfrom = NULL) {
  U <- ag$U; w <- ag$w; M <- ctx$M
  nstart <- if (G == 0) 1 else control$nstart
  inits <- lapply(seq_len(nstart), function(s) {
    ini <- switch(control$init, mixed = c("partition", "centers")[2 - s %% 2],
                  sem = "partition", control$init)
    i0 <- .randominit(U, w, G, noise, lambdamodel, ctx, ini)
    if (control$init == "sem" && G > 0)
      i0 <- .semone(U, w, G, noise, lambdamodel, ctx, i0,
                    control$semiter)[c("centers", "lambda", "p")]
    i0
  })
  if (!is.null(startfrom) && G >= 1) {
    if (!inherits(startfrom, "mixdbm") || startfrom$G != G - 1 ||
        startfrom$dist != ctx$dist || startfrom$noise != noise)
      stop("'startfrom' must be a mixdbm fit with G - 1 components, ",
           "the same distance and the same noise option")
    inits <- c(inits, .nestedinits(startfrom, U, w, G, noise, lambdamodel))
  }
  fits <- vector("list", length(inits))
  for (s in seq_along(inits)) {
    fits[[s]] <- .emone(U, w, lambdamodel, ctx, inits[[s]], control$maxiter, control$tol)
    if (control$verbose) cat(sprintf("start %2d: loglik = %.4f (%d it.)\n",
                                     s, fits[[s]]$loglik, fits[[s]]$iter))
  }
  lls <- vapply(fits, `[[`, 0, "loglik")
  ## si tiene il miglior avvio non degenere; se lo sono tutti, il migliore
  deg <- vapply(fits, function(f) length(.degenerate(f, n, ctx, control)) > 0, NA)
  best <- fits[[if (all(deg)) which.max(lls) else which(!deg)[which.max(lls[!deg])]]]

  ## gradi di liberta': proporzioni + precisioni (+ centerdf per centro)
  nlambda <- if (G == 0) 0 else if (lambdamodel == "free") G else 1
  df <- (G + noise - 1) + nlambda + control$centerdf * G
  bic <- 2 * best$loglik - df * log(n)
  zlogz <- best$z * log(best$z); zlogz[!is.finite(zlogz)] <- 0
  icl <- bic + 2 * sum(w * zlogz)

  z <- best$z[ag$idx, , drop = FALSE]
  cn <- c(if (G > 0) paste0("C", seq_len(G)), if (noise) "Noise")
  centers <- best$centers
  if (G > 1) {                                     # ordina per proporzione
    o <- order(best$p[seq_len(G)], decreasing = TRUE)
    centers <- centers[o, , drop = FALSE]
    best$lambda <- best$lambda[o]
    best$p <- c(best$p[o], if (noise) best$p[G + 1])
    z <- z[, c(o, if (noise) G + 1), drop = FALSE]
  }
  colnames(z) <- cn
  dimnames(centers) <- list(cn[seq_len(G)], labs)
  best$p <- stats::setNames(best$p, cn)
  best$lambda <- stats::setNames(best$lambda, cn[seq_len(G)])
  cluster <- max.col(z, ties.method = "first")
  degenerate <- names(best$lambda)[.degenerate(best, n, ctx, control)]

  structure(list(
    centers = centers, ocenters = if (G > 0) .toorderings(centers, labs),
    lambda = best$lambda, p = best$p, z = z, cluster = cluster,
    loglik = best$loglik, df = df, BIC = bic, ICL = icl, n = n, M = M,
    G = G, noise = noise, dist = ctx$dist, lambdamodel = lambdamodel,
    model = .modelname(ctx$dist, lambdamodel, G, noise),
    iter = best$iter, converged = best$converged, lltrace = best$lltrace,
    loglikstarts = lls, degstarts = deg, degenerate = degenerate,
    control = control),
    class = "mixdbm")
}

## [IT] .degeneratemsg: testo dell'avviso sulle componenti degeneri.
#' @keywords internal
#' @noRd
.degeneratemsg <- function(fit) {
  g <- fit$degenerate
  sprintf(paste0("%s: all %d starts gave a degenerate solution. Component(s) %s ",
                 "(expected size %s) have precision at lambdamax = %g and are ",
                 "concentrated on a few identical rankings; the log-likelihood, BIC ",
                 "and ICL are not reliable. Consider fewer components."),
          fit$model, length(fit$degstarts), paste(g, collapse = ", "),
          paste(format(fit$n * fit$p[g], digits = 2), collapse = ", "),
          fit$control$lambdamax)
}
