## -----------------------------------------------------------------------------
## [IT] .kendallW: matrice delle preferenze pesate W[a, b] = somma dei pesi
## dei ranking che mettono l'oggetto a prima di b. L'obiettivo di Kendall
## per un candidato R e' la somma di W[a, b] sulle coppie che R ordina al
## contrario (.kendallobj).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.kendallW <- function(U, wz) {
  ## colonna b: W[a, b] = sum(wz * (U[, a] < U[, b])) per tutti gli a
  vapply(seq_len(ncol(U)), function(b) colSums(wz * (U < U[, b])), numeric(ncol(U)))
}

#' @keywords internal
#' @noRd
.kendallobj <- function(R, W) sum(W[outer(R, R, ">")])

## -----------------------------------------------------------------------------
## [IT] .neighbours: intorno di un ranking per la ricerca locale: tutti gli
## scambi di due oggetti e tutti gli spostamenti di un oggetto in un'altra
## posizione dell'ordinamento (le mosse naturali di Cayley e di Ulam).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.neighbours <- function(R) {
  M <- length(R)
  o <- order(R)
  out <- list()
  for (i in 1:(M - 1)) for (j in (i + 1):M) {
    s <- R; s[c(i, j)] <- R[c(j, i)]; out[[length(out) + 1]] <- s
  }
  for (a in 1:M) for (b in 1:M) if (abs(a - b) > 1) {
    oo <- append(o[-a], o[a], after = b - 1)
    out[[length(out) + 1]] <- order(oo)
  }
  do.call(rbind, out)
}

## -----------------------------------------------------------------------------
## [IT] .hungarian: problema di assegnamento (algoritmo ungherese,
## Kuhn-Munkres, versione O(M^3) con potenziali). Data una matrice quadrata
## C[j, r] = costo di dare la posizione r all'oggetto j, restituisce il
## ranking R che minimizza sum_j C[j, R[j]]. E' il centro ESATTO per footrule
## e Hamming, la cui somma pesata delle distanze si scompone oggetto per
## oggetto.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.hungarian <- function(C) {
  n <- nrow(C)
  u <- numeric(n + 1); v <- numeric(n + 1)
  p <- integer(n + 1); way <- integer(n + 1)       # indici 0..n -> 1..n+1
  for (i in 1:n) {
    p[1] <- i; j0 <- 0L
    minv <- rep(Inf, n + 1); used <- rep(FALSE, n + 1)
    repeat {
      used[j0 + 1] <- TRUE
      i0 <- p[j0 + 1]
      free <- which(!used[-1])                     # colonne j = 1..n libere
      cur <- C[i0, free] - u[i0 + 1] - v[free + 1]
      upd <- cur < minv[free + 1]
      minv[free[upd] + 1] <- cur[upd]
      way[free[upd] + 1] <- j0
      k <- which.min(minv[free + 1])
      delta <- minv[free[k] + 1]; j1 <- free[k]
      uj <- which(used) - 1                        # colonne usate (0..n)
      u[p[uj + 1] + 1] <- u[p[uj + 1] + 1] + delta
      v[uj + 1] <- v[uj + 1] - delta
      minv[!used] <- minv[!used] - delta
      j0 <- j1
      if (p[j0 + 1] == 0) break
    }
    repeat {
      j1 <- way[j0 + 1]; p[j0 + 1] <- p[j1 + 1]; j0 <- j1
      if (j0 == 0) break
    }
  }
  R <- integer(n)
  for (j in 1:n) R[p[j + 1]] <- j
  R
}

## -----------------------------------------------------------------------------
## [IT] .centerupdate: passo M per il centro di una componente. Dati i
## ranking distinti U, i pesi wz e il centro corrente Rold, restituisce un
## centro con somma pesata delle distanze non superiore (l'EM resta monotono;
## a parita' di obiettivo si tiene Rold, per stabilita'). Metodi:
##   spearmansq      : Borda pesato (ordinare le posizioni medie), esatto;
##   footrule/hamming: problema di assegnamento (.hungarian), esatto;
##   kendall, M > 8  : branch-and-bound di Emond & Mason con
##                     ConsRank::consrank(algorithm = "BB", full = TRUE), esatto;
##   M piccolo       : ricerca esaustiva su tutte le M! permutazioni (Kendall
##                     M <= 8, altrimenti M <= 7), esatta;
##   altrimenti      : ricerca locale (scambi e spostamenti) da piu' partenze;
##                     per Cayley gli scambi si valutano tutti in O(n M^2).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.centerupdate <- function(U, wz, dist, Rold, ctx) {
  M <- ncol(U)
  obj <- function(R) sum(wz * .dist2center(U, R, dist))
  better <- function(R) if (obj(R) < obj(Rold) - 1e-9 * (1 + abs(obj(Rold)))) R else Rold
  if (dist == "spearmansq")                        # Borda pesato
    return(better(rank(colSums(wz * U), ties.method = "first")))
  if (dist %in% c("footrule", "hamming")) {        # assegnamento
    C <- outer(seq_len(M), seq_len(M), Vectorize(function(j, r)
      if (dist == "footrule") sum(wz * abs(U[, j] - r)) else sum(wz * (U[, j] != r))))
    return(better(.hungarian(C)))
  }
  if (dist == "kendall" && isTRUE(ctx$bb)) {
    ## Pesi riscalati a interi: con pesi non interi il branch-and-bound di
    ## ConsRank 2.1.4 puo' fallire (per arrotondamenti tutti i rami di un giro
    ## vengono potati e BBconsensus2 si ferma con "argomento di lunghezza 0").
    wi <- round(wz / max(wz) * 1e8)
    keep <- wi > 0
    S <- tryCatch(
      ConsRank::consrank(U[keep, , drop = FALSE], wk = wi[keep],
                         algorithm = "BB", full = TRUE, ps = FALSE)$Consensus,
      error = function(e) NULL)
    if (is.null(S)) {
      warning("ConsRank branch-and-bound failed; local search used instead")
      return(.centerupdate(U, wz, dist, Rold, utils::modifyList(ctx, list(bb = FALSE))))
    }
    S <- matrix(as.integer(S), ncol = M)
    if (any(apply(S, 1, function(s) all(s == Rold)))) return(Rold)
    W <- .kendallW(U, wz)
    if (.kendallobj(S[1, ], W) < .kendallobj(Rold, W) - 1e-9) return(S[1, ])
    return(Rold)
  }
  ## obiettivo per molti candidati insieme (righe di P)
  if (dist == "kendall") {
    W <- .kendallW(U, wz)
    objmany <- function(P) {
      v <- numeric(nrow(P))
      for (a in 1:M) for (b in 1:M) if (a != b) v <- v + W[a, b] * (P[, a] > P[, b])
      v
    }
  } else if (dist == "cayley") {
    objmany <- function(P) {
      Y <- do.call(rbind, lapply(seq_len(nrow(P)), function(k)
        U[, order(P[k, ]), drop = FALSE]))
      colSums(matrix(wz * (M - .ncycles(Y)), nrow(U)))
    }
  } else {
    objmany <- function(P) apply(P, 1, obj)
  }
  if (ctx$exhaustive) {
    vals <- if (!is.null(ctx$Dall)) {
      drop(crossprod(ctx$Dall, wz))
    } else if (!is.null(ctx$Kall)) {
      ## coppia a < b: costa W[b, a] se R mette a prima di b, W[a, b] altrimenti
      wab <- W[ctx$pairs]; wba <- W[ctx$pairs[, 2:1, drop = FALSE]]
      sum(wba) + drop(ctx$Kall %*% (wab - wba))
    } else objmany(ctx$perms)
    best <- which.min(vals)
    old <- vals[match(paste(Rold, collapse = ","), ctx$permkey)]
    return(if (vals[best] < old - 1e-12) ctx$perms[best, ] else Rold)
  }
  ## ricerca locale (EM generalizzato: mai peggio di Rold). Partenze: centro
  ## corrente, Borda pesato, i 5 ranking osservati con peso maggiore.
  starts <- rbind(Rold, rank(colSums(wz * U), ties.method = "first"),
                  U[utils::head(order(-wz), 5), , drop = FALSE])
  sv <- objmany(starts)
  R <- starts[which.min(sv), ]; cur <- min(sv)
  if (dist == "cayley") {
    ## Scambiare le posizioni di due oggetti a, b compone Y con la
    ## trasposizione (a b): d scende di 1 se a e b stanno nello stesso ciclo
    ## di Y, sale di 1 altrimenti. Tutti gli scambi valutati in O(n M^2).
    sw <- sum(wz)
    repeat {
      L <- .ncycles(U[, order(R), drop = FALSE], labels = TRUE)
      best <- 0; ab <- NULL
      for (a in 1:(M - 1)) for (b in (a + 1):M) {
        delta <- sw - 2 * sum(wz[L[, a] == L[, b]])
        if (delta < best - 1e-12) { best <- delta; ab <- c(a, b) }
      }
      if (is.null(ab)) break
      o <- order(R)[ab]
      R[o] <- R[rev(o)]
    }
    return(better(unname(R)))
  }
  repeat {
    Nb <- .neighbours(R)
    v <- objmany(Nb)
    if (min(v) < cur - 1e-12) { R <- Nb[which.min(v), ]; cur <- min(v) } else break
  }
  better(unname(R))
}

## -----------------------------------------------------------------------------
## [IT] .makectx: prepara il "contesto" della stima del centro, cioe' le scelte
## che non cambiano durante l'EM: metodo (esaustivo / branch-and-bound /
## locale), elenco delle permutazioni e, se la memoria lo consente, la matrice
## delle distanze di tutti i ranking distinti da tutte le permutazioni (cosi'
## il passo M esaustivo e' un solo prodotto matrice-vettore).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.makectx <- function(U, M, dist, centermethod = "auto", lambdamax = 50) {
  centermethod <- match.arg(centermethod, c("auto", "exhaustive", "local", "bb"))
  if (centermethod == "bb" && dist != "kendall")
    stop("centermethod = \"bb\" is available for the Kendall distance only")
  if (centermethod == "auto" && dist == "kendall" && M > 8) centermethod <- "bb"
  exhaustive <- switch(centermethod, exhaustive = TRUE, local = FALSE, bb = FALSE,
                       auto = (dist == "kendall" && M <= 8) || M <= 7)
  if (dist %in% c("spearmansq", "footrule", "hamming")) exhaustive <- FALSE  # esatti
  if (exhaustive && M > 9) stop("exhaustive search of the center requires M <= 9")
  perms <- Dall <- Kall <- pairs <- NULL
  if (exhaustive) {
    perms <- .allperms(M)
    if (dist == "kendall") {
      ## Kall[k, q] = 1 se la permutazione k mette l'oggetto pairs[q, 1] dopo
      ## pairs[q, 2] (coppie a < b): l'obiettivo di Kendall di tutte le M!
      ## permutazioni diventa un prodotto matrice-vettore (.centerupdate).
      pairs <- which(upper.tri(diag(M)), arr.ind = TRUE)
      Kall <- 1 * (perms[, pairs[, 1], drop = FALSE] > perms[, pairs[, 2], drop = FALSE])
    }
    if (dist != "kendall" && as.numeric(nrow(U)) * nrow(perms) <= 2e7)
      Dall <- vapply(seq_len(nrow(perms)), function(k)
        .dist2center(U, perms[k, ], dist), numeric(nrow(U)))
  }
  list(M = M, dist = dist, exhaustive = exhaustive, bb = centermethod == "bb",
       perms = perms, Kall = Kall, pairs = pairs,
       permkey = if (exhaustive) apply(perms, 1, paste, collapse = ","),
       Dall = if (!is.null(Dall)) matrix(Dall, nrow(U)) else NULL,
       lambdamax = lambdamax)
}

## -----------------------------------------------------------------------------
## [IT] .centerestimate: centro di un insieme pesato di ranking (non dentro
## l'EM). Si parte dal Borda pesato e si applica .centerupdate finche' il
## centro non cambia piu' (con i metodi esatti basta un passo).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.centerestimate <- function(U, w, dist, ctx) {
  R <- rank(colSums(w * U), ties.method = "first")
  repeat {
    Rn <- .centerupdate(U, w, dist, R, ctx)
    if (all(Rn == R)) break
    R <- Rn
  }
  R
}
