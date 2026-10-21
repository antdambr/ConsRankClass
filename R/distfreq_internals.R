## [IT] Ambiente di cache: le distribuzioni N(d) (e le partizioni per Ulam) si
## calcolano una sola volta per (distanza, M) in ogni sessione di R.
.mixdbmcache <- new.env(parent = emptyenv())

## -----------------------------------------------------------------------------
## [IT] .permcostdp: programmazione dinamica su sottoinsiemi (bitmask) per le
## distanze ADDITIVE sulle posizioni, d(Y, e) = sum_k cost[k, Y[k]]
## (footrule: |k - v|; Spearman al quadrato: (k - v)^2). Si riempiono le
## posizioni 1, 2, ... in ordine; lo stato e' l'insieme dei valori gia' usati
## (mask) e f[[mask]] e' il vettore dei conteggi per ogni costo parziale.
## Costo O(2^M M S), con S il costo massimo: per questo M <= 14.
## Con keep = TRUE conserva tutta la tabella (serve per simulare, sezione 4).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.permcostdp <- function(M, cost, keep = FALSE) {
  if (M > 14) stop("exact computation for this distance requires M <= 14 ",
                   "(2^M states); use dist = \"kendall\" or \"cayley\"")
  maxS <- sum(apply(cost, 1, max))
  f <- vector("list", 2^M)
  f[[1]] <- c(1, numeric(maxS))
  pc <- integer(2^M)                             # numero di bit di ogni mask
  if (M > 0) for (mask in 1:(2^M - 1)) pc[mask + 1] <- pc[(mask %/% 2) + 1] + mask %% 2
  for (mask in 0:(2^M - 2)) {
    cur <- f[[mask + 1]]
    if (is.null(cur)) next
    k <- pc[mask + 1] + 1                         # posizione da riempire
    for (v in seq_len(M)) {
      b <- 2^(v - 1)
      if ((mask %/% b) %% 2 == 1) next
      a <- cost[k, v]
      sh <- if (a == 0) cur else c(numeric(a), cur[seq_len(maxS + 1 - a)])
      nm <- mask + b + 1
      f[[nm]] <- if (is.null(f[[nm]])) sh else f[[nm]] + sh
    }
    if (!keep) f[mask + 1] <- list(NULL)
  }
  list(N = f[[2^M]], f = if (keep) f, maxS = maxS)
}

## -----------------------------------------------------------------------------
## [IT] .costmatrix: matrice cost[k, v] delle distanze additive (footrule e
## Spearman, quest'ultima prima della radice).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.costmatrix <- function(M, dist) {
  D <- outer(seq_len(M), seq_len(M), "-")
  if (dist == "footrule") abs(D) else D^2
}

## -----------------------------------------------------------------------------
## [IT] .npartitions e .partitions: numero e elenco delle partizioni intere di
## M (per Ulam). L'elenco e' generato con l'algoritmo ZS1 di Zoghbi &
## Stojmenovic (1998), in ordine anti-lessicografico; ogni partizione e' un
## vettore decrescente di interi positivi.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.npartitions <- function(M) {
  p <- c(1, numeric(M))
  for (k in seq_len(M)) for (s in k:M) p[s + 1] <- p[s + 1] + p[s - k + 1]
  p[M + 1]
}

#' @keywords internal
#' @noRd
.partitions <- function(M) {
  M <- as.integer(M)
  out <- vector("list", .npartitions(M))
  x <- rep(1L, M); x[1] <- M; m <- 1L; h <- 1L
  cnt <- 1L; out[[1]] <- x[1]
  while (x[1] != 1L) {
    if (x[h] == 2L) {
      m <- m + 1L; x[h] <- 1L; h <- h - 1L
    } else {
      r <- x[h] - 1L; t <- m - h + 1L; x[h] <- r
      while (t >= r) { h <- h + 1L; x[h] <- r; t <- t - r }
      if (t == 0L) {
        m <- h
      } else {
        m <- h + 1L
        if (t > 1L) { h <- h + 1L; x[h] <- t }
      }
    }
    cnt <- cnt + 1L
    out[[cnt]] <- x[seq_len(m)]
  }
  out
}

## -----------------------------------------------------------------------------
## [IT] .ulamtab: per ogni partizione lambda di M calcola log f_lambda, dove
## f_lambda = numero di tableau di Young standard di forma lambda (formula
## degli uncini: f = M! / prod(uncini)). Per la corrispondenza di
## Robinson-Schensted, il numero di permutazioni la cui sottosequenza
## crescente massima e' lunga k e' sum_{lambda: lambda_1 = k} f_lambda^2.
## Risultato in cache (serve anche per simulare, sezione 4).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.ulamtab <- function(M) {
  if (M > 50) stop("Ulam distance distribution available for M <= 50")
  key <- paste("ulamtab", M)
  if (!is.null(.mixdbmcache[[key]])) return(.mixdbmcache[[key]])
  parts <- .partitions(M)
  logf <- vapply(parts, function(lam) {
    i <- rep(seq_along(lam), lam)
    j <- sequence(lam)
    conj <- vapply(seq_len(lam[1]), function(c) sum(lam >= c), 0)
    hook <- lam[i] - j + conj[j] - i + 1
    lgamma(M + 1) - sum(log(hook))
  }, 0)
  res <- list(parts = parts, l1 = vapply(parts, function(p) as.integer(p[1]), 0L), logf = logf)
  assign(key, res, envir = .mixdbmcache)
  res
}

## -----------------------------------------------------------------------------
## [IT] .logderange: log del numero di dismutazioni (permutazioni senza punti
## fissi) di d elementi: D_d = d! sum_{k=0}^d (-1)^k / k!. Per Hamming
## N(d) = choose(M, d) D_d. Per d = 1 D_d = 0 (log = -Inf).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.logderange <- function(d) {
  vapply(d, function(x) {
    k <- 0:x
    s <- sum((-1)^k / factorial(pmin(k, 170)) * (k <= 170))
    if (s <= 0) -Inf else lgamma(x + 1) + log(s)
  }, 0)
}

## -----------------------------------------------------------------------------
## [IT] .distfreq: distribuzione esatta della distanza dall'identita' sulle M!
## permutazioni: valori d, conteggi N (se rappresentabili) e log N. Metodi:
##   kendall  : numeri di Mahonian (prodotto dei polinomi 1 + q + ... + q^(k-1));
##   cayley   : numeri di Stirling di prima specie senza segno (N(d) =
##              permutazioni con M - d cicli);
##   hamming  : choose(M, d) * dismutazioni(d), qualunque M;
##   ulam     : somme di f_lambda^2 sulle partizioni (M <= 50);
##   footrule, spearman, spearmansq: programmazione dinamica (M <= 14).
## Solo i valori con N > 0 vengono restituiti.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.distfreq <- function(M, dist) {
  key <- paste("freq", dist, M)
  if (!is.null(.mixdbmcache[[key]])) return(.mixdbmcache[[key]])
  res <- switch(dist,
    kendall = {
      if (M > 170) stop("Kendall distribution: M <= 170")
      N <- 1
      if (M > 1) for (k in 2:M) {
        L <- length(N); new <- numeric(L + k - 1)
        for (i in 0:(k - 1)) new[(1:L) + i] <- new[(1:L) + i] + N
        N <- new
      }
      list(d = seq_along(N) - 1, logN = log(N))
    },
    cayley = {
      if (M > 170) stop("Cayley distribution: M <= 170")
      s <- 1
      if (M > 1) for (k in 1:(M - 1)) s <- c(k * s, 0) + c(0, s)
      list(d = 0:(M - 1), logN = log(rev(s)))     # s[c] = perm. con c cicli
    },
    hamming = {
      d <- 0:M
      list(d = d, logN = lchoose(M, d) + .logderange(d))
    },
    ulam = {
      tab <- .ulamtab(M)
      k <- seq_len(M)
      logS <- vapply(k, function(kk) .logsumexp(2 * tab$logf[tab$l1 == kk]), 0)
      list(d = rev(M - k), logN = rev(logS))       # d = M - k
    },
    {                                              # footrule, spearman(sq)
      dp <- .permcostdp(M, .costmatrix(M, if (dist == "footrule") "footrule" else "spearman"))
      s <- seq_along(dp$N) - 1
      list(d = if (dist == "spearman") sqrt(s) else s, logN = log(dp$N), s = s)
    })
  keep <- is.finite(res$logN)
  res <- lapply(res, function(v) v[keep])
  res$N <- exp(res$logN)
  res$M <- M; res$dist <- dist
  assign(key, res, envir = .mixdbmcache)
  res
}
