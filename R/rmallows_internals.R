## -----------------------------------------------------------------------------
## [IT] .randsyt: tableau di Young standard UNIFORME di forma lam, con la
## "passeggiata sugli uncini" di Greene, Nijenhuis & Wilf (1979): si sceglie
## una cella a caso e ci si sposta a caso nel suo uncino finche' non si
## arriva a un angolo, dove si mette il numero piu' grande rimasto.
## Restituisce il tableau come lista di righe.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.randsyt <- function(lam) {
  N <- sum(lam); sh <- lam
  Tb <- lapply(lam, function(l) integer(l))
  for (t in N:1) {
    rows <- rep(seq_along(sh), sh); cols <- sequence(sh)
    u <- sample.int(length(rows), 1)
    i <- rows[u]; j <- cols[u]
    repeat {
      arm <- sh[i] - j; leg <- sum(sh >= j) - i
      if (arm + leg == 0) break
      v <- sample.int(arm + leg, 1)
      if (v <= arm) j <- j + v else i <- i + (v - arm)
    }
    Tb[[i]][j] <- t
    sh[i] <- sh[i] - 1L
  }
  Tb
}

## -----------------------------------------------------------------------------
## [IT] .invrsk: corrispondenza di Robinson-Schensted inversa: da una coppia
## (P, Q) di tableau standard della stessa forma ricostruisce la permutazione
## w (come sequenza). La lunghezza della sottosequenza crescente massima di w
## e' la lunghezza della prima riga della forma.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.invrsk <- function(P, Q) {
  N <- sum(lengths(P)); w <- integer(N)
  for (t in N:1) {
    r <- which(vapply(Q, function(q) length(q) > 0 && q[length(q)] == t, TRUE))
    L <- length(P[[r]]); x <- P[[r]][L]
    P[[r]] <- P[[r]][-L]; Q[[r]] <- Q[[r]][-L]
    if (r > 1) for (rr in (r - 1):1) {
      row <- P[[rr]]
      pos <- max(which(row < x))
      y <- row[pos]; row[pos] <- x; P[[rr]] <- row; x <- y
    }
    w[t] <- x
  }
  w
}

## -----------------------------------------------------------------------------
## [IT] .rdistance: estrae n valori della distanza dalla sua distribuzione
## sotto il modello di Mallows, P(d) proporzionale a N(d) exp(-lambda d).
## Restituisce gli INDICI nella tabella .distfreq(M, dist).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.rdistance <- function(n, M, lambda, dist) {
  dd <- .distfreq(M, dist)
  a <- dd$logN - lambda * dd$d
  sample.int(length(a), n, replace = TRUE, prob = exp(a - max(a)))
}

## -----------------------------------------------------------------------------
## [IT] .rdbmidentity: n permutazioni (in notazione ranking) con probabilita'
## proporzionale a exp(-lambda d(sigma, identita')). Metodi ESATTI:
##   kendall : vettore delle inversioni, componenti indipendenti con legge
##             geometrica troncata (Fligner & Verducci, 1986);
##   cayley  : il modello e' la legge di Ewens con theta = e^lambda, simulata
##             con il "ristorante cinese";
##   hamming : si estrae d, poi d posizioni a caso e una dismutazione uniforme
##             di esse (tutte le permutazioni con lo stesso d sono equiprobabili);
##   ulam    : si estrae d, poi la forma lambda (prob. prop. a f_lambda^2, prima
##             riga = M - d), due tableau uniformi e la RS inversa;
##   footrule, spearman(sq): si estrae il costo s, poi una permutazione
##             uniforme di costo s, posizione per posizione, con i conteggi
##             della programmazione dinamica (M <= 14).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.rdbmidentity <- function(n, M, lambda, dist) {
  O <- matrix(0L, n, M)
  if (n == 0) return(O)
  if (dist == "kendall") {
    for (i in seq_len(n)) {
      rem <- seq_len(M); o <- integer(M)
      for (j in seq_len(M)) {
        k <- M - j + 1
        v <- sample.int(k, 1, prob = exp(-lambda * (0:(k - 1)))) - 1
        o[j] <- rem[v + 1]; rem <- rem[-(v + 1)]
      }
      O[i, ] <- order(o)
    }
    return(O)
  }
  if (dist == "cayley") {
    th <- exp(lambda)
    for (i in seq_len(n)) {
      s <- 1L
      if (M > 1) for (k in 2:M) {
        s[k] <- k
        if (stats::runif(1) > th / (th + k - 1)) {
          j <- sample.int(k - 1, 1)
          s[k] <- s[j]; s[j] <- k
        }
      }
      O[i, ] <- s
    }
    return(O)
  }
  dd <- .distfreq(M, dist)
  id <- .rdistance(n, M, lambda, dist)
  if (dist == "hamming") {
    for (i in seq_len(n)) {
      d <- dd$d[id[i]]; y <- seq_len(M)
      if (d > 0) {
        S <- if (d == M) seq_len(M) else sort(sample.int(M, d))
        repeat { p <- sample.int(d); if (all(p != seq_len(d))) break }
        y[S] <- S[p]
      }
      O[i, ] <- y
    }
    return(O)
  }
  if (dist == "ulam") {
    tab <- .ulamtab(M)
    for (i in seq_len(n)) {
      k <- M - dd$d[id[i]]
      cand <- which(tab$l1 == k)
      lw <- 2 * tab$logf[cand]
      lam <- tab$parts[[cand[sample.int(length(cand), 1, prob = exp(lw - max(lw)))]]]
      O[i, ] <- .invrsk(.randsyt(lam), .randsyt(lam))
    }
    return(O)
  }
  ## footrule, spearman, spearmansq: conteggi completi della DP
  cost <- .costmatrix(M, if (dist == "footrule") "footrule" else "spearman")
  dp <- .permcostdp(M, cost, keep = TRUE)
  full <- 2^M - 1
  for (i in seq_len(n)) {
    s <- dd$s[id[i]]; mask <- full; y <- integer(M)
    for (k in M:1) {
      vv <- which(bitwAnd(mask, 2^(seq_len(M) - 1)) > 0)
      cnt <- vapply(vv, function(v) {
        ss <- s - cost[k, v]
        if (ss < 0) return(0)
        fm <- dp$f[[mask - 2^(v - 1) + 1]]
        if (is.null(fm)) 0 else fm[ss + 1]
      }, 0)
      v <- vv[sample.int(length(vv), 1, prob = cnt)]
      y[k] <- v; s <- s - cost[k, v]; mask <- mask - 2^(v - 1)
    }
    O[i, ] <- y
  }
  O
}
