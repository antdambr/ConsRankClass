## -----------------------------------------------------------------------------
## [IT] .ncycles: numero di cicli di ogni riga di Y, vista come permutazione
## j -> Y[, j]. Vettorizzata su tutte le righe: per ogni elemento si segue la
## sua orbita per M - 1 passi tenendo il minimo incontrato; ogni ciclo ha un
## solo "capo" (il suo elemento minimo). Con labels = TRUE restituisce, per
## ogni elemento, l'etichetta del ciclo a cui appartiene (il suo minimo):
## serve alla ricerca locale veloce del centro di Cayley.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.ncycles <- function(Y, labels = FALSE) {
  n <- nrow(Y); M <- ncol(Y)
  y <- as.integer(Y)                       # per colonne: Y[i, j] = y[i + (j-1) n]
  base <- rep.int(seq_len(n), M)
  start <- rep(seq_len(M), each = n)
  cur <- mn <- start
  for (s in seq_len(M - 1)) {
    cur <- y[base + (cur - 1L) * n]
    mn <- pmin.int(mn, cur)
  }
  if (labels) return(matrix(mn, n, M))
  .rowSums(mn == start, n, M)
}

## -----------------------------------------------------------------------------
## [IT] .lis: lunghezza della piu' lunga sottosequenza crescente di ogni riga
## di Y (algoritmo "patience sorting", vettorizzato sulle righe: Tp[i, ]
## contiene le cime delle pile, crescenti; ogni nuovo valore va sulla prima
## pila con cima maggiore). Serve per la distanza di Ulam.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.lis <- function(Y) {
  n <- nrow(Y); M <- ncol(Y)
  Tp <- matrix(Inf, n, M)
  for (j in seq_len(M)) {
    x <- Y[, j]
    pos <- .rowSums(Tp < x, n, M) + 1L
    Tp[cbind(seq_len(n), pos)] <- x
  }
  .rowSums(is.finite(Tp), n, M)
}

## -----------------------------------------------------------------------------
## [IT] .distid: distanza di ogni riga di Y (permutazione relativa) dalla
## permutazione identita'. E' il "motore" di tutte le distanze:
##   kendall    = numero di inversioni (coppie k < l con Y[k] > Y[l]);
##   cayley     = M - numero di cicli (minimo numero di scambi);
##   hamming    = numero di posizioni con Y[k] != k;
##   ulam       = M - piu' lunga sottosequenza crescente (minimo numero di
##                spostamenti "togli e reinserisci" di un oggetto);
##   footrule   = sum |Y[k] - k|;
##   spearman   = sqrt(sum (Y[k] - k)^2);
##   spearmansq = sum (Y[k] - k)^2.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.distid <- function(Y, dist) {
  n <- nrow(Y); M <- ncol(Y)
  K <- matrix(seq_len(M), n, M, byrow = TRUE)
  switch(dist,
    kendall = {
      d <- numeric(n)
      if (M > 1) for (k in 1:(M - 1)) for (l in (k + 1):M)
        d <- d + (Y[, k] > Y[, l])
      d
    },
    cayley     = M - .ncycles(Y),
    hamming    = .rowSums(Y != K, n, M),
    ulam       = M - .lis(Y),
    footrule   = .rowSums(abs(Y - K), n, M),
    spearman   = sqrt(.rowSums((Y - K)^2, n, M)),
    spearmansq = .rowSums((Y - K)^2, n, M))
}

## -----------------------------------------------------------------------------
## [IT] .dist2center: distanza di ogni riga di X dal singolo ranking R.
## Si rietichettano gli oggetti nell'ordine di R (Y = X[, order(R)]), cosi'
## che R diventi l'identita', e si usa .distid (invarianza a destra).
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.dist2center <- function(X, R, dist) {
  .distid(X[, order(R), drop = FALSE], dist)
}
