## =============================================================================
## mixdbm_new.R
##
## Distances between rankings, Mallows (distance-based) models and mixtures of
## distance-based models (Murphy & Martin, 2003) for complete rankings.
## Written to be merged into the package ConsRankClass (which Depends on
## ConsRank): roxygen2 documentation, ConsRankClass naming style (lower case,
## no underscores), S3 methods print/summary/logLik/predict, control object
## as in cca()/ccacontrol().
##
## [IT] Il file contiene tre livelli, ciascuno costruito sul precedente:
##   1. DISTANZE tra ranking completi (Kendall, Cayley, Hamming, Ulam,
##      footrule, Spearman e Spearman al quadrato). Sono usabili da sole e da
##      qualunque metodo del package (anche non parametrico).
##   2. MODELLO DI MALLOWS (modello basato su una distanza) con un solo centro:
##      costante di normalizzazione, valore atteso della distanza, stima di
##      lambda e del centro, densita', simulazione (l'equivalente essenziale di
##      PerMallows, ma con pesi e con tutte le distanze sopra).
##   3. MISTURE di modelli di Mallows (Murphy & Martin, 2003): EM con piu'
##      avvii, vincoli sulle precisioni, componente uniforme, BIC e ICL.
##
## Notazione: un ranking r = (r_1, ..., r_M) contiene in r_j la POSIZIONE
## dell'oggetto j (1 = preferito). I dati sono una matrice n x M, un ranking
## per riga, come in ConsRank/ConsRankClass. Pareggi non ammessi (per ora).
##
## [IT] Idea chiave usata ovunque: tutte le distanze implementate sono
## invarianti a destra, cioe' d(r, R) dipende solo dalla permutazione "relativa"
## Y = r[order(R)] (Y[k] = posizione in r dell'oggetto che R mette al k-esimo
## posto). Quindi d(r, R) = d(Y, identita') e la costante di normalizzazione
## del modello di Mallows non dipende dal centro R: basta conoscere la
## distribuzione di d(sigma, identita') sulle M! permutazioni.
##
## Per aggiungere una distanza bisogna intervenire SOLO nella sezione 2
## (".distid", ".distfreq" e, se esistono formule chiuse, ".logpsi",
## ".expectd", ".vard"), nella sezione 4 (".rdbmidentity", campionamento) e
## nella sezione 5 (".centerupdate", stima del centro). Tutto il resto
## (Mallows, misture, selezione dei modelli) funziona senza modifiche.
##
## Dipendenze: solo R base (stats) e, per il centro di Kendall con M > 8,
## ConsRank::consrank (branch-and-bound di Emond & Mason, 2002).
## =============================================================================



## -----------------------------------------------------------------------------
## [IT] .rankdistances: elenco delle distanze disponibili, nell'ordine in cui
## compaiono negli argomenti "dist". "spearman" e' la radice della somma dei
## quadrati (come nell'eq. (1) di Murphy & Martin, 2003); "spearmansq" e' la
## somma dei quadrati senza radice (la "Spearman's rho distance" di molti
## autori, ad es. Critchlow, Marden, e quella usata dal pacchetto MSmix).
## -----------------------------------------------------------------------------
.rankdistances <- c("kendall", "cayley", "hamming", "ulam",
                    "footrule", "spearman", "spearmansq")

## -----------------------------------------------------------------------------
## [IT] .matchdist: controlla e completa il nome di una distanza (accetta
## abbreviazioni non ambigue, ad es. "ken"). Restituisce il nome completo.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.matchdist <- function(dist) {
  if (length(dist) != 1 || !is.character(dist))
    stop("'dist' must be one of: ", paste(.rankdistances, collapse = ", "))
  match.arg(tolower(dist), .rankdistances)
}

## -----------------------------------------------------------------------------
## [IT] .checkrank: verifica che X sia una matrice (o un vettore, trattato come
## una riga) di ranking COMPLETI, cioe' che ogni riga sia una permutazione di
## 1..M. Se trova pareggi lo dice esplicitamente (i dati con pareggi non sono
## ancora gestiti). Restituisce la matrice con valori interi e con i nomi
## delle colonne (gli oggetti) conservati.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.checkrank <- function(X, what = "X") {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (is.null(dim(X))) X <- matrix(X, 1, dimnames = list(NULL, names(X)))
  if (!is.numeric(X)) stop("'", what, "' must be numeric")
  if (anyNA(X)) stop("'", what, "' contains missing values")
  M <- ncol(X)
  ok <- apply(X, 1, function(r) all(sort(r) == seq_len(M)))
  if (!all(ok)) {
    bad <- which(!ok)
    tied <- any(apply(X[bad, , drop = FALSE], 1, anyDuplicated) > 0)
    stop("rows ", paste(utils::head(bad), collapse = ", "),
         if (length(bad) > 6) ", ..." else "",
         " of '", what, "' are not complete rankings of 1..", M,
         if (tied) " (rankings with ties are not supported yet)" else "")
  }
  storage.mode(X) <- "integer"
  X
}

## -----------------------------------------------------------------------------
## [IT] .aggregaterank: raggruppa i ranking identici. Restituisce U (ranking
## distinti), w (frequenze o somma dei pesi) e idx (per ogni osservazione
## originale, la riga di U corrispondente). Tutti i calcoli (EM compreso) si
## fanno sui ranking distinti: con M piccolo e n grande (ad es. APA: 5738
## osservazioni, 120 ranking distinti) il guadagno di tempo e' enorme.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.aggregaterank <- function(X, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, nrow(X))
  key <- apply(X, 1, paste, collapse = ",")
  first <- !duplicated(key)
  idx <- match(key, key[first])
  w <- as.numeric(tapply(weights, factor(idx, levels = seq_len(sum(first))), sum))
  list(U = X[first, , drop = FALSE], w = w, idx = idx)
}

## -----------------------------------------------------------------------------
## [IT] .allperms: tutte le M! permutazioni di 1..M, una per riga (M <= 10).
## Costruzione per inserimento: le permutazioni di 1..k si ottengono inserendo
## k in ogni posizione delle permutazioni di 1..k-1. Usata per la ricerca
## esaustiva del centro e per le verifiche.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.allperms <- function(M) {
  if (M > 10) stop("too many permutations: M must be <= 10")
  P <- matrix(1L, 1, 1)
  if (M == 1) return(P)
  for (k in 2:M) {
    P <- do.call(rbind, lapply(seq_len(k), function(pos)
      cbind(P[, seq_len(pos - 1), drop = FALSE], k,
            P[, seq_len(k - 1) >= pos, drop = FALSE])))
  }
  unname(P)
}

## -----------------------------------------------------------------------------
## [IT] .logsumexp e .logsumexprows: log(sum(exp(a))) calcolato in modo stabile
## (si sottrae il massimo), per un vettore e per ogni riga di una matrice.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.logsumexp <- function(a) {
  m <- max(a)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(a - m)))
}

#' @keywords internal
#' @noRd
.logsumexprows <- function(A) {
  m <- apply(A, 1, max)
  m[!is.finite(m)] <- 0
  m + log(rowSums(exp(A - m)))
}

## -----------------------------------------------------------------------------
## [IT] .objnames: nomi degli oggetti (colonne). Se mancano usa "O1", "O2", ...
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.objnames <- function(X) {
  if (!is.null(colnames(X))) colnames(X) else paste0("O", seq_len(ncol(X)))
}

## -----------------------------------------------------------------------------
## [IT] .toorderings: trasforma una matrice di ranking (centri) nella matrice
## delle corrispondenti ORDINAZIONI scritte con i nomi degli oggetti, dal
## preferito al meno preferito (come "oclc" in cca()). Serve per la stampa.
## -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.toorderings <- function(R, labs) {
  R <- matrix(R, ncol = length(labs))
  O <- matrix(t(apply(R, 1, function(r) labs[order(r)])), nrow(R))
  dimnames(O) <- list(rownames(R), paste0("pos", seq_len(ncol(R))))
  O
}

## [IT] .distlabel: nome leggibile della distanza (per le stampe).
#' @keywords internal
#' @noRd
.distlabel <- function(dist) {
  c(kendall = "Kendall", cayley = "Cayley", hamming = "Hamming", ulam = "Ulam",
    footrule = "Footrule", spearman = "Spearman",
    spearmansq = "Spearman-squared")[[dist]]
}
