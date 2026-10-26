#compact drawing of a ranktree (plot.ranktree with style = "compact"): the tree
#grows from left to right, in black and white; the condition of each branch is
#written on the branch, the terminal nodes are aligned on the right with their
#class, size and Tau_x

#' @importFrom graphics par plot.new plot.window points segments strwidth text

plotcompacttree <- function(Tree, cex = 1){

  nd <- Tree$node
  num <- unlist(nd$number) #number of each node (root 1, children 2t and 2t + 1)
  term <- as.logical(unlist(nd$terminal)) #TRUE for the terminal nodes
  pos <- function(nn) which(num == nn) #position of a node in the lists of Tree$node
  if (length(num) == 1) {stop("The tree has only the root node: there is nothing to draw.")}

  #depth of each node (root = 0): number of halvings of its number down to 1
  depth <- sapply(num, function(nn) { d <- 0; while (nn > 1) { nn <- nn %/% 2; d <- d + 1 }; d })
  #vertical position: the terminal nodes from top to bottom in the order of the
  #tree (left child first), each internal node in the middle of its two children
  yy <- rep(NA_real_, length(num))
  nextleaf <- 0
  sety <- function(nn) {
    p <- pos(nn)
    if (term[p]) { nextleaf <<- nextleaf + 1; yy[p] <<- -nextleaf; return(yy[p]) } #terminal: next row
    yy[p] <<- mean(c(sety(2 * nn), sety(2 * nn + 1))) #internal: middle of the children
    yy[p]
  }
  invisible(sety(1))
  D <- max(depth) #depth of the deepest terminal node
  xx <- ifelse(term, D, depth) #terminal nodes aligned on the right

  #condition of the branch from the father to node nn
  branchlabel <- function(nn) {
    p <- pos(nn); pf <- pos(nn %/% 2) #the node and its father
    v <- nd$varsplit[[pf]] #split variable of the father
    cut <- nd$cutspli[[pf]] #cut point (numeric) or set of categories of the left child
    left <- nn %% 2 == 0 #left child = even number
    if (is.numeric(cut)) { #numeric predictor: the bound, as an expression (sign "less or equal")
      cv <- format(round(cut, 3))
      return(if (left) bquote(.(v) <= .(cv)) else bquote(.(v) > .(cv)))
    }
    here <- function(q) levels(droplevels(as.factor(Tree$X[[v]][nd$idatnode[[q]]]))) #categories in a node
    set <- here(p) #categories of the variable in this node
    sib <- here(pos(if (left) nn + 1 else nn - 1)) #categories in the sibling
    if (length(sib) < length(set)) { #fewer excluded categories: write those
      return(paste0(v, ": not ", paste(sib, collapse = ", ")))
    }
    paste0(v, ": ", paste(set, collapse = ", "))
  }

  #class of a terminal node as an ordering (bucket of tied objects in parentheses);
  #with a top-k consensus only the first k objects
  classlabel <- function(p) {
    y <- nd$nclass[[p]] #ranks of the objects (named)
    if (is.null(names(y))) {names(y) <- colnames(Tree$Y)}
    ranks <- sort(unique(y)) #positions of the buckets
    k <- Tree$control$topk #NULL if the consensus is not a top-k list
    out <- character(0); nobj <- 0
    for (r in ranks) {
      b <- names(y)[y == r] #objects of this bucket
      out <- c(out, if (length(b) > 1) paste0("(", paste(b, collapse = ", "), ")") else b)
      nobj <- nobj + length(b)
      if (!is.null(k) && nobj >= k) break #top-k: stop after the first k objects
    }
    paste(out, collapse = " > ")
  }

  #text of the terminal nodes and of the root
  lv <- which(term) #positions of the terminal nodes
  cls <- sapply(lv, classlabel) #their classes
  stat <- lapply(lv, function(p) #size and Tau_x
    bquote(m == .(format(nd$size[[p]], big.mark = ",")) * "," ~~ tau[X] == .(sprintf("%.3f", nd$tau[[p]]))))
  rootlab <- paste0("m = ", format(nd$size[[1]], big.mark = ","))

  #the frame: horizontal scale such that the labels of the root (on the left) and
  #of the terminal nodes (on the right) fit in the plotting region
  op <- par(mar = c(0.2, 0.2, 0.2, 0.2)); on.exit(par(op))
  plot.new()
  W <- par("pin")[1] #width of the plotting region in inches
  cl <- 0.75 * cex; cs <- 0.65 * cex; cb <- 0.6 * cex #text sizes: classes, statistics, branches
  wl <- max(strwidth(cls, units = "inches", cex = cl),
            sapply(stat, function(s) strwidth(s, units = "inches", cex = cs))) #widest label of the terminal nodes
  wr <- strwidth(rootlab, units = "inches", cex = cs) #label of the root
  pad <- 0.1 #inches between a node and its label
  u <- max((W - wl - wr - 3 * pad) / D, 0.2) #inches per level of the tree
  plot.window(xlim = c(-(wr + pad) / u, D + (wl + 2 * pad) / u), ylim = c(min(yy) - 0.5, -0.4),
              xaxs = "i")

  #the branches, as elbows from the father, with their conditions
  for (p in seq_along(num)[num > 1]) {
    pf <- pos(num[p] %/% 2) #the father
    segments(xx[pf], yy[pf], xx[pf], yy[p]) #vertical part at the father
    segments(xx[pf], yy[p], xx[p], yy[p]) #horizontal part to the node
    lab <- branchlabel(num[p])
    if (is.character(lab)) { #categories: wrapped to the length of the segment
      room <- 0.9 * (xx[p] - xx[pf]) * u #length of the segment in inches
      wd <- nchar(lab) + 1 #wrapping width in characters (strwrap breaks lines shorter than it), reduced until the lines fit
      repeat {
        lines <- strwrap(lab, width = wd)
        if (max(strwidth(lines, units = "inches", cex = cb)) <= room || wd <= 8) break
        wd <- wd - 1
      }
      lab <- paste(lines, collapse = "\n")
    }
    text((xx[pf] + xx[p]) / 2, yy[p], lab, pos = 3, offset = 0.2, cex = cb) #centred above the segment
  }

  #the nodes and their labels
  points(xx[!term], yy[!term], pch = 19, cex = 0.7) #internal nodes
  text(xx[1], yy[1], rootlab, pos = 2, cex = cs) #size of the root, on its left
  points(xx[lv], yy[lv], pch = 15, cex = 0.8) #terminal nodes
  dx <- pad / u #gap between a terminal node and its labels
  text(xx[lv] + dx, yy[lv] + 0.13, cls, adj = 0, cex = cl) #the class, above
  for (i in seq_along(lv)) { #size and Tau_x, below
    text(xx[lv[i]] + dx, yy[lv[i]] - 0.17, stat[[i]], adj = 0, cex = cs)
  }
  invisible(NULL)

}
