# ConsRankClass 1.1.0

## Contributors

* Luca Coraggio (contributor) is the author of the descriptive tools
  (`marginals()`, `pairmatrix()`, `rankproj()`, `objclust()`, `judgemds()`, with
  their `print` and `plot` methods) and, with Antonio D'Ambrosio, of `mdscale()`.

## New features

* `plot.ranktree()`: new option `style = "compact"` (default `"classic"`, the
  drawing used so far, unchanged), a black and white drawing of the tree from
  left to right, with the condition of each branch written on it and the
  terminal nodes aligned on the right with their class (with a top-k consensus
  only the first k objects), size and Tau_x; option `cex` for the size of its
  text.
* `validatetree()`: new argument `se` (default 0) for the s-SE rule. The
  standard error of the Tau_x of each subtree is returned (`tause`): with
  `unranked = "last"` it is the standard deviation of the Tau_x of the judges
  divided by sqrt(m); with `unranked = "unknown"` (Tau_x on the observed pairs, a
  ratio) it is the first-order standard error of a ratio estimator. The subtree
  chosen by the rule (`best_tause`) is the smallest one whose Tau_x is at least
  the largest Tau_x minus `se` times its standard error; with `se = 0` it is
  `best_tau`. With cross-validation the standard error is an approximation, as
  in Breiman et al. (1984). With `model = "mallows"` the same rule is applied to
  the deviance (`best_deviancese`; `best_deviance1se` is unchanged). All the
  previous components are unchanged.
* `plot.kunfolding()`: new option `bw = TRUE` for a black and white plot
  (individuals as grey points, groups given by `cols` as different symbols with
  a legend, items in black with their labels above them). The default plot is
  unchanged.
* `plot.rankproj()`: new options `labpos` (position of the names of the objects
  with respect to the tips of their poles: below, left, above, right; for all
  the objects or for some of them) and `laboffset` (shifts of single names), to
  separate overlapping names. The default plot is unchanged.

* Distances between complete rankings: `rankdist()`, `kendalld()`, `cayleyd()`,
  `hammingd()`, `ulamd()`, `footruled()`, `spearmand()`, and their exact
  distribution under the uniform model, `distfreq()`.
* Mallows (distance-based) model for all these distances, with exact
  normalising constants: `mallows()` (maximum likelihood estimation, with
  `print`, `summary` and `logLik` methods), `mallowscenter()`, `mallowsnorm()`,
  `mallowsexpect()`, `mallowslambda()`, `rmallows()` (exact simulation) and
  `dmallows()`.
* Mixtures of distance-based models (Murphy and Martin, 2003): `mixdbm()`,
  fitted by the EM algorithm with several starting points (with `print`,
  `summary`, `logLik` and `predict` methods), `mixdbmcontrol()`,
  `mixdbmselect()` (model selection by BIC and ICL) and `rmixdbm()`.
* `ranktree()`: new option `model` of `ranktreecontrol()`. `model = "kemeny"`
  (default) is the recursive partitioning method of D'Ambrosio and Heiser
  (2016), unchanged; `model = "mallows"` is the distance-based tree of Lee and
  Yu (2010): a Mallows model with Kendall distance is fitted within each node
  and the impurity is its deviance. Complete rankings and top-k partial
  rankings ("pick k out of M", e.g. the `EVS` data, which can be used directly
  with `NaN` for the unranked items) are allowed; the unranked objects are in
  unknown order and the likelihood of a top-k ranking, i.e. the sum over its
  compatible complete rankings, is computed exactly in closed form. Rankings
  with ties among the ranked objects are not allowed. Other new options of
  `ranktreecontrol()`: `impurity` (`model = "kemeny"` only; `"pairwise"`, the
  impurity used so far, is at present the only choice), `dist` (`"kendall"`),
  `minprop` (minimum size of a child node as a proportion of its father,
  default 0.1) and `lambdamax` (these three for `model = "mallows"` only). The
  options that do not concern the chosen model are ignored, with a warning if
  set explicitly.
* With `model = "mallows"`, `print()` and `summary()` report the measures of
  the Mallows models instead of those based on Kemeny distances: log-likelihood,
  deviance, number of parameters, AIC and BIC of the tree; modal ranking,
  precision lambda, phi = exp(-lambda), expected Kendall distance from the modal
  ranking and log-likelihood of each leaf (also in `nodepath()` and
  `treepaths()`). `validatetree()` returns the mean deviance of the subtrees on
  the test set or in cross-validation, its standard error, the best subtree and
  the subtree chosen by the 1-SE rule, and `plot(..., plot.type = "pruningseq")`
  shows the deviance along the pruning sequence.
* `ranktreecontrol(prunecost = "error")`: the cost-complexity pruning can use
  the Kemeny distance between the rankings and the median ranking of each node
  (the measure used by `validatetree()`) instead of the impurity (default,
  unchanged).
* `ranktreecontrol(algorithm = "BOLT")` (`model = "kemeny"`): the median
  ranking of the nodes can be computed with the BOLT heuristic of ConsRank
  (Romano, D'Ambrosio and Conversano), fast also with many objects; its options
  `otp` and `maxjudges` are new options of `ranktreecontrol()`. With
  `impurity = "pairwise"` the splits do not depend on the median rankings, so
  the tree is the same with any algorithm: only the classes of the nodes and
  the measures based on them (error, Tau_X, `prunecost = "error"`,
  `validatetree()`) change. The name of the algorithm is now checked, ignoring
  the case. The package requires ConsRank 3.1 or later.
* `ranktreecontrol(impurity = "median")` (`model = "kemeny"`): the impurity of
  a node is the sum of the Kemeny distances of its rankings from its median
  ranking (computed with `algorithm`), instead of the sum of the distances
  between all the pairs of rankings divided by the size of the node
  (`"pairwise"`, default, unchanged). The best split maximises the decrease of
  this sum divided by the sample size. The median rankings of the two children
  of every candidate split are computed from their combined input matrices,
  obtained from the counts of the pairwise preferences (cumulative sums over
  the cutting points of ordered predictors), and the impurity is computed
  exactly from the same counts; the class of each node is the median ranking
  used for its impurity. Rankings with missing values are not allowed. Each
  candidate split requires two median rankings, so the tree is much slower to
  grow than with `"pairwise"`, mainly with continuous predictors with many
  distinct values. For them the new option `maxcuts` of `ranktreecontrol()`
  (default `NULL`, all the cutting points) limits the cutting points of an
  ordered predictor in a node to those leaving in the left child the
  proportions 1/(K+1), ..., K/(K+1) of the judges (in a simulation with 300
  judges, 20 objects and a continuous predictor, from 68 s to 4-7 s with
  K = 10-20).
* Constrained consensus in `ranktree()` (`model = "kemeny"`): new options
  `nbuckets` (median bucket order with a fixed number of buckets),
  `bucketsizes` (median bucket order with a pattern of bucket sizes, as in
  `ConsRank::mcbo()`; requires `algorithm = "BOLT"`) and `topk` (median top-k
  list, with `ConsRank::consranktopk()`: exact, or BOLT with
  `algorithm = "BOLT"`) of `ranktreecontrol()`; at most one of them, and not
  with `full = TRUE`. The class of each node is the constrained consensus,
  computed from the combined input matrix of the node. With
  `impurity = "pairwise"` the tree is unchanged and only the classes (and the
  measures based on them) change; with `impurity = "median"` the impurity is
  the sum of the distances from the constrained consensus, so the constraint
  drives the splits. `print()` reports the constraint.
* Objects not ranked by the judges in `ranktree()` (`model = "kemeny"`): new
  option `unranked` of `ranktreecontrol()`. With `"last"` (default) the
  objects coded as `NA` are tied in the position after the last object ranked
  by the judge, as in D'Ambrosio and Heiser (2016) (before, rankings with `NA`
  gave an error; without `NA` the trees are unchanged). With `"unknown"` (only
  with `impurity = "median"`) they come after the ranked objects but in unknown
  order: a pair of two objects not ranked by a judge is not observed and counts
  neither as a preference nor as a tie. For each pair only the judges who rank
  at least one of its objects enter the combined input matrix of a node, and
  the impurity, `tau`, `error` (weighted by the observed pairs), `prunecost =
  "error"` and `validatetree()` are computed on the observed pairs only. For
  top-k lists the combined input matrix is the average of the matrices of all
  the completions of the lists. This convention is not in D'Ambrosio and
  Heiser (2016). `print()` reports it.

### Differences from Lee and Yu (2010) in `model = "mallows"`

* The center and the precision of each node are exact maximum likelihood
  estimates of the model for top-k rankings. With the same data (the `EVS`
  data, whose training set is the complement of their test set, Table 8) the
  log-likelihoods are larger than those of their Table 6 (Kendall, root node:
  -3275 against -3387), and the test deviance of the one-node model is smaller
  than that of the uniform model, which is not the case for their values: the
  published figures cannot be reproduced.
* Only the Kendall distance is available. Weighted distances, phi-component
  models and the other distances of their Table 6 are not.
* The tree is grown with the stopping rules of `ranktree()` (`num`,
  `decrmin`), not with their fixed minimum node size of 200 (which can be set
  with `num = 200`); `decrmin` is the decrease in deviance divided by the
  sample size. Their minimum child size of 1/10 of the father node is the
  default of `minprop`.
* The deviance returned by `validatetree()` is -2 times the mean
  log-likelihood per ranking.
* The number of parameters (for AIC and BIC) follows their Table 6: one
  precision for each leaf and one for each split; the modal rankings are not
  counted. BIC is not in their paper.

* Descriptive tools, each a function returning an object with `print` and
  `plot` methods (lattice graphics, grey levels): `marginals()`, the marginal
  matrix of the positions of the objects with their mean ranks (bubble plot),
  and `pairmatrix()`, the matrix of pairwise preferences (heatmap) with the
  counts of preferences and ties kept separate, the proportions with ties
  counting one half (default), of strict preferences or of ties, and the
  normalised Borda scores (Marden, 1995, Sections 2.3-2.4). Rankings with ties
  are allowed (a bucket shares its positions in the marginal matrix) and so are
  top-k partial rankings (unranked objects tied in the last position), with
  optional weights. The package now imports lattice, grid and grDevices.
* `rankproj()`: projection of the rankings on the principal planes of the
  permutation polytope (the biplot of Cohen and Mallows, 1980; Marden, 1995,
  Section 2.2), i.e. the eigenvectors of the second moments of the rankings
  about the center of the polytope (not about the mean ranking), with the
  poles of the objects. The scores are scaled so that permutations lie on the
  unit sphere (rankings with ties inside), and the poles are projected unit
  vectors, whose length shows how well each object is represented. Ties
  (midranks), top-k rankings and weights are allowed. `plot()` draws the
  projected polytope-o-gram (symbols with area proportional to the number of
  judges, poles as arrows, unit circle), optionally by groups of judges (e.g.
  leaves of a tree or clusters) and with the labels of the rankings.
* `objclust()`: hierarchical clustering of the objects (Marden, 1995, Example
  2.18), with the mean squared difference of their ranks as distance and
  complete linkage by default; `plot()` draws the dendrogram, optionally with
  the groups of a cut. On the sports data the sequence of merges is the one
  of Marden (1995, Figure 2.10).
* `judgemds()`: multidimensional scaling of the judges with `mdscale()` on the
  Kemeny distances between their distinct rankings, each pair weighted by the
  product of the frequencies of the two rankings (`weighted = TRUE`), so that
  the configuration represents the judges; `plot()` draws the distinct
  rankings with area proportional to their frequency, optionally by groups and
  with labels.
* `mdscale()`: metric and non-metric multidimensional scaling with the six
  criteria of the MATLAB function `mdscale` (`"sstress"`, the default;
  `"stress"`, i.e. Kruskal's stress-1; `"metricstress"`; `"metricsstress"`;
  `"sammon"`; `"strain"`), weights, missing dissimilarities (`NA`), similarity
  matrices, starts `"cmdscale"`, `"random"` or user (matrix or array) and
  replicates, with `print`, `summary` and `plot` (configuration or Shepard
  diagram) methods; options in `mdscalecontrol()`. The engine is written in C++
  (the package now has compiled code and imports Rcpp). On the same data and from
  the same start the value of each criterion is the one computed by MATLAB
  (R2023a) to about 1e-16. In a comparison with MATLAB on the same random
  starts (20-30 per data set), the mean of the final criterion was equal or
  lower in almost all cases (differences are due to different local minima),
  and `mdscale()` was 50-200 times faster.
* `kunfolding()`: new option `engine` of `mdscontrol()`: `"smacof"` (default,
  unchanged) or `"mdscale"`; with `"mdscale"` the options `criterion` and `nrep`
  (replicates) of `mdscontrol()` can be used. With `engine = "smacof"` the
  results are identical to those of version 1.0.2. The output has two new
  components, `engine` and `criterion`.

### Differences from the MATLAB function `mdscale`

* The default criterion is `"sstress"`; the default of MATLAB is `"stress"`
  (the kind of criterion already minimised by the default engine of
  `kunfolding()`, `smacofSym`). With `engine = "mdscale"`, `kunfolding()` with
  `transf = "primary"` or `"secondary"` minimises `"sstress"` too, and with
  `transf = "ratio"` or `"none"` it minimises `"metricsstress"`.

* The criteria are minimised with a limited-memory BFGS method with a line
  search for the strong Wolfe conditions, on the criterion as a function of the
  configuration only (the disparities are recomputed at every evaluation and the
  gradient includes their variation, which is zero for `"stress"`). MATLAB uses
  a different gradient method, so the two functions can stop at different local
  minima.
* For the non-metric criteria, which do not depend on the scale, the
  configuration is brought back to its initial size after each iteration, and
  the final configuration is scaled so that the weighted sum of squared
  distances equals the weighted sum of squared dissimilarities.
* Stopping rules: the relative change of the criterion is measured over the
  last `window` iterations (default 10; `window = 1` is the rule of MATLAB);
  the norm of the gradient is multiplied by the norm of the configuration and
  divided by the criterion, and the norm of the step is divided by the norm of
  the configuration, so that the rules do not depend on the unit of the
  dissimilarities (with Kemeny distances the absolute rules of MATLAB stop too
  early). The default values of `maxiter`, `tolfun` and `tolx` are those of
  MATLAB.
* Option `ties = "secondary"` (tied dissimilarities get the same disparity);
  the default `"primary"` is the treatment of MATLAB.
* Coincident points (zero distances, e.g. identical rankings with Kemeny
  distances) are allowed: MATLAB stops with an error ("Points in the
  configuration have co-located").
* The sign of each axis of the final configuration is chosen so that its
  largest absolute coordinate is positive.
* The random start is normal with variance such that the expected squared
  distance equals the mean squared dissimilarity, drawn with the random number
  generator of R.

## Speed

* `ranktree()` with `model = "kemeny"` computes the impurity of all the
  candidate splits from the pairwise preferences of the judges (cumulative
  sums), without the matrices of Kemeny distances: the trees are identical to
  those of version 1.0.2 and are built up to about 60 times faster (cross-validation
  with `validatetree()` too). Rankings with missing values use the previous
  code.
* `fuzzyconcordance()`: the expected value of the normalized degree of
  concordance under permutation is computed exactly in O(N log N) time
  (`nperms = NULL`, new default), instead of by Monte Carlo (still available
  with `nperms` set to a number).

## Bug fixes

* `fuzzyconcordance()`: `as-matrix` replaced by `as.matrix`.
* `ranktree()` failed when no split of the root node was accepted (tree with a
  single node), and so did `summary()`, `nodepath()` and `treepaths()` on such
  a tree (also when obtained with `getsubtree()`). `validatetree()` now stops
  with a clear message for such a tree.
* `mdscontrol()`: `userinit` was accepted but not returned, so
  `kunfolding()` with `init = "user"` failed; now the starting configuration is
  used (by both engines).
* `kunfolding()`: the component `avgspearman` was the mean Pearson correlation
  (the same value as `avgedpearson`); it is now the mean Spearman correlation,
  as documented (the value printed in the recovery measures was already right).
* `kunfolding()`: `transf = "tertiary"` was silently fitted with the primary
  approach to ties; now the tertiary approach of `smacofSym` is used. With
  `engine = "mdscale"` it stops with a message (not available).
* `kunfolding()`: `model = "metric"` with `transf = "primary"` (the default) or
  `"secondary"` stopped with the error "object 'stype' not found"; now the
  transformation wins, as the code intended, and an ordinal MDS is fitted
  (with `engine = "mdscale"` too, criterion `"sstress"`).
* `kunfolding()`: the messages were printed also with `printscr = FALSE`; now
  nothing is printed.
* `kunfolding()`: with a judge who ties all the items the Pearson, Spearman and
  Kendall correlations of that judge are not defined, and the averaged
  correlations (`avgedpearson`, `avgspearman`, `avgkendall` and the recovery
  measures) were `NA`, with warnings from `cor()`. Now these averages are
  computed on the other judges, with one warning saying how many judges were
  excluded; such judges stay in the solution and in the averaged Tau_X. Without
  such judges the results do not change.
* `predict.ranktree()`: the column `leafid` of `info` was the position of the
  leaf in the internal list of the nodes (e.g. 10 for the node 14); it is now the
  node number, as in `print()`, `summary()` and `nodepath()`. The predicted
  rankings and orderings do not change.
* `cca()`: the result was printed inside the function, also when assigned;
  now it is printed only when shown (`print.cca()`). The result does not change.
* `summary.ranktree()`: `NULL` was printed at the end of the summary; now the
  paths of the leaves are returned invisibly.
* Typos: "itrems" in `print.kunfolding()`; the row "Sperman" (now "Spearman") of
  the component `resumeaug` of `kunfolding()` and `res` of `augmatrix()`.
* `plot.kunfolding()`: `cols` given as a factor gave the warning "'+' not
  meaningful for factors" and no colours; now the groups of a factor (or of a
  character vector) are used. Numeric `cols` as before.
* Help of `EVS`: 1911 respondents (not 3584), the `NaN` of the rankings and
  their treatment in `ranktree()`, country names.

# ConsRankClass 1.0.2

* Version on CRAN before the changes above.
