// mdscale_cpp.cpp
// Engine of mdscale(): metric and non-metric multidimensional scaling with the six
// criteria of MATLAB's mdscale (stress, sstress, metricstress, metricsstress,
// sammon, strain). The criteria are written from their definitions (Kruskal, 1964;
// Takane, Young and de Leeuw, 1977; Sammon, 1969; Torgerson, 1958) and are minimised
// as functions of the configuration only (for the non-metric criteria the disparities
// are recomputed at every evaluation), with a limited-memory BFGS method and a line
// search for the strong Wolfe conditions (Nocedal and Wright, 2006).
// Author: Antonio D'Ambrosio

#include <Rcpp.h>          // R/C++ interface
#include <vector>          // std::vector
#include <algorithm>       // std::sort, std::stable_sort
#include <cmath>           // std::sqrt, std::fabs
#include <numeric>         // std::iota
using namespace Rcpp;

// Codes of the criteria (the same order is used on the R side)
// 1 = stress, 2 = sstress, 3 = metricstress, 4 = metricsstress, 5 = sammon, 6 = strain

// ---------------------------------------------------------------------------
// Weighted pool-adjacent-violators algorithm.
// y: values in the order in which they must be non-decreasing; w: weights.
// The least squares monotone fit is written in fit; blk receives, for each
// position, the index of the block (set of pooled values) it belongs to.
// ---------------------------------------------------------------------------
static void pava(const std::vector<double>& y, const std::vector<double>& w,
                 std::vector<double>& fit, std::vector<int>& blk) {
  int n = y.size();                                  // number of values
  std::vector<double> val(n), wt(n);                 // value and weight of each block
  std::vector<int> len(n);                           // number of values in each block
  int nb = 0;                                        // current number of blocks
  for (int i = 0; i < n; i++) {                      // scan the values from left to right
    val[nb] = y[i]; wt[nb] = w[i]; len[nb] = 1;      // the new value is a new block
    nb++;                                            // one more block
    while (nb > 1 && val[nb - 2] > val[nb - 1]) {    // while the last two blocks violate the order
      double ws = wt[nb - 2] + wt[nb - 1];           // weight of the pooled block
      val[nb - 2] = (wt[nb - 2] * val[nb - 2] + wt[nb - 1] * val[nb - 1]) / ws; // weighted mean
      wt[nb - 2] = ws;                               // weight of the pooled block
      len[nb - 2] += len[nb - 1];                    // size of the pooled block
      nb--;                                          // the two blocks are now one
    }
  }
  int k = 0;                                         // position in the output
  for (int b = 0; b < nb; b++)                       // for each block
    for (int l = 0; l < len[b]; l++) { fit[k] = val[b]; blk[k] = b; k++; } // value and block index
}

// ---------------------------------------------------------------------------
// Data shared by the evaluations of the criterion
// ---------------------------------------------------------------------------
struct MdsData {
  int N, p, crit, ties;                 // points, dimensions, criterion, ties (1 primary, 2 secondary)
  int K;                                // number of active pairs (non-missing, positive weight)
  std::vector<int> ii, jj;              // the two points of each active pair (0-based)
  std::vector<double> dis, w;           // dissimilarity and weight of each active pair
  std::vector<int> ord;                 // active pairs sorted by increasing dissimilarity
  std::vector<int> bstart;              // start of each block of tied dissimilarities in ord (plus K at the end)
  double sumwD2, sumwD4, sumwD;         // constants of the metric criteria
  NumericMatrix B2, W1;                 // strain only: -D^2/2 and weights with unit diagonal
  std::vector<double> rw;               // strain only: row sums of the weights (diagonal excluded)
  std::vector<double> dh;               // disparities of the active pairs (last evaluation)
  std::vector<int> blk;                 // monotone regression block of each active pair (last evaluation)
  int nblk;                             // number of monotone regression blocks
};

// ---------------------------------------------------------------------------
// Monotone regression of the distances on the order of the dissimilarities
// (non-metric criteria). The disparities are written in M.dh and the block of each
// pair in M.blk. Primary approach: tied dissimilarities may get different disparities
// (within a tie block the pairs are put in increasing order of distance). Secondary
// approach: tied dissimilarities get the same disparity (the tie block is pooled first).
// ---------------------------------------------------------------------------
static void monotone(MdsData& M, const std::vector<double>& d) {
  int K = M.K;                                         // number of active pairs
  int ntie = M.bstart.size() - 1;                      // number of tie blocks
  if (M.ties == 1) {                                   // primary approach
    std::vector<int> o(M.ord);                         // working copy of the order
    for (int b = 0; b < ntie; b++) {                   // for each tie block
      int s = M.bstart[b], e = M.bstart[b + 1];        // its range in the order
      if (e - s > 1)                                   // only blocks with ties need sorting
        std::sort(o.begin() + s, o.begin() + e,        // sort the block by distance
                  [&d](int a, int c) { return d[a] < d[c]; });
    }
    std::vector<double> y(K), wv(K), f(K);             // values, weights and fit in that order
    std::vector<int> bb(K);                            // block of each position
    for (int k = 0; k < K; k++) { y[k] = d[o[k]]; wv[k] = M.w[o[k]]; } // fill them
    pava(y, wv, f, bb);                                // monotone regression
    for (int k = 0; k < K; k++) { M.dh[o[k]] = f[k]; M.blk[o[k]] = bb[k]; } // back to the pair index
    M.nblk = (K > 0) ? bb[K - 1] + 1 : 0;              // number of blocks
  } else {                                             // secondary approach
    std::vector<double> y(ntie), wv(ntie), f(ntie);    // one value per tie block
    std::vector<int> bb(ntie);                         // block of each tie block
    for (int b = 0; b < ntie; b++) {                   // for each tie block
      double sw = 0.0, sy = 0.0;                       // sum of weights and of weighted distances
      for (int k = M.bstart[b]; k < M.bstart[b + 1]; k++) { // members of the tie block
        sw += M.w[M.ord[k]]; sy += M.w[M.ord[k]] * d[M.ord[k]]; // accumulate
      }
      y[b] = (sw > 0) ? sy / sw : 0.0;                 // weighted mean of the tie block
      wv[b] = sw;                                      // weight of the tie block
    }
    pava(y, wv, f, bb);                                // monotone regression on the tie blocks
    for (int b = 0; b < ntie; b++)                     // for each tie block
      for (int k = M.bstart[b]; k < M.bstart[b + 1]; k++) { // its members share
        M.dh[M.ord[k]] = f[b]; M.blk[M.ord[k]] = bb[b];     // disparity and block
      }
    M.nblk = (ntie > 0) ? bb[ntie - 1] + 1 : 0;        // number of blocks
  }
}

// ---------------------------------------------------------------------------
// Criterion and gradient at the configuration X (N x p, column-major in x).
// Returns the value of the criterion and writes its gradient in g.
// ---------------------------------------------------------------------------
static double evalcrit(MdsData& M, const std::vector<double>& x, std::vector<double>& g) {
  int N = M.N, p = M.p, K = M.K;                       // sizes
  std::fill(g.begin(), g.end(), 0.0);                  // gradient starts at zero
  if (M.crit == 6) {                                   // ---- strain
    double sr = 0.0;                                   // total of the row sums
    for (int i = 0; i < N; i++) sr += M.rw[i];         // accumulate
    std::vector<double> mu(p, 0.0);                    // weighted mean of each column
    for (int a = 0; a < p; a++) {                      // for each dimension
      for (int i = 0; i < N; i++) mu[a] += M.rw[i] * x[i + a * N]; // weighted sum
      mu[a] /= sr;                                     // weighted mean
    }
    std::vector<double> y(N * p);                      // centred configuration
    for (int a = 0; a < p; a++)                        // for each dimension
      for (int i = 0; i < N; i++) y[i + a * N] = x[i + a * N] - mu[a]; // centre
    double F = 0.0;                                    // sum of weighted squared residuals
    std::vector<double> EY(N * p, 0.0);                // E %*% Y, E = W1 * (YY' - B2)
    for (int i = 0; i < N; i++) {                      // rows
      for (int j = 0; j < N; j++) {                    // columns
        if (M.W1(i, j) == 0.0) continue;               // zero weight: no contribution
        double gij = 0.0;                              // inner product of points i and j
        for (int a = 0; a < p; a++) gij += y[i + a * N] * y[j + a * N]; // accumulate
        double r = gij - M.B2(i, j);                   // residual
        F += M.W1(i, j) * r * r;                       // weighted squared residual
        double e = M.W1(i, j) * r;                     // element of E
        for (int a = 0; a < p; a++) EY[i + a * N] += e * y[j + a * N]; // E %*% Y
      }
    }
    double S = std::sqrt(F);                           // criterion
    if (S <= 0) return 0.0;                            // perfect fit: zero gradient
    std::vector<double> csum(p, 0.0);                  // column sums of E %*% Y
    for (int a = 0; a < p; a++)                        // for each dimension
      for (int i = 0; i < N; i++) csum[a] += EY[i + a * N]; // accumulate
    for (int a = 0; a < p; a++)                        // gradient with respect to X:
      for (int i = 0; i < N; i++)                      // C'(4 E Y) / (2 S), C = I - 1 r'/sum(r)
        g[i + a * N] = (4.0 * (EY[i + a * N] - M.rw[i] * csum[a] / sr)) / (2.0 * S);
    return S;                                          // value of the criterion
  }
  std::vector<double> d(K);                            // distances of the active pairs
  for (int k = 0; k < K; k++) {                        // for each active pair
    double s2 = 0.0;                                   // squared distance
    for (int a = 0; a < p; a++) {                      // for each dimension
      double df = x[M.ii[k] + a * N] - x[M.jj[k] + a * N]; // coordinate difference
      s2 += df * df;                                   // accumulate
    }
    d[k] = std::sqrt(s2);                              // distance
  }
  std::vector<double> coef(K, 0.0);                    // derivative of the criterion w.r.t. d (or d^2)
  bool sq = (M.crit == 4);                             // metric sstress: derivative w.r.t. d^2
  double S = 0.0;                                      // value of the criterion
  if (M.crit == 1) {                                   // ---- non-metric stress (Kruskal's stress-1)
    monotone(M, d);                                    // disparities
    double T = 0.0, U = 0.0;                           // raw stress and normalising sum
    for (int k = 0; k < K; k++) {                      // for each pair
      T += M.w[k] * (M.dh[k] - d[k]) * (M.dh[k] - d[k]); // weighted squared residual
      U += M.w[k] * d[k] * d[k];                       // weighted squared distance
    }
    if (U <= 0) return 0.0;                            // all points coincide
    S = std::sqrt(T / U);                              // stress-1
    if (S <= 0) return 0.0;                            // perfect fit
    // the disparities minimise the raw stress for the given distances, so their own
    // variation does not contribute to the gradient (within each block the residuals sum to 0)
    for (int k = 0; k < K; k++)                        // dS/dd
      coef[k] = M.w[k] * ((d[k] - M.dh[k]) - (T / U) * d[k]) / (S * U);
  } else if (M.crit == 2) {                            // ---- non-metric sstress
    // as in MATLAB, the disparities are the monotone regression of the distances
    // (the same as for stress) and the criterion compares their squares with d^2
    monotone(M, d);                                    // disparities
    double T = 0.0, U = 0.0;                           // raw sstress and normalising sum
    for (int k = 0; k < K; k++) {                      // for each pair
      double r = M.dh[k] * M.dh[k] - d[k] * d[k];      // residual on squares
      T += M.w[k] * r * r;                             // weighted squared residual
      U += M.w[k] * std::pow(d[k], 4);                 // weighted fourth power of distances
    }
    if (U <= 0) return 0.0;                            // all points coincide
    S = std::sqrt(T / U);                              // sstress
    if (S <= 0) return 0.0;                            // perfect fit
    std::vector<double> A(M.nblk, 0.0), WB(M.nblk, 0.0); // per block: sum of dS/ddh and of weights
    for (int k = 0; k < K; k++) {                      // for each pair
      double r = M.dh[k] * M.dh[k] - d[k] * d[k];      // residual on squares
      A[M.blk[k]] += 2.0 * M.w[k] * r * M.dh[k] / (S * U); // dS/ddh_k accumulated in its block
      WB[M.blk[k]] += M.w[k];                          // weight of the block
      // direct part: dS/dd_k = dS/d(d_k^2) * 2 d_k
      coef[k] = 2.0 * d[k] * M.w[k] * (-r - (T / U) * d[k] * d[k]) / (S * U);
    }
    // indirect part: each disparity is the weighted mean of the distances of its block,
    // so ddh_k/dd_l = w_l / W_B for k, l in the same block B
    for (int k = 0; k < K; k++) coef[k] += M.w[k] / WB[M.blk[k]] * A[M.blk[k]];
  } else if (M.crit == 3) {                            // ---- metric stress
    double T = 0.0;                                    // raw stress
    for (int k = 0; k < K; k++) T += M.w[k] * (M.dis[k] - d[k]) * (M.dis[k] - d[k]); // accumulate
    S = std::sqrt(T / M.sumwD2);                       // normalised by the dissimilarities
    if (S <= 0) return 0.0;                            // perfect fit
    for (int k = 0; k < K; k++) coef[k] = M.w[k] * (d[k] - M.dis[k]) / (S * M.sumwD2); // dS/dd
  } else if (M.crit == 4) {                            // ---- metric sstress
    double T = 0.0;                                    // raw sstress
    for (int k = 0; k < K; k++) {                      // for each pair
      double r = d[k] * d[k] - M.dis[k] * M.dis[k];    // residual on squares
      T += M.w[k] * r * r;                             // accumulate
    }
    S = std::sqrt(T / M.sumwD4);                       // normalised by the dissimilarities
    if (S <= 0) return 0.0;                            // perfect fit
    for (int k = 0; k < K; k++)                        // dS/d(d^2)
      coef[k] = M.w[k] * (d[k] * d[k] - M.dis[k] * M.dis[k]) / (S * M.sumwD4);
  } else if (M.crit == 5) {                            // ---- Sammon's mapping
    double T = 0.0;                                    // weighted relative squared residuals
    for (int k = 0; k < K; k++) T += M.w[k] * (M.dis[k] - d[k]) * (M.dis[k] - d[k]) / M.dis[k]; // accumulate
    S = T / M.sumwD;                                   // Sammon's criterion
    for (int k = 0; k < K; k++) coef[k] = 2.0 * M.w[k] * (d[k] - M.dis[k]) / (M.dis[k] * M.sumwD); // dS/dd
  }
  for (int k = 0; k < K; k++) {                        // chain rule to the coordinates
    int i = M.ii[k], j = M.jj[k];                      // the two points
    double c;                                          // factor multiplying x_i - x_j
    if (sq) c = 2.0 * coef[k];                         // d(d^2)/dx_i = 2 (x_i - x_j)
    else {                                             // dd/dx_i = (x_i - x_j)/d
      if (d[k] <= 0) continue;                         // coincident points: no direction
      c = coef[k] / d[k];                              // factor
    }
    for (int a = 0; a < p; a++) {                      // for each dimension
      double df = x[i + a * N] - x[j + a * N];         // coordinate difference
      g[i + a * N] += c * df;                          // point i
      g[j + a * N] -= c * df;                          // point j, opposite sign
    }
  }
  return S;                                            // value of the criterion
}

// inner product of two vectors
static double dot(const std::vector<double>& a, const std::vector<double>& b) {
  double s = 0.0;                                      // inner product
  for (size_t i = 0; i < a.size(); i++) s += a[i] * b[i]; // accumulate
  return s;                                            // result
}

// ---------------------------------------------------------------------------
// Main engine: minimises the criterion from the starting configuration X0.
// dis, w, ii, jj: dissimilarities, weights and points (0-based) of the active pairs;
// B2, W1, rw: strain only; maxiter, tolfun, tolx: stopping rules; display: 0 off,
// 1 final, 2 iterations; mhist: number of pairs kept in the L-BFGS memory.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List mdscale_engine(NumericVector dis, NumericVector w, IntegerVector ii, IntegerVector jj,
                    int N, NumericMatrix X0, int crit, int ties,
                    NumericMatrix B2, NumericMatrix W1, NumericVector rw,
                    int maxiter, double tolfun, double tolx, int display, int mhist = 10, int relwin = 1) {
  MdsData M;                                            // shared data
  M.N = N; M.p = X0.ncol(); M.crit = crit; M.ties = ties; // sizes and options
  M.K = dis.size();                                     // number of active pairs
  M.ii = as<std::vector<int> >(ii);                     // first point of each pair
  M.jj = as<std::vector<int> >(jj);                     // second point of each pair
  M.dis = as<std::vector<double> >(dis);                // dissimilarities
  M.w = as<std::vector<double> >(w);                    // weights
  M.dh.assign(M.K, 0.0);                                // disparities
  M.blk.assign(M.K, 0); M.nblk = 0;                     // monotone regression blocks
  M.sumwD2 = 0.0; M.sumwD4 = 0.0; M.sumwD = 0.0;        // constants of the metric criteria
  for (int k = 0; k < M.K; k++) {                       // for each pair
    M.sumwD2 += M.w[k] * M.dis[k] * M.dis[k];           // sum of w D^2
    M.sumwD4 += M.w[k] * std::pow(M.dis[k], 4);         // sum of w D^4
    M.sumwD += M.w[k] * M.dis[k];                       // sum of w D
  }
  M.ord.resize(M.K);                                    // order of the pairs
  std::iota(M.ord.begin(), M.ord.end(), 0);             // 0, 1, ..., K-1
  std::stable_sort(M.ord.begin(), M.ord.end(),          // sorted by dissimilarity
                   [&M](int a, int b) { return M.dis[a] < M.dis[b]; });
  M.bstart.clear();                                     // starts of the tie blocks
  for (int k = 0; k < M.K; k++)                         // a block starts where the value changes
    if (k == 0 || M.dis[M.ord[k]] != M.dis[M.ord[k - 1]]) M.bstart.push_back(k);
  M.bstart.push_back(M.K);                              // end of the last block
  M.B2 = B2; M.W1 = W1;                                 // strain matrices
  M.rw = as<std::vector<double> >(rw);                  // strain row sums

  int n = N * M.p;                                      // number of free coordinates
  bool nonmetric = (crit == 1 || crit == 2);            // scale-invariant criteria
  std::vector<double> x(X0.begin(), X0.end());          // current configuration
  std::vector<double> g(n), xt(n), gt(n), dvec(n);      // gradient, trial point, trial gradient, direction
  std::vector<double> xb(n), gb(n);                     // accepted trial point and its gradient
  double f = evalcrit(M, x, g);                         // criterion at the start
  std::vector<std::vector<double> > S, Y;               // L-BFGS memory: steps and gradient changes
  std::vector<double> rho;                              // 1 / (y's) for each stored pair
  int iter = 0;                                         // iteration counter
  std::vector<double> fhist(1, f);                      // criterion at the start and after each iteration
  int reason = 0;                                       // 0 = iteration limit
  double gnorm = std::sqrt(dot(g, g));                  // norm of the gradient
  double xnorm0 = std::sqrt(dot(x, x));                 // size of the starting configuration
  if (display == 2) {                                   // iterative display, header
    Rprintf("\n   Iteration       Criterion   Norm of gradient    Norm of step   Line search\n");
    Rprintf("  %10d  %14.6g  %16.6g\n", 0, f, gnorm);   // starting values
  }
  if (f < tolfun) reason = 1;                           // already below the tolerance
  else if (gnorm * xnorm0 < tolfun * f) reason = 3;    // already stationary (scale-free gradient)

  // value and directional derivative of the criterion along dvec at step a;
  // the trial point and its gradient are left in xt and gt
  auto phi = [&](double a, double& dphi) {
    for (int t = 0; t < n; t++) xt[t] = x[t] + a * dvec[t]; // trial point
    double v = evalcrit(M, xt, gt);                     // criterion and gradient there
    dphi = dot(gt, dvec);                               // directional derivative
    return v;                                           // value
  };

  while (reason == 0 && iter < maxiter) {               // main loop
    iter++;                                             // new iteration
    // --- search direction: two-loop recursion of L-BFGS
    std::vector<double> q(g);                           // start from the gradient
    int m = S.size();                                   // stored pairs
    std::vector<double> alpha(m);                       // coefficients of the first loop
    for (int l = m - 1; l >= 0; l--) {                  // newest to oldest
      alpha[l] = rho[l] * dot(S[l], q);                 // coefficient
      for (int t = 0; t < n; t++) q[t] -= alpha[l] * Y[l][t]; // update q
    }
    double gamma = 1.0;                                 // scaling of the initial Hessian
    if (m > 0) gamma = dot(S[m - 1], Y[m - 1]) / dot(Y[m - 1], Y[m - 1]); // usual L-BFGS scaling
    for (int t = 0; t < n; t++) q[t] *= gamma;          // apply it
    for (int l = 0; l < m; l++) {                       // oldest to newest
      double beta = rho[l] * dot(Y[l], q);              // coefficient
      for (int t = 0; t < n; t++) q[t] += S[l][t] * (alpha[l] - beta); // update q
    }
    for (int t = 0; t < n; t++) dvec[t] = -q[t];        // direction = -H g
    double gd = dot(g, dvec);                           // directional derivative at the current point
    if (!(gd < 0)) {                                    // not a descent direction
      S.clear(); Y.clear(); rho.clear(); m = 0;         // forget the memory
      for (int t = 0; t < n; t++) dvec[t] = -g[t];      // steepest descent
      gd = -dot(g, g);                                  // its directional derivative
    }
    double a = 1.0;                                     // first trial step (the natural L-BFGS step)
    if (m == 0) {                                       // steepest descent: scale the first step
      double xnorm = std::sqrt(dot(x, x));              // size of the current configuration
      double scale = (xnorm > 0) ? xnorm : ((xnorm0 > 0) ? xnorm0 : 1.0); // fall back on the start
      a = 0.1 * scale / std::sqrt(dot(dvec, dvec));     // move by 10% of the configuration size
    }
    // --- line search for the strong Wolfe conditions (Nocedal and Wright, 2006, Alg. 3.5-3.6)
    const double c1 = 1e-4, c2 = 0.9;                   // sufficient decrease and curvature constants
    double f0 = f;                                      // value at the current point
    double aprev = 0.0, fprev = f0, dprev = gd;         // previous trial step
    double astar = 0.0, fstar = f0;                     // accepted step and its value
    int ls = 0;                                         // number of evaluations
    double alo = 0, ahi = 0, flo = 0, dlo = 0, fhi = 0; // bracket for the zoom phase
    bool zoom = false, done = false;                    // whether the zoom phase is needed / finished
    while (ls < 30) {                                   // bracketing phase
      double da; double fa = phi(a, da); ls++;          // evaluate the trial step
      if (fa > f0 + c1 * a * gd || (ls > 1 && fa >= fprev)) { // too long: bracket found
        alo = aprev; flo = fprev; dlo = dprev; ahi = a; fhi = fa; zoom = true; break;
      }
      xb = xt; gb = gt; astar = a; fstar = fa;          // sufficient decrease: keep this step
      if (std::fabs(da) <= -c2 * gd) { done = true; break; } // strong Wolfe conditions hold
      if (da >= 0) {                                    // past a minimum: bracket found
        alo = a; flo = fa; dlo = da; ahi = aprev; fhi = fprev; zoom = true; break;
      }
      aprev = a; fprev = fa; dprev = da;                // remember it
      a *= 2.0;                                         // and try a longer step
    }
    if (zoom && !done) {                                // zoom phase between alo and ahi
      while (ls < 30) {
        double dlt = ahi - alo;                         // signed width of the bracket
        double den = 2.0 * (fhi - flo - dlo * dlt);     // curvature of the interpolating quadratic
        double aj = (den > 0) ? alo - dlo * dlt * dlt / den : alo + 0.5 * dlt; // its minimiser
        double lo = std::min(alo, ahi), hi = std::max(alo, ahi); // interval
        if (!(aj > lo + 0.1 * (hi - lo) && aj < hi - 0.1 * (hi - lo))) aj = 0.5 * (alo + ahi); // safeguard: bisection
        double dj; double fj = phi(aj, dj); ls++;       // evaluate
        if (fj > f0 + c1 * aj * gd || fj >= flo) {      // not enough decrease: shrink from above
          ahi = aj; fhi = fj;
        } else {
          xb = xt; gb = gt; astar = aj; fstar = fj;     // sufficient decrease: keep it
          if (std::fabs(dj) <= -c2 * gd) break;         // strong Wolfe conditions hold
          if (dj * (ahi - alo) >= 0) { ahi = alo; fhi = flo; } // keep a minimum inside the bracket
          alo = aj; flo = fj; dlo = dj;                 // move the lower end
        }
        if (std::fabs(ahi - alo) < 1e-14 * std::max(1.0, std::fabs(alo))) break; // bracket collapsed
      }
    }
    if (astar <= 0.0) {                                 // no step with sufficient decrease
      if (m > 0) {                                      // with memory: retry with steepest descent
        S.clear(); Y.clear(); rho.clear();              // forget the memory
        iter--;                                         // this iteration does not count
        continue;                                       // new attempt
      }
      reason = 4;                                       // the step cannot be reduced further
      break;                                            // stop
    }
    std::vector<double> s(n), y(n);                     // step and gradient change
    for (int t = 0; t < n; t++) { s[t] = xb[t] - x[t]; y[t] = gb[t] - g[t]; } // compute them

    x = xb; g = gb; f = fstar;                          // accept the step
    double snorm = std::sqrt(dot(s, s));                // norm of the step
    if (nonmetric && xnorm0 > 0) {                      // non-metric criteria are scale invariant:
      double xn = std::sqrt(dot(x, x));                 // the configuration could slowly shrink or grow,
      double c = xnorm0 / xn;                           // so it is brought back to the size of the start
      for (int t = 0; t < n; t++) { x[t] *= c; g[t] /= c; s[t] *= c; y[t] /= c; } // f(cX) = f(X), grad(cX) = grad(X)/c
      for (size_t l = 0; l < S.size(); l++)             // the stored pairs are mapped to the new scale
        for (int t = 0; t < n; t++) { S[l][t] *= c; Y[l][t] /= c; } // (s'y, hence rho, is unchanged)
    }
    double sy = dot(s, y);                              // curvature
    if (sy > 1e-12 * std::sqrt(dot(s, s) * dot(y, y))) { // keep the pair only if the curvature is positive
      if ((int) S.size() == mhist) {                    // memory full
        S.erase(S.begin()); Y.erase(Y.begin()); rho.erase(rho.begin()); // drop the oldest
      }
      S.push_back(s); Y.push_back(y); rho.push_back(1.0 / sy); // store the new pair
    }
    gnorm = std::sqrt(dot(g, g));                       // norm of the new gradient
    if (display == 2) Rprintf("  %10d  %14.6g  %16.6g  %14.6g  %10d\n", iter, f, gnorm, snorm, ls); // display
    fhist.push_back(f);                                 // store the new value of the criterion
    int nf = fhist.size();                              // values stored (start included)
    double fref = (nf - 1 >= relwin) ? fhist[nf - 1 - relwin] : -1.0; // value relwin iterations ago
    if (f < tolfun) reason = 1;                         // criterion below tolerance
    else if (fref > 0 && std::fabs(fref - f) < tolfun * fref) reason = 2; // small relative change over relwin iterations
    else if (gnorm * std::sqrt(dot(x, x)) < tolfun * f) reason = 3; // small gradient: scale-free measure, since
                                                        // the gradient scales as 1/scale of the data and the configuration
    else if (snorm < tolx * std::sqrt(dot(x, x))) reason = 4; // small step, relative to the size of the configuration
    R_CheckUserInterrupt();                             // allow the user to interrupt
  }
  f = evalcrit(M, x, g);                                // criterion and disparities at the final point
  NumericMatrix X(N, M.p);                              // final configuration
  std::copy(x.begin(), x.end(), X.begin());             // copy it
  return List::create(Named("conf") = X,                // configuration
                      Named("stress") = f,              // final value of the criterion
                      Named("dhat") = wrap(M.dh),       // disparities of the active pairs
                      Named("iter") = iter,             // iterations
                      Named("reason") = reason);        // stopping rule
}

// ---------------------------------------------------------------------------
// Criterion and gradient at a given configuration (used by the tests).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List mdscale_evalgrad(NumericVector dis, NumericVector w, IntegerVector ii, IntegerVector jj,
                      int N, NumericMatrix X, int crit, int ties,
                      NumericMatrix B2, NumericMatrix W1, NumericVector rw) {
  MdsData M;                                            // shared data, filled as in the engine
  M.N = N; M.p = X.ncol(); M.crit = crit; M.ties = ties; M.K = dis.size();
  M.ii = as<std::vector<int> >(ii); M.jj = as<std::vector<int> >(jj);
  M.dis = as<std::vector<double> >(dis); M.w = as<std::vector<double> >(w);
  M.dh.assign(M.K, 0.0); M.blk.assign(M.K, 0); M.nblk = 0;
  M.sumwD2 = 0.0; M.sumwD4 = 0.0; M.sumwD = 0.0;
  for (int k = 0; k < M.K; k++) {
    M.sumwD2 += M.w[k] * M.dis[k] * M.dis[k]; M.sumwD4 += M.w[k] * std::pow(M.dis[k], 4); M.sumwD += M.w[k] * M.dis[k];
  }
  M.ord.resize(M.K); std::iota(M.ord.begin(), M.ord.end(), 0);
  std::stable_sort(M.ord.begin(), M.ord.end(), [&M](int a, int b) { return M.dis[a] < M.dis[b]; });
  for (int k = 0; k < M.K; k++) if (k == 0 || M.dis[M.ord[k]] != M.dis[M.ord[k - 1]]) M.bstart.push_back(k);
  M.bstart.push_back(M.K);
  M.B2 = B2; M.W1 = W1; M.rw = as<std::vector<double> >(rw);
  std::vector<double> x(X.begin(), X.end()), g(x.size());
  double f = evalcrit(M, x, g);
  NumericMatrix G(N, M.p); std::copy(g.begin(), g.end(), G.begin());
  return List::create(Named("f") = f, Named("grad") = G, Named("dhat") = wrap(M.dh));
}
