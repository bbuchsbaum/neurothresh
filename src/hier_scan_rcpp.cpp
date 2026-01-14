// src/hier_scan_rcpp.cpp
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <climits>

using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::IntegerMatrix;
using Rcpp::List;

// ------------------------------------------------------------
// 1) Build octree split info in mask-space
// idx_parent: indices into mask-space vectors (R 1-based).
// x,y,z: integer coordinate vectors in mask-space (same length as Z_vec/pi_vec).
// bbox: integer length-6 = (x0,x1,y0,y1,z0,z1), inclusive coords.
// pi_vec: prior weights in mask-space.
// min_pi_mass: drop children with extremely tiny prior mass.
// ------------------------------------------------------------

// [[Rcpp::export]]
List octree_split_info_cpp(const IntegerVector& idx_parent,
                           const IntegerVector& bbox,
                           const IntegerVector& x,
                           const IntegerVector& y,
                           const IntegerVector& z,
                           const NumericVector& pi_vec,
                           const double min_pi_mass = 1e-10) {
  const int L = idx_parent.size();
  if (bbox.size() != 6) Rcpp::stop("bbox must have length 6");
  if (x.size() != pi_vec.size() || y.size() != pi_vec.size() || z.size() != pi_vec.size())
    Rcpp::stop("x,y,z must match length of pi_vec");

  const int x0 = bbox[0], x1 = bbox[1];
  const int y0 = bbox[2], y1 = bbox[3];
  const int z0 = bbox[4], z1 = bbox[5];

  const int xm = (x0 + x1) / 2;
  const int ym = (y0 + y1) / 2;
  const int zm = (z0 + z1) / 2;

  // 8 children in the dyadic split
  double den1_8[8]; double den2_8[8];
  for (int j = 0; j < 8; j++) { den1_8[j] = 0.0; den2_8[j] = 0.0; }

  IntegerVector child8(L);        // 0..7
  NumericVector w_parent(L);      // pi weights aligned to idx_parent

  for (int t = 0; t < L; t++) {
    const int ix = idx_parent[t] - 1; // to 0-based
    if (ix < 0 || ix >= pi_vec.size()) Rcpp::stop("idx_parent out of range");

    const double w = pi_vec[ix];
    w_parent[t] = w;

    int c = 0;
    if (x[ix] > xm) c += 1;
    if (y[ix] > ym) c += 2;
    if (z[ix] > zm) c += 4;

    child8[t] = c;

    den1_8[c] += w;
    den2_8[c] += w * w;
  }

  // Determine non-empty children (by prior mass)
  std::vector<int> keep8;
  keep8.reserve(8);
  for (int c = 0; c < 8; c++) {
    if (den1_8[c] >= min_pi_mass) keep8.push_back(c);
  }
  const int m = (int)keep8.size();
  if (m == 0) {
    return List::create(
      Rcpp::Named("m") = 0,
      Rcpp::Named("child_id") = IntegerVector(L, -1),
      Rcpp::Named("w_parent") = w_parent,
      Rcpp::Named("log_den1") = NumericVector(0),
      Rcpp::Named("inv_sqrt_den2") = NumericVector(0),
      Rcpp::Named("keep8") = IntegerVector(0)
    );
  }

  // Map 0..7 -> 0..m-1 or -1
  int map8toM[8];
  for (int c = 0; c < 8; c++) map8toM[c] = -1;
  for (int j = 0; j < m; j++) map8toM[ keep8[j] ] = j;

  IntegerVector child_id(L);              // -1 or 0..m-1
  for (int t = 0; t < L; t++) {
    child_id[t] = map8toM[ child8[t] ];
  }

  NumericVector log_den1(m);
  NumericVector inv_sqrt_den2(m);
  IntegerVector keep8_out(m);

  for (int j = 0; j < m; j++) {
    const int c = keep8[j];
    keep8_out[j] = c + 1; // return 1..8 for readability in R
    log_den1[j] = std::log(den1_8[c]);
    inv_sqrt_den2[j] = 1.0 / std::sqrt(den2_8[c]);
  }

  return List::create(
    Rcpp::Named("m") = m,
    Rcpp::Named("child_id") = child_id,
    Rcpp::Named("w_parent") = w_parent,
    Rcpp::Named("log_den1") = log_den1,
    Rcpp::Named("inv_sqrt_den2") = inv_sqrt_den2,
    Rcpp::Named("keep8") = keep8_out
  );
}

// ------------------------------------------------------------
// 2) One-pass child scoring for a single Z_vec (one permutation)
// child_id is aligned with idx_parent and w_parent,
// and takes values -1 (dropped) or 0..m-1.
// kappa_grid_pos should be > 0 values only.
// do_abs optionally applies fabs() to z.
// ------------------------------------------------------------

// [[Rcpp::export]]
NumericVector score_children_onepass_cpp(const NumericVector& Z_vec,
                                        const IntegerVector& idx_parent,
                                        const NumericVector& w_parent,
                                        const IntegerVector& child_id,
                                        const NumericVector& log_den1,
                                        const NumericVector& inv_sqrt_den2,
                                        const NumericVector& kappa_grid_pos,
                                        const bool do_abs = false) {
  const int L = idx_parent.size();
  const int m = log_den1.size();
  const int K = kappa_grid_pos.size();

  if (w_parent.size() != L || child_id.size() != L)
    Rcpp::stop("w_parent and child_id must match idx_parent length");
  if (inv_sqrt_den2.size() != m)
    Rcpp::stop("inv_sqrt_den2 must match log_den1 length");

  // kappa=0 accumulators
  std::vector<double> num0(m, 0.0);

  // kappa>0 streaming log-sum-exp accumulators
  const double NEG_INF = -std::numeric_limits<double>::infinity();
  std::vector<double> maxa(m * K, NEG_INF);
  std::vector<double> sumexp(m * K, 0.0);

  for (int t = 0; t < L; t++) {
    const int j = child_id[t];
    if (j < 0 || j >= m) continue;

    const int ix = idx_parent[t] - 1;
    double zval = Z_vec[ix];
    if (Rcpp::NumericVector::is_na(zval)) continue;
    if (do_abs) zval = std::fabs(zval);

    const double w = w_parent[t];

    // kappa=0
    num0[j] += w * zval;

    // kappa>0
    for (int kk = 0; kk < K; kk++) {
      const double kappa = kappa_grid_pos[kk];
      const double a = kappa * zval;
      const int pos = j + m * kk;

      const double mx = maxa[pos];
      if (a > mx) {
        // sumexp <- sumexp * exp(mx-a) + w  (since exp(a-a)=1 for the new max element)
        if (mx == NEG_INF) {
          sumexp[pos] = w;
        } else {
          sumexp[pos] = sumexp[pos] * std::exp(mx - a) + w;
        }
        maxa[pos] = a;
      } else {
        sumexp[pos] += w * std::exp(a - mx);
      }
    }
  }

  NumericVector out(m);

  for (int j = 0; j < m; j++) {
    // U0 = sum(pi Z) / sqrt(sum(pi^2))
    const double U0 = num0[j] * inv_sqrt_den2[j];

    double best_soft = NEG_INF;
    for (int kk = 0; kk < K; kk++) {
      const double kappa = kappa_grid_pos[kk];
      const int pos = j + m * kk;

      if (sumexp[pos] <= 0.0) continue; // should not happen if child non-empty

      const double logsum = std::log(sumexp[pos]) + maxa[pos];
      const double S = (logsum - log_den1[j]) / kappa;

      if (S > best_soft) best_soft = S;
    }

    out[j] = (best_soft > U0) ? best_soft : U0;
  }

  return out;
}

// ------------------------------------------------------------
// 3) Split idx_parent into idx_child lists (for recursion)
// Returns a list of IntegerVectors of mask-space indices.
// ------------------------------------------------------------

// [[Rcpp::export]]
List split_indices_cpp(const IntegerVector& idx_parent,
                       const IntegerVector& child_id,
                       const int m) {
  const int L = idx_parent.size();
  if (child_id.size() != L) Rcpp::stop("child_id must match idx_parent length");
  if (m <= 0) return List::create();

  std::vector<int> counts(m, 0);
  for (int t = 0; t < L; t++) {
    const int j = child_id[t];
    if (j >= 0 && j < m) counts[j]++;
  }

  // allocate
  std::vector< IntegerVector > children;
  children.reserve(m);
  for (int j = 0; j < m; j++) children.push_back(IntegerVector(counts[j]));

  // fill
  std::vector<int> pos(m, 0);
  for (int t = 0; t < L; t++) {
    const int j = child_id[t];
    if (j >= 0 && j < m) {
      children[j][ pos[j]++ ] = idx_parent[t];
    }
  }

  // wrap into List
  List out(m);
  for (int j = 0; j < m; j++) out[j] = children[j];
  return out;
}
