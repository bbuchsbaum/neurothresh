// src/spatial_rcpp.cpp
//
// R bindings for the header-only spatial core (inst/include/neurothresh/
// spatial_core.h): TFCE and cluster labelling on masked 3D grids.
#include <Rcpp.h>
// Relative path so no Makevars include flag is needed (src/Makevars is
// git-ignored in this repository).
#include "../inst/include/neurothresh/spatial_core.h"

#include <vector>

using Rcpp::IntegerVector;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

namespace {

void read_dims(const IntegerVector& dims, int out[3]) {
  if (dims.size() != 3) Rcpp::stop("dims must have length 3");
  for (int d = 0; d < 3; ++d) {
    if (dims[d] == NA_INTEGER || dims[d] < 1) Rcpp::stop("dims must be positive");
    out[d] = dims[d];
  }
}

ntcore::Graph make_graph(const IntegerVector& dims, const IntegerVector& mask_idx,
                         int connectivity) {
  int d[3];
  read_dims(dims, d);
  std::vector<int> vox(mask_idx.size());
  for (R_xlen_t i = 0; i < mask_idx.size(); ++i) {
    if (mask_idx[i] == NA_INTEGER) Rcpp::stop("mask_idx contains NA");
    vox[i] = mask_idx[i] - 1;  // R is 1-based
  }
  try {
    return ntcore::build_grid_graph(d, vox, connectivity);
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  }
  return ntcore::Graph();  // not reached
}

void check_tfce_params(double H, double E, double dh) {
  if (!(H > 0.0) || !(E > 0.0) || !(dh > 0.0)) {
    Rcpp::stop("H, E and dh must be positive");
  }
}

}  // namespace

// TFCE of in-mask values.
//
// z: values at the mask voxels (length(mask_idx)).
// tail: 1 = positive, -1 = negative, 0 = two-sided (per-voxel maximum of the
//   two tails).
// Returns the TFCE values at the mask voxels.
// [[Rcpp::export(rng = false)]]
NumericVector tfce_masked_cpp(const NumericVector& z, const IntegerVector& dims,
                              const IntegerVector& mask_idx, double H, double E,
                              double dh, int connectivity, int tail) {
  check_tfce_params(H, E, dh);
  if (z.size() != mask_idx.size()) Rcpp::stop("length(z) must equal length(mask_idx)");
  ntcore::Graph g = make_graph(dims, mask_idx, connectivity);
  const int n = g.n;
  ntcore::TfceWork w;
  w.resize(n);
  NumericVector out(n);
  std::vector<double> x(z.begin(), z.end());
  if (tail == 1 || tail == 0) {
    ntcore::tfce(g, x.data(), H, E, dh, out.begin(), w);
  }
  if (tail == -1 || tail == 0) {
    std::vector<double> neg(n), tneg(n);
    for (int v = 0; v < n; ++v) neg[v] = -x[v];
    ntcore::tfce(g, neg.data(), H, E, dh, tneg.data(), w);
    if (tail == -1) {
      for (int v = 0; v < n; ++v) out[v] = tneg[v];
    } else {
      for (int v = 0; v < n; ++v) out[v] = std::max(out[v], tneg[v]);
    }
  }
  return out;
}

// Maximum TFCE of each column of Z (mask voxels x draws). Builds the
// neighbour graph once for the whole batch.
// [[Rcpp::export(rng = false)]]
NumericVector tfce_max_batch_cpp(const NumericMatrix& Z, const IntegerVector& dims,
                                 const IntegerVector& mask_idx, double H, double E,
                                 double dh, int connectivity, int tail) {
  check_tfce_params(H, E, dh);
  if (Z.nrow() != mask_idx.size()) Rcpp::stop("nrow(Z) must equal length(mask_idx)");
  ntcore::Graph g = make_graph(dims, mask_idx, connectivity);
  const int n = g.n;
  const int B = Z.ncol();
  ntcore::TfceWork w;
  w.resize(n);
  std::vector<double> x(n), out(n);
  NumericVector res(B);
  for (int b = 0; b < B; ++b) {
    const double* col = &Z(0, b);
    double mx = 0.0;
    if (tail == 1 || tail == 0) {
      mx = ntcore::tfce(g, col, H, E, dh, out.data(), w);
    }
    if (tail == -1 || tail == 0) {
      for (int v = 0; v < n; ++v) x[v] = -col[v];
      mx = std::max(mx, ntcore::tfce(g, x.data(), H, E, dh, out.data(), w));
    }
    res[b] = mx;
    if ((b & 63) == 63) Rcpp::checkUserInterrupt();
  }
  return res;
}

// Connected components of {x >= u} among the mask voxels.
//
// Returns labels (per mask voxel, 0 = none, clusters numbered by their
// smallest mask position), sizes, excess mass sum(x - u) and the mask
// position (1-based) of each cluster's peak.
// [[Rcpp::export(rng = false)]]
Rcpp::List cluster_label_cpp(const NumericVector& x, const IntegerVector& dims,
                             const IntegerVector& mask_idx, double u,
                             int connectivity) {
  if (x.size() != mask_idx.size()) Rcpp::stop("length(x) must equal length(mask_idx)");
  if (!std::isfinite(u)) Rcpp::stop("threshold must be finite");
  ntcore::Graph g = make_graph(dims, mask_idx, connectivity);
  ntcore::ClusterWork w;
  w.resize(g.n);
  IntegerVector labels(g.n);
  std::vector<int> sizes, peak;
  std::vector<double> mass;
  std::vector<double> xv(x.begin(), x.end());
  ntcore::cluster_label(g, xv.data(), u, w, labels.begin(), &sizes, &mass, &peak);
  IntegerVector peak1(peak.size());
  for (std::size_t c = 0; c < peak.size(); ++c) peak1[c] = peak[c] + 1;
  return Rcpp::List::create(
    Rcpp::Named("labels") = labels,
    Rcpp::Named("sizes") = IntegerVector(sizes.begin(), sizes.end()),
    Rcpp::Named("mass") = NumericVector(mass.begin(), mass.end()),
    Rcpp::Named("peak") = peak1
  );
}
