// src/acf_estimate_rcpp.cpp
#include <Rcpp.h>
#include <cmath>
#include <cstdint>

using Rcpp::IntegerVector;
using Rcpp::LogicalVector;
using Rcpp::NumericVector;
using Rcpp::List;

namespace {

inline int idx3(const int x, const int y, const int z, const int nx, const int ny) {
  return x + nx * (y + ny * z);
}

inline int sgn_clip(const double v, const double clip) {
  if (!std::isfinite(v)) return 0;
  double x = v;
  if (clip > 0.0) {
    if (x > clip) x = clip;
    else if (x < -clip) x = -clip;
  }
  if (x > 0.0) return 1;
  if (x < 0.0) return -1;
  return 1; // treat 0 as positive
}

} // namespace

// Robust sign-correlation at a given integer lag.
// Returns mean(sign(Z(s))*sign(Z(s+h))) over s where both voxels are in mask.
//
// [[Rcpp::export]]
List acf_sign_lag_mean_cpp(const NumericVector& z_full,
                           const LogicalVector& mask_full,
                           const IntegerVector& dims,
                           const IntegerVector& lag_xyz,
                           const double clip = 3.0) {
  if (dims.size() != 3) Rcpp::stop("dims must have length 3");
  if (lag_xyz.size() != 3) Rcpp::stop("lag_xyz must have length 3");
  const int nx = dims[0], ny = dims[1], nz = dims[2];
  const int n = nx * ny * nz;
  if (z_full.size() != n) Rcpp::stop("z_full length must equal prod(dims)");
  if (mask_full.size() != n) Rcpp::stop("mask_full length must equal prod(dims)");

  const int dx = lag_xyz[0];
  const int dy = lag_xyz[1];
  const int dz = lag_xyz[2];

  double sum = 0.0;
  int64_t npairs = 0;

  // Iterate over valid origin voxels so (x+dx,y+dy,z+dz) stays in bounds.
  const int x0 = (dx >= 0) ? 0 : -dx;
  const int y0 = (dy >= 0) ? 0 : -dy;
  const int z0 = (dz >= 0) ? 0 : -dz;
  const int x1 = (dx >= 0) ? (nx - dx) : nx;
  const int y1 = (dy >= 0) ? (ny - dy) : ny;
  const int z1 = (dz >= 0) ? (nz - dz) : nz;

  for (int z = z0; z < z1; z++) {
    for (int y = y0; y < y1; y++) {
      for (int x = x0; x < x1; x++) {
        const int i0 = idx3(x, y, z, nx, ny);
        if (!mask_full[i0]) continue;
        const int i1 = idx3(x + dx, y + dy, z + dz, nx, ny);
        if (!mask_full[i1]) continue;

        const int s0 = sgn_clip(z_full[i0], clip);
        const int s1 = sgn_clip(z_full[i1], clip);
        if (s0 == 0 || s1 == 0) continue;

        sum += (double)(s0 * s1);
        npairs++;
      }
    }
  }

  const double mean = (npairs > 0) ? (sum / (double)npairs) : NA_REAL;
  return List::create(
    Rcpp::Named("mean") = mean,
    Rcpp::Named("n_pairs") = (double)npairs
  );
}

