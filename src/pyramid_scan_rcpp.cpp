// src/pyramid_scan_rcpp.cpp
#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <cstdint>
#include <limits>

using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::List;

namespace {

constexpr uint64_t SHIFT_Y = 21;
constexpr uint64_t SHIFT_Z = 42;
constexpr uint64_t MASK_21 = (1ULL << 21) - 1ULL;

inline uint64_t pack3(const uint32_t ix, const uint32_t iy, const uint32_t iz) {
  return (uint64_t)ix | ((uint64_t)iy << SHIFT_Y) | ((uint64_t)iz << SHIFT_Z);
}

inline void unpack3(const uint64_t key, uint32_t& ix, uint32_t& iy, uint32_t& iz) {
  ix = (uint32_t)(key & MASK_21);
  iy = (uint32_t)((key >> SHIFT_Y) & MASK_21);
  iz = (uint32_t)((key >> SHIFT_Z) & MASK_21);
}

inline int ilog2_pow2(const int n) {
  // assumes n is a power of two and n >= 1
  int k = 0;
  int v = n;
  while (v > 1) {
    v >>= 1;
    k++;
  }
  return k;
}

struct DenPyramid {
  int side = 0;
  int n_levels = 0; // includes level 0 leaves
  std::vector< std::unordered_map<uint64_t, double> > inv_sqrt_den2;
};

} // namespace

// [[Rcpp::export]]
SEXP build_den_pyramid_xptr_cpp(const IntegerVector& x0,
                                const IntegerVector& y0,
                                const IntegerVector& z0,
                                const NumericVector& pi_vec,
                                const int side) {
  const int n = pi_vec.size();
  if (x0.size() != n || y0.size() != n || z0.size() != n) {
    Rcpp::stop("x0,y0,z0 must match length of pi_vec");
  }
  if (side <= 0) Rcpp::stop("side must be positive");

  const int L = ilog2_pow2(side);

  DenPyramid* den = new DenPyramid();
  den->side = side;
  den->n_levels = L + 1;
  den->inv_sqrt_den2.resize(den->n_levels);

  std::unordered_map<uint64_t, double> cur_sum;
  cur_sum.reserve((size_t)(n * 1.3) + 16);

  for (int i = 0; i < n; i++) {
    const double w = pi_vec[i];
    if (!std::isfinite(w) || w <= 0.0) continue;

    const int ix = x0[i];
    const int iy = y0[i];
    const int iz = z0[i];

    if (ix < 0 || iy < 0 || iz < 0 || ix >= side || iy >= side || iz >= side) {
      Rcpp::stop("x0,y0,z0 out of range for side");
    }

    const uint64_t key = pack3((uint32_t)ix, (uint32_t)iy, (uint32_t)iz);
    cur_sum[key] += w * w;
  }

  for (int l = 0; l <= L; l++) {
    auto& inv_map = den->inv_sqrt_den2[l];
    inv_map.reserve(cur_sum.size() + 16);

    for (const auto& kv : cur_sum) {
      const double d2 = kv.second;
      if (d2 <= 0.0) continue;
      inv_map.emplace(kv.first, 1.0 / std::sqrt(d2));
    }

    if (l == L) break;

    std::unordered_map<uint64_t, double> next_sum;
    next_sum.reserve((size_t)(cur_sum.size() / 6) + 16);

    for (const auto& kv : cur_sum) {
      uint32_t ix, iy, iz;
      unpack3(kv.first, ix, iy, iz);
      const uint64_t parent = pack3(ix >> 1, iy >> 1, iz >> 1);
      next_sum[parent] += kv.second;
    }

    cur_sum.swap(next_sum);
  }

  Rcpp::XPtr<DenPyramid> xp(den, true);
  return xp;
}

namespace {

inline double pyramid_max_u0_impl(const NumericVector& Z_vec,
                                  const NumericVector& pi_vec,
                                  const IntegerVector& x0,
                                  const IntegerVector& y0,
                                  const IntegerVector& z0,
                                  DenPyramid* den,
                                  const bool do_abs,
                                  const bool do_signflip) {
  const int n = Z_vec.size();
  if (pi_vec.size() != n || x0.size() != n || y0.size() != n || z0.size() != n) {
    Rcpp::stop("Z_vec, pi_vec, x0,y0,z0 must have the same length");
  }

  std::unordered_map<uint64_t, double> cur_sum;
  cur_sum.reserve((size_t)(n * 1.3) + 16);

  for (int i = 0; i < n; i++) {
    const double w = pi_vec[i];
    double z = Z_vec[i];
    if (!std::isfinite(w) || w <= 0.0) continue;
    if (!std::isfinite(z)) continue;

    if (do_signflip) {
      z *= (R::unif_rand() < 0.5) ? -1.0 : 1.0;
    }
    if (do_abs) z = std::fabs(z);

    const uint64_t key = pack3((uint32_t)x0[i], (uint32_t)y0[i], (uint32_t)z0[i]);
    cur_sum[key] += w * z;
  }

  double max_score = -std::numeric_limits<double>::infinity();

  for (int l = 0; l < den->n_levels; l++) {
    const auto& inv_map = den->inv_sqrt_den2[l];

    for (const auto& kv : cur_sum) {
      auto it = inv_map.find(kv.first);
      if (it == inv_map.end()) continue;
      const double score = kv.second * it->second;
      if (score > max_score) max_score = score;
    }

    if (l + 1 >= den->n_levels) break;

    std::unordered_map<uint64_t, double> next_sum;
    next_sum.reserve((size_t)(cur_sum.size() / 6) + 16);

    for (const auto& kv : cur_sum) {
      uint32_t ix, iy, iz;
      unpack3(kv.first, ix, iy, iz);
      const uint64_t parent = pack3(ix >> 1, iy >> 1, iz >> 1);
      next_sum[parent] += kv.second;
    }

    cur_sum.swap(next_sum);
  }

  return max_score;
}

} // namespace

// [[Rcpp::export]]
double pyramid_max_u0_cpp(const NumericVector& Z_vec,
                          const NumericVector& pi_vec,
                          const IntegerVector& x0,
                          const IntegerVector& y0,
                          const IntegerVector& z0,
                          SEXP den_xptr,
                          const bool do_abs = false) {
  Rcpp::XPtr<DenPyramid> den(den_xptr);
  Rcpp::RNGScope scope;
  return pyramid_max_u0_impl(Z_vec, pi_vec, x0, y0, z0, den.get(), do_abs, false);
}

// [[Rcpp::export]]
double pyramid_max_u0_signflip_cpp(const NumericVector& Z_vec,
                                   const NumericVector& pi_vec,
                                   const IntegerVector& x0,
                                   const IntegerVector& y0,
                                   const IntegerVector& z0,
                                   SEXP den_xptr,
                                   const bool do_abs = false) {
  Rcpp::XPtr<DenPyramid> den(den_xptr);
  Rcpp::RNGScope scope;
  return pyramid_max_u0_impl(Z_vec, pi_vec, x0, y0, z0, den.get(), do_abs, true);
}

// [[Rcpp::export]]
List pyramid_scores_u0_cpp(const NumericVector& Z_vec,
                           const NumericVector& pi_vec,
                           const IntegerVector& x0,
                           const IntegerVector& y0,
                           const IntegerVector& z0,
                           SEXP den_xptr,
                           const bool do_abs = false) {
  Rcpp::XPtr<DenPyramid> den(den_xptr);
  Rcpp::RNGScope scope;

  const int n = Z_vec.size();
  if (pi_vec.size() != n || x0.size() != n || y0.size() != n || z0.size() != n) {
    Rcpp::stop("Z_vec, pi_vec, x0,y0,z0 must have the same length");
  }

  std::unordered_map<uint64_t, double> cur_sum;
  cur_sum.reserve((size_t)(n * 1.3) + 16);

  for (int i = 0; i < n; i++) {
    const double w = pi_vec[i];
    double z = Z_vec[i];
    if (!std::isfinite(w) || w <= 0.0) continue;
    if (!std::isfinite(z)) continue;
    if (do_abs) z = std::fabs(z);

    const uint64_t key = pack3((uint32_t)x0[i], (uint32_t)y0[i], (uint32_t)z0[i]);
    cur_sum[key] += w * z;
  }

  std::vector<int> out_level;
  std::vector<int> out_i;
  std::vector<int> out_j;
  std::vector<int> out_k;
  std::vector<double> out_score;

  out_level.reserve((size_t)(n * 1.2) + 64);
  out_i.reserve(out_level.capacity());
  out_j.reserve(out_level.capacity());
  out_k.reserve(out_level.capacity());
  out_score.reserve(out_level.capacity());

  double max_score = -std::numeric_limits<double>::infinity();

  for (int l = 0; l < den->n_levels; l++) {
    const auto& inv_map = den->inv_sqrt_den2[l];

    for (const auto& kv : cur_sum) {
      auto it = inv_map.find(kv.first);
      if (it == inv_map.end()) continue;

      const double score = kv.second * it->second;
      if (score > max_score) max_score = score;

      uint32_t ix, iy, iz;
      unpack3(kv.first, ix, iy, iz);

      out_level.push_back(l);
      out_i.push_back((int)ix + 1);
      out_j.push_back((int)iy + 1);
      out_k.push_back((int)iz + 1);
      out_score.push_back(score);
    }

    if (l + 1 >= den->n_levels) break;

    std::unordered_map<uint64_t, double> next_sum;
    next_sum.reserve((size_t)(cur_sum.size() / 6) + 16);

    for (const auto& kv : cur_sum) {
      uint32_t ix, iy, iz;
      unpack3(kv.first, ix, iy, iz);
      const uint64_t parent = pack3(ix >> 1, iy >> 1, iz >> 1);
      next_sum[parent] += kv.second;
    }

    cur_sum.swap(next_sum);
  }

  List nodes = List::create(
    Rcpp::Named("level") = out_level,
    Rcpp::Named("i") = out_i,
    Rcpp::Named("j") = out_j,
    Rcpp::Named("k") = out_k,
    Rcpp::Named("score") = out_score
  );
  nodes.attr("class") = "data.frame";
  nodes.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -(int)out_score.size());

  return List::create(
    Rcpp::Named("max_score") = max_score,
    Rcpp::Named("nodes") = nodes,
    Rcpp::Named("side") = den->side,
    Rcpp::Named("n_levels") = den->n_levels
  );
}
