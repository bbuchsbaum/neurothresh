// src/smooth3d_rcpp.cpp
#include <Rcpp.h>
#include <vector>
#include <cmath>

using Rcpp::IntegerVector;
using Rcpp::NumericVector;

namespace {

std::vector<double> gaussian_kernel_1d(const double sigma, const int radius) {
  const int K = 2 * radius + 1;
  std::vector<double> k(K);
  if (sigma <= 0.0) {
    k[radius] = 1.0;
    return k;
  }
  const double inv2s2 = 1.0 / (2.0 * sigma * sigma);
  double s = 0.0;
  for (int i = -radius; i <= radius; i++) {
    const double v = std::exp(- (double)(i * i) * inv2s2);
    k[i + radius] = v;
    s += v;
  }
  if (s > 0.0) {
    for (int i = 0; i < K; i++) k[i] /= s;
  }
  return k;
}

inline int idx3(const int x, const int y, const int z, const int nx, const int ny) {
  return x + nx * (y + ny * z);
}

} // namespace

// [[Rcpp::export]]
NumericVector gaussian_smooth3d_cpp(const NumericVector& vol,
                                   const IntegerVector& dims,
                                   const NumericVector& sigma_xyz,
                                   const int radius_mult = 3) {
  if (dims.size() != 3) Rcpp::stop("dims must have length 3");
  const int nx = dims[0], ny = dims[1], nz = dims[2];
  const int n = nx * ny * nz;
  if (vol.size() != n) Rcpp::stop("vol length must equal prod(dims)");
  if (sigma_xyz.size() != 3) Rcpp::stop("sigma_xyz must have length 3");
  if (radius_mult < 1) Rcpp::stop("radius_mult must be >= 1");

  const double sx = sigma_xyz[0];
  const double sy = sigma_xyz[1];
  const double sz = sigma_xyz[2];

  const int rx = (sx > 0.0) ? (int)std::ceil(radius_mult * sx) : 0;
  const int ry = (sy > 0.0) ? (int)std::ceil(radius_mult * sy) : 0;
  const int rz = (sz > 0.0) ? (int)std::ceil(radius_mult * sz) : 0;

  const std::vector<double> kx = gaussian_kernel_1d(sx, rx);
  const std::vector<double> ky = gaussian_kernel_1d(sy, ry);
  const std::vector<double> kz = gaussian_kernel_1d(sz, rz);

  std::vector<double> tmp1(n, 0.0);
  std::vector<double> tmp2(n, 0.0);

  // X pass
  for (int z = 0; z < nz; z++) {
    for (int y = 0; y < ny; y++) {
      for (int x = 0; x < nx; x++) {
        double acc = 0.0;
        for (int dx = -rx; dx <= rx; dx++) {
          const int xx = x + dx;
          if (xx < 0 || xx >= nx) continue;
          acc += kx[dx + rx] * vol[idx3(xx, y, z, nx, ny)];
        }
        tmp1[idx3(x, y, z, nx, ny)] = acc;
      }
    }
  }

  // Y pass
  for (int z = 0; z < nz; z++) {
    for (int x = 0; x < nx; x++) {
      for (int y = 0; y < ny; y++) {
        double acc = 0.0;
        for (int dy = -ry; dy <= ry; dy++) {
          const int yy = y + dy;
          if (yy < 0 || yy >= ny) continue;
          acc += ky[dy + ry] * tmp1[idx3(x, yy, z, nx, ny)];
        }
        tmp2[idx3(x, y, z, nx, ny)] = acc;
      }
    }
  }

  // Z pass
  NumericVector out(n);
  for (int y = 0; y < ny; y++) {
    for (int x = 0; x < nx; x++) {
      for (int z = 0; z < nz; z++) {
        double acc = 0.0;
        for (int dz = -rz; dz <= rz; dz++) {
          const int zz = z + dz;
          if (zz < 0 || zz >= nz) continue;
          acc += kz[dz + rz] * tmp2[idx3(x, y, zz, nx, ny)];
        }
        out[idx3(x, y, z, nx, ny)] = acc;
      }
    }
  }

  return out;
}

