// src/exp_smooth3d_rcpp.cpp
#include <Rcpp.h>
#include <vector>
#include <cmath>

using Rcpp::IntegerVector;
using Rcpp::NumericVector;

namespace {

inline int idx3(const int x, const int y, const int z, const int nx, const int ny) {
  return x + nx * (y + ny * z);
}

std::vector<double> exp_kernel_1d(const double lambda, const int radius) {
  const int K = 2 * radius + 1;
  std::vector<double> k(K);
  if (lambda <= 0.0) {
    k[radius] = 1.0;
    return k;
  }
  double s = 0.0;
  for (int i = -radius; i <= radius; i++) {
    const double v = std::exp(-std::fabs((double)i) / lambda);
    k[i + radius] = v;
    s += v;
  }
  if (s > 0.0) {
    for (int i = 0; i < K; i++) k[i] /= s;
  }
  return k;
}

} // namespace

// [[Rcpp::export]]
NumericVector exp_smooth3d_cpp(const NumericVector& vol,
                              const IntegerVector& dims,
                              const NumericVector& lambda_xyz,
                              const int radius_mult = 6) {
  if (dims.size() != 3) Rcpp::stop("dims must have length 3");
  const int nx = dims[0], ny = dims[1], nz = dims[2];
  const int n = nx * ny * nz;
  if (vol.size() != n) Rcpp::stop("vol length must equal prod(dims)");
  if (lambda_xyz.size() != 3) Rcpp::stop("lambda_xyz must have length 3");
  if (radius_mult < 1) Rcpp::stop("radius_mult must be >= 1");

  const double lx = lambda_xyz[0];
  const double ly = lambda_xyz[1];
  const double lz = lambda_xyz[2];

  const int rx = (lx > 0.0) ? (int)std::ceil(radius_mult * lx) : 0;
  const int ry = (ly > 0.0) ? (int)std::ceil(radius_mult * ly) : 0;
  const int rz = (lz > 0.0) ? (int)std::ceil(radius_mult * lz) : 0;

  const std::vector<double> kx = exp_kernel_1d(lx, rx);
  const std::vector<double> ky = exp_kernel_1d(ly, ry);
  const std::vector<double> kz = exp_kernel_1d(lz, rz);

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

