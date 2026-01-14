// src/bbox_helpers.cpp
#include <Rcpp.h>
#include <climits>
#include <vector>

using Rcpp::IntegerVector;
using Rcpp::IntegerMatrix;
using Rcpp::NumericVector;
using Rcpp::List;

// ------------------------------------------------------------
// bbox_from_indices_cpp(): Bbox for Any Mask-Space Index List
// ------------------------------------------------------------

// [[Rcpp::export]]
IntegerVector bbox_from_indices_cpp(const IntegerVector& idx,
                                   const IntegerVector& x,
                                   const IntegerVector& y,
                                   const IntegerVector& z) {
  const int n = idx.size();
  if (x.size() != y.size() || x.size() != z.size())
    Rcpp::stop("x,y,z must have same length");

  if (n == 0) {
    return IntegerVector::create(NA_INTEGER, NA_INTEGER,
                                 NA_INTEGER, NA_INTEGER,
                                 NA_INTEGER, NA_INTEGER);
  }

  int minx = INT_MAX, maxx = INT_MIN;
  int miny = INT_MAX, maxy = INT_MIN;
  int minz = INT_MAX, maxz = INT_MIN;

  for (int i = 0; i < n; i++) {
    int ix = idx[i] - 1;  // mask-space index, R->C++
    if (ix < 0 || ix >= x.size()) Rcpp::stop("idx out of range");

    const int xi = x[ix];
    const int yi = y[ix];
    const int zi = z[ix];

    if (xi < minx) minx = xi; if (xi > maxx) maxx = xi;
    if (yi < miny) miny = yi; if (yi > maxy) maxy = yi;
    if (zi < minz) minz = zi; if (zi > maxz) maxz = zi;
  }

  return IntegerVector::create(minx, maxx, miny, maxy, minz, maxz);
}

// ------------------------------------------------------------
// split_indices_bbox_cpp(): Split Parent into Children and Compute Each Child Bbox
// ------------------------------------------------------------

// [[Rcpp::export]]
List split_indices_bbox_cpp(const IntegerVector& idx_parent,
                            const IntegerVector& child_id, // -1 or 0..m-1
                            const int m,
                            const IntegerVector& x,
                            const IntegerVector& y,
                            const IntegerVector& z) {
  const int L = idx_parent.size();
  if (child_id.size() != L) Rcpp::stop("child_id must match idx_parent length");
  if (m <= 0) return List::create(Rcpp::Named("idx") = List(0),
                                  Rcpp::Named("bbox") = IntegerMatrix(0, 6));

  if (x.size() != y.size() || x.size() != z.size())
    Rcpp::stop("x,y,z must have same length");

  // First pass: counts + bbox mins/maxs
  std::vector<int> counts(m, 0);

  std::vector<int> minx(m, INT_MAX), maxx(m, INT_MIN);
  std::vector<int> miny(m, INT_MAX), maxy(m, INT_MIN);
  std::vector<int> minz(m, INT_MAX), maxz(m, INT_MIN);

  for (int t = 0; t < L; t++) {
    const int j = child_id[t];
    if (j < 0 || j >= m) continue;

    const int ix = idx_parent[t] - 1;
    if (ix < 0 || ix >= x.size()) Rcpp::stop("idx_parent out of range");

    counts[j]++;

    const int xi = x[ix], yi = y[ix], zi = z[ix];
    if (xi < minx[j]) minx[j] = xi; if (xi > maxx[j]) maxx[j] = xi;
    if (yi < miny[j]) miny[j] = yi; if (yi > maxy[j]) maxy[j] = yi;
    if (zi < minz[j]) minz[j] = zi; if (zi > maxz[j]) maxz[j] = zi;
  }

  // Allocate child index vectors
  std::vector< IntegerVector > children;
  children.reserve(m);
  for (int j = 0; j < m; j++) children.push_back(IntegerVector(counts[j]));

  // Second pass: fill child indices
  std::vector<int> pos(m, 0);
  for (int t = 0; t < L; t++) {
    const int j = child_id[t];
    if (j < 0 || j >= m) continue;
    children[j][ pos[j]++ ] = idx_parent[t];
  }

  // Pack bboxes into an m x 6 matrix
  IntegerMatrix bbox(m, 6);
  for (int j = 0; j < m; j++) {
    if (counts[j] == 0) {
      bbox(j,0)=NA_INTEGER; bbox(j,1)=NA_INTEGER;
      bbox(j,2)=NA_INTEGER; bbox(j,3)=NA_INTEGER;
      bbox(j,4)=NA_INTEGER; bbox(j,5)=NA_INTEGER;
    } else {
      bbox(j,0)=minx[j]; bbox(j,1)=maxx[j];
      bbox(j,2)=miny[j]; bbox(j,3)=maxy[j];
      bbox(j,4)=minz[j]; bbox(j,5)=maxz[j];
    }
  }

  // Wrap children into an R list
  List idx_list(m);
  for (int j = 0; j < m; j++) idx_list[j] = children[j];

  return List::create(Rcpp::Named("idx") = idx_list,
                      Rcpp::Named("bbox") = bbox);
}

// ------------------------------------------------------------
// split_indices_bbox_mass_cpp(): Split + Bbox + Prior Masses
// ------------------------------------------------------------

// [[Rcpp::export]]
List split_indices_bbox_mass_cpp(const IntegerVector& idx_parent,
                                 const IntegerVector& child_id, // -1 or 0..m-1
                                 const int m,
                                 const IntegerVector& x,
                                 const IntegerVector& y,
                                 const IntegerVector& z,
                                 const NumericVector& pi_vec) {

  const int L = idx_parent.size();
  if (child_id.size() != L) Rcpp::stop("child_id must match idx_parent length");
  if (m <= 0) {
    return List::create(Rcpp::Named("idx") = List(0),
                        Rcpp::Named("bbox") = IntegerMatrix(0, 6),
                        Rcpp::Named("den1") = NumericVector(0),
                        Rcpp::Named("den2") = NumericVector(0));
  }
  if (x.size() != y.size() || x.size() != z.size() || x.size() != pi_vec.size())
    Rcpp::stop("x,y,z,pi_vec must have the same length (mask-space)");

  // First pass: counts + bbox mins/maxs + den1/den2
  std::vector<int> counts(m, 0);

  std::vector<int> minx(m, INT_MAX), maxx(m, INT_MIN);
  std::vector<int> miny(m, INT_MAX), maxy(m, INT_MIN);
  std::vector<int> minz(m, INT_MAX), maxz(m, INT_MIN);

  std::vector<double> den1(m, 0.0);
  std::vector<double> den2(m, 0.0);

  for (int t = 0; t < L; t++) {
    const int j = child_id[t];
    if (j < 0 || j >= m) continue;

    const int ix = idx_parent[t] - 1; // R->C++ index shift
    if (ix < 0 || ix >= x.size()) Rcpp::stop("idx_parent out of range");

    counts[j]++;

    // bbox
    const int xi = x[ix], yi = y[ix], zi = z[ix];
    if (xi < minx[j]) minx[j] = xi; if (xi > maxx[j]) maxx[j] = xi;
    if (yi < miny[j]) miny[j] = yi; if (yi > maxy[j]) maxy[j] = yi;
    if (zi < minz[j]) minz[j] = zi; if (zi > maxz[j]) maxz[j] = zi;

    // masses
    const double w = pi_vec[ix];
    den1[j] += w;
    den2[j] += w * w;
  }

  // Allocate child index vectors
  std::vector< IntegerVector > children;
  children.reserve(m);
  for (int j = 0; j < m; j++) children.push_back(IntegerVector(counts[j]));

  // Second pass: fill child indices
  std::vector<int> pos(m, 0);
  for (int t = 0; t < L; t++) {
    const int j = child_id[t];
    if (j < 0 || j >= m) continue;
    children[j][ pos[j]++ ] = idx_parent[t];
  }

  // Pack bboxes into an m x 6 matrix
  IntegerMatrix bbox(m, 6);
  for (int j = 0; j < m; j++) {
    if (counts[j] == 0) {
      bbox(j,0)=NA_INTEGER; bbox(j,1)=NA_INTEGER;
      bbox(j,2)=NA_INTEGER; bbox(j,3)=NA_INTEGER;
      bbox(j,4)=NA_INTEGER; bbox(j,5)=NA_INTEGER;
    } else {
      bbox(j,0)=minx[j]; bbox(j,1)=maxx[j];
      bbox(j,2)=miny[j]; bbox(j,3)=maxy[j];
      bbox(j,4)=minz[j]; bbox(j,5)=maxz[j];
    }
  }

  // Wrap idx list
  List idx_list(m);
  for (int j = 0; j < m; j++) idx_list[j] = children[j];

  // Wrap den1/den2
  NumericVector den1_out(m), den2_out(m);
  for (int j = 0; j < m; j++) { den1_out[j] = den1[j]; den2_out[j] = den2[j]; }

  return List::create(Rcpp::Named("idx") = idx_list,
                      Rcpp::Named("bbox") = bbox,
                      Rcpp::Named("den1") = den1_out,
                      Rcpp::Named("den2") = den2_out);
}
