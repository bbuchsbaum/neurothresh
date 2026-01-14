# Rcpp Implementation Specification

## Source Files

All C++ code is in `src/`. Production code is provided in discussion notes - **USE VERBATIM**.

| File | Functions | Source |
|------|-----------|--------|
| `bbox_helpers.cpp` | bbox_from_indices_cpp, split_indices_bbox_cpp, split_indices_bbox_mass_cpp | part-05 |
| `split_indices.cpp` | split_indices_cpp | part-04 |
| `score_children.cpp` | octree_split_info_cpp, score_children_onepass_cpp | part-04 |

---

## C++ Function Signatures

### From part-04-onepass-rcpp.md

```cpp
// [[Rcpp::export]]
List octree_split_info_cpp(
  IntegerVector indices,    // 1-based voxel indices
  IntegerVector dims        // volume dimensions [nx, ny, nz]
);
// Returns: List(octants=IntegerVector, bbox=IntegerVector)

// [[Rcpp::export]]
NumericVector score_children_onepass_cpp(
  IntegerVector indices,    // parent voxel indices
  IntegerVector octants,    // octant assignment (0-7) per voxel
  NumericVector z_vals,     // Z-scores at indices
  NumericVector prior_vals, // prior weights at indices
  double kappa              // temperature parameter
);
// Returns: NumericVector of length 8 (scores for each child)

// [[Rcpp::export]]
List split_indices_cpp(
  IntegerVector indices,
  IntegerVector dims
);
// Returns: List of 8 IntegerVectors (child indices)
```

### From part-05-bbox-helpers-and-baselines.md

```cpp
// [[Rcpp::export]]
IntegerVector bbox_from_indices_cpp(
  IntegerVector indices,    // 1-based voxel indices
  IntegerVector dims        // volume dimensions
);
// Returns: IntegerVector(i0, i1, j0, j1, k0, k1) - 0-based bbox

// [[Rcpp::export]]
List split_indices_bbox_cpp(
  IntegerVector indices,
  IntegerVector dims
);
// Returns: List of 8 IntegerVectors

// [[Rcpp::export]]
List split_indices_bbox_mass_cpp(
  IntegerVector indices,
  IntegerVector dims,
  NumericVector prior_vec   // prior weight per voxel in volume
);
// Returns: List(children=List of 8, masses=NumericVector of 8)
```

---

## Implementation Notes

### Index Convention
- R uses 1-based indices
- C++ internally converts to 0-based
- All returned indices are 1-based for R compatibility

### Coordinate Conversion
```cpp
// 1-based index to 0-based (i,j,k)
int idx0 = idx - 1;  // convert to 0-based
int i = idx0 % nx;
int j = (idx0 / nx) % ny;
int k = idx0 / (nx * ny);

// 0-based (i,j,k) to 1-based index
int idx = i + j * nx + k * nx * ny + 1;
```

### Octant Assignment
```cpp
// Given midpoints (mi, mj, mk) and voxel (i, j, k):
int octant = (i >= mi ? 1 : 0) +
             (j >= mj ? 2 : 0) +
             (k >= mk ? 4 : 0);
// octant ∈ {0, 1, 2, 3, 4, 5, 6, 7}
```

### Score Accumulation
```cpp
// One-pass scoring
NumericVector sums(8, 0.0);
for (int v = 0; v < n; v++) {
  int oct = octants[v];
  sums[oct] += prior_vals[v] * exp(kappa * z_vals[v]);
}
// Return log(sums)
```

---

## Makevars

```makefile
PKG_CXXFLAGS = -I../inst/include
CXX_STD = CXX11
```

---

## NAMESPACE Entries

After running `Rcpp::compileAttributes()`:

```
useDynLib(neurothresh, .registration=TRUE)
importFrom(Rcpp, evalCpp)
```

---

## Testing C++ Functions

```r
# In tests/testthat/test-rcpp.R
test_that("bbox_from_indices_cpp works", {
  dims <- c(10L, 10L, 10L)
  indices <- c(1L, 10L, 100L, 500L)  # 1-based
  bbox <- bbox_from_indices_cpp(indices, dims)
  expect_length(bbox, 6)
  expect_true(all(bbox >= 0))  # 0-based bbox
})

test_that("score_children_onepass_cpp sums correctly", {
  indices <- 1:8
  octants <- 0:7  # each voxel in different octant
  z_vals <- rep(1.0, 8)
  prior_vals <- rep(1.0, 8)
  kappa <- 1.0
  scores <- score_children_onepass_cpp(indices, octants, z_vals, prior_vals, kappa)
  expect_equal(scores, rep(log(exp(1.0)), 8))  # log(e^1) = 1
})
```
