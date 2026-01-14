# Part 5: Bbox Helpers and Standard Thresholding Methods (RFT, TFCE, Cluster-FDR)

## Additional Rcpp Helpers: Bounding Box and Combined Split

These micro-optimizations matter because, unlike classic GRF/RFT peak-height inference (where the max tail is approximated analytically via expected Euler characteristic of excursion sets), our method is permutation-calibrated; so you pay for repeated set scoring and any R-level subsetting overhead shows up immediately.

---

## 1) Rcpp Exports for Bounding Box

### (A) bbox_from_indices_cpp(): Bbox for Any Mask-Space Index List

```cpp
#include <Rcpp.h>
#include <climits>
#include <cmath>

// [[Rcpp::plugins(cpp11)]]
using Rcpp::IntegerVector;

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
```

**Use case:**
- Root parcel bbox
- Any arbitrary set bbox
- Sanity checks in tests

### (B) split_indices_bbox_cpp(): Split Parent into Children and Compute Each Child Bbox

This is the production one: it avoids both R subsetting and repeated bbox computation.

```cpp
#include <Rcpp.h>
#include <climits>

// [[Rcpp::plugins(cpp11)]]
using Rcpp::IntegerVector;
using Rcpp::IntegerMatrix;
using Rcpp::List;

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
```

---

## 2) R Wrappers

```r
bbox_from_indices <- function(idx, x, y, z) {
  bbox_from_indices_cpp(idx, x, y, z)
}

split_indices_bbox <- function(idx_parent, split, x, y, z) {
  # split is the list returned by octree_split_info_cpp:
  # split$child_id is -1 or 0..m-1; split$m is m
  split_indices_bbox_cpp(idx_parent, split$child_id, split$m, x, y, z)
}
```

---

## 3) Updated Recursion: No More min(x[idx]) in R

Replace:

```r
idx_children <- split_indices_cpp(...)
bbox_child <- c(min(x[idx_child]), max(...), ...)
```

With:

```r
# After WY step-down among children:
res <- split_indices_bbox(idx_parent, split, x, y, z)
idx_children <- res$idx         # list length m
bbox_children <- res$bbox       # integer matrix m x 6

for (k in seq_along(rej)) {
  j <- rej[k]                       # child index 1..m in R
  idx_child <- idx_children[[j]]
  bbox_child <- bbox_children[j, ]  # c(x0,x1,y0,y1,z0,z1)
  hits <- descend_fast(idx_child, bbox_child, alpha_budget=..., ...)
}
```

For the root parcel bbox:

```r
idx_parent <- parcel_idx_list[[p]]            # mask-space indices
bbox <- bbox_from_indices(idx_parent, x, y, z)
```

---

## 4) What This Buys You

**Before:** For each recursion step you were doing, in R:
- Build idx_child
- Compute `min(x[idx_child])`, `max(x[idx_child])` (subsetting allocates)
- Repeat for y and z

That's 6 allocations + 6 passes per child just for bbox, per node, per recursion path.

**After:** For each split you do:
- One C++ pass through parent voxels to compute all child bboxes
- One C++ pass to split indices (already needed)

So the bbox cost is essentially "free" relative to the scoring.

---

## Extended Helper: split_indices_bbox_mass_cpp()

This helper also accumulates `den1_child[j] = sum_{v in child j} pi(v)` and `den2_child[j] = sum_{v in child j} pi(v)^2` in the same pass.

```cpp
// [[Rcpp::export]]
Rcpp::List split_indices_bbox_mass_cpp(const Rcpp::IntegerVector& idx_parent,
                                       const Rcpp::IntegerVector& child_id, // -1 or 0..m-1
                                       const int m,
                                       const Rcpp::IntegerVector& x,
                                       const Rcpp::IntegerVector& y,
                                       const Rcpp::IntegerVector& z,
                                       const Rcpp::NumericVector& pi_vec) {
  using Rcpp::IntegerVector;
  using Rcpp::IntegerMatrix;
  using Rcpp::NumericVector;
  using Rcpp::List;

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
```

### R Wrapper

```r
split_indices_bbox_mass <- function(idx_parent, split, x, y, z, pi_vec) {
  split_indices_bbox_mass_cpp(idx_parent, split$child_id, split$m, x, y, z, pi_vec)
}
```

### Updated Recursion: α Allocation with No R-Level Sums

Replace:

```r
idx_children <- split_child_indices(...)
child_mass <- vapply(idx_children, function(ii) sum(pi_vec[ii]), numeric(1))
w <- child_mass[rej] / sum(child_mass[rej])
```

With:

```r
res <- split_indices_bbox_mass(idx_parent, split, x, y, z, pi_vec)

idx_children  <- res$idx
bbox_children <- res$bbox
den1_child    <- res$den1        # sum(pi) per child
den2_child    <- res$den2        # sum(pi^2) per child (optional use)

# allocate remaining descendant alpha among rejected kids:
w <- den1_child[rej] / sum(den1_child[rej])

for (k in seq_along(rej)) {
  j <- rej[k]
  idx_child  <- idx_children[[j]]
  bbox_child <- bbox_children[j, ]   # c(x0,x1,y0,y1,z0,z1)

  hits <- descend_fast(idx_child, bbox_child,
                       alpha_budget = alpha_desc * w[k],
                       ...)
}
```

---

# Standard Thresholding Methods (Baselines)

## Shared Conventions

**Inputs:**
- `stat`: DenseNeuroVol containing a statistic field (typically Z, T, or -log10(p))
- `mask`: optional LogicalNeuroVol defining analysis domain
- `connect`: "26-connect" | "18-connect" | "6-connect" for 3D adjacency
- `tail`: "pos" | "neg" | "two"

**Clustering primitive (neuroim2):**

Use `conn_comp(stat, threshold=..., cluster_table=TRUE, local_maxima=TRUE, connect=...)`, which returns:
- `index`: ClusteredNeuroVol of cluster labels
- `size`: NeuroVol of cluster sizes
- `voxels`: list of voxel coordinates per cluster
- `cluster_table`: dataframe with cluster stats
- `local_maxima`: matrix of local maxima coords (optional)

---

## 1) Standard GRF/RFT Thresholding (SPM-Style Baseline)

This is the classic analytic approximation that controls FWER by approximating `P(max_{s∈S} Z(s) > u)` using the expected Euler characteristic (EC) of the excursion set above threshold u. The key RFT heuristic: for high thresholds, the expected EC approximates the tail probability of the maximum.

Worsley's chapter summarizes this as:
```
P(max Z > z) ≈ Σ_{d=0}^{D} Resels_d · EC_d(z)
```

### 1A) Peak-Height (Voxelwise) FWER Control

Goal: find u_α such that `P(max Z > u_α) ≈ α`, then threshold voxels with Z ≥ u_α.

**Practical note:** RFT is not "unsmoothed peak thresholding." It assumes a smooth field (often after Gaussian smoothing), and uses smoothness/FWHM (often estimated from residuals).

```r
rft_peak_fwer <- function(stat, mask=NULL, alpha=0.05,
                         fwhm_mm, # scalar or length-3
                         df=Inf, tail="pos") {
  if (is.null(mask)) mask <- as.mask(is.finite(stat))
  Z <- stat
  if (tail == "neg") Z <- -Z
  if (tail == "two") stop("use two-sided wrapper; see below")

  # 1) Estimate search geometry or LKCs/resels
  # Minimal: Resels_3 = V / FWHM^3 (stationary isotropic approx)
  V_mm3 <- voxel_volume_mm3(Z) * nvox(mask)
  Resel3 <- V_mm3 / prod(fwhm_mm)

  # 2) Define RFT p-value approximation for max statistic
  p_rft_max <- function(u) {
    # EC density for Gaussian field (D=3) or t-field
    # p ≈ Resel3 * rho3(u) + ... (lower-d terms omitted for brevity)
    Resel3 * rho3_gaussian(u)
  }

  # 3) Also compute Bonferroni bound and take min threshold if desired
  u_bonf <- qnorm(1 - alpha / nvox(mask))  # for Gaussian Z
  u_rft  <- uniroot(function(u) p_rft_max(u) - alpha, c(2, 10))$root
  u_star <- min(u_bonf, u_rft)

  sig_mask <- as.mask((Z >= u_star) & mask)
  list(u=u_star, sig_mask=sig_mask)
}

# Two-sided wrapper
rft_peak_fwer_two <- function(stat, ...) {
  res_pos <- rft_peak_fwer(stat, tail="pos", ...)
  res_neg <- rft_peak_fwer(stat, tail="neg", ...)
  sig <- as.mask(res_pos$sig_mask | res_neg$sig_mask)
  list(u_pos=res_pos$u, u_neg=res_neg$u, sig_mask=sig)
}
```

### 1B) Cluster-Extent (Clusterwise) FWER Control

Goal: threshold at some cluster-forming height z_0 (often ~3 in Gaussian Z fields), form clusters, and test clusters by their extent.

Worsley gives:
- Extent of a single cluster: `P(S > s) ≈ exp{-z_0 (s/c)^{2/D}}`
- Correction via expected number of clusters E(K): `P(max S > s) ≈ E(K) P(S > s)`

```r
rft_cluster_fwer <- function(stat, mask=NULL,
                             z0=3.0, alpha=0.05,
                             fwhm_mm, connect="26-connect",
                             tail="pos") {
  if (is.null(mask)) mask <- as.mask(is.finite(stat))
  Z <- stat
  if (tail == "neg") Z <- -Z
  if (tail == "two") stop("wrap two-sided externally")

  # 1) Find supra-threshold clusters at z0
  comps <- neuroim2::conn_comp(Z * mask, threshold=z0,
                               cluster_table=TRUE, local_maxima=TRUE,
                               connect=connect)

  tab <- comps$cluster_table
  if (nrow(tab) == 0) return(list(sig_mask=as.mask(FALSE * mask), table=tab))

  # 2) Convert cluster size from voxels -> volume or resels
  voxvol <- voxel_volume_mm3(Z)
  s_mm3 <- tab$N * voxvol
  s_resel <- s_mm3 / prod(fwhm_mm)  # cluster "resels"

  # 3) Uncorrected p for each cluster's extent under RFT model
  c_const <- rft_c_const(fwhm_mm, z0, D=3)
  p_unc <- exp(-z0 * (s_resel / c_const)^(2/3))

  # 4) Corrected cluster p-values (FWER) via expected number of clusters E(K)
  EK <- rft_expected_num_clusters(z0, mask, fwhm_mm)
  p_fwer <- pmin(1, EK * p_unc)

  tab$p_unc <- p_unc
  tab$p_fwer <- p_fwer

  # 5) Significant clusters -> voxel mask
  sig_ids <- tab$index[p_fwer <= alpha]
  sig_mask <- mask_from_cluster_ids(comps$index, sig_ids)

  list(sig_mask=sig_mask, table=tab, comps=comps)
}
```

---

## 2) FSL TFCE (Threshold-Free Cluster Enhancement)

TFCE replaces the "choose a cluster-forming threshold z_0" step with an integral over thresholds. Smith & Nichols define TFCE at voxel p as:

```
TFCE(p) = ∫_{h=h_0}^{h_p} e(h)^E · h^H dh
```

where e(h) is the extent of the cluster containing p when thresholding at level h. Typical parameters: E=0.5, H=2. Approximate the integral with steps dh (e.g., dh=0.1 for t/Z images).

For inference, TFCE is typically paired with permutation testing: build the null distribution of the maximum TFCE value across voxels and threshold the observed TFCE map at the appropriate percentile for FWER control.

### 2A) TFCE Transform (No Inference Yet)

```r
tfce_transform <- function(stat, mask=NULL,
                           E=0.5, H=2.0, dh=0.1,
                           connect="26-connect",
                           tail="pos") {
  if (is.null(mask)) mask <- as.mask(is.finite(stat))
  Z <- stat
  if (tail == "neg") Z <- -Z
  if (tail == "two") stop("call twice and combine")

  # Only positive support
  Zp <- pmax(Z, 0) * mask
  hmax <- max(Zp, na.rm=TRUE)
  hs <- seq(0, hmax, by=dh)

  TF <- zero_like(Z)  # DenseNeuroVol same space

  # Iterate thresholds
  for (h in hs) {
    comps <- neuroim2::conn_comp(Zp, threshold=h,
                                 cluster_table=TRUE,
                                 local_maxima=FALSE,
                                 connect=connect)
    tab <- comps$cluster_table
    if (nrow(tab) == 0) next

    # cluster extent e(h): use physical volume for resolution-invariance
    voxvol <- voxel_volume_mm3(Z)
    extent <- tab$N * voxvol  # e(h) in mm^3

    # For each cluster, add extent^E * h^H * dh to voxels in that cluster
    for (r in seq_len(nrow(tab))) {
      cid <- tab$index[r]
      vox_idx <- cluster_indices(comps$index, cid)
      TF[vox_idx] <- TF[vox_idx] + (extent[r]^E) * (h^H) * dh
    }
  }
  TF
}

tfce_two_sided <- function(stat, ...) {
  TF_pos <- tfce_transform(stat, tail="pos", ...)
  TF_neg <- tfce_transform(stat, tail="neg", ...)
  list(TF_pos=TF_pos, TF_neg=TF_neg)
}
```

### 2B) TFCE Inference (FWER via Permutation)

```r
tfce_fwer <- function(stat, mask, perms,
                      E=0.5, H=2, dh=0.1,
                      connect="26-connect",
                      alpha=0.05) {

  TF_obs <- tfce_transform(stat, mask, E, H, dh, connect, tail="pos")

  max_null <- numeric(length(perms))
  for (b in seq_along(perms)) {
    stat_b <- apply_permutation(stat, perms[[b]])
    TF_b <- tfce_transform(stat_b, mask, E, H, dh, connect, tail="pos")
    max_null[b] <- max(TF_b[mask], na.rm=TRUE)
  }

  thr <- quantile(max_null, probs=1-alpha, na.rm=TRUE)
  sig_mask <- as.mask((TF_obs >= thr) & mask)

  # Optionally produce voxelwise FWER-corrected p-values
  # p_corr(v) = (1 + #{b: max_null[b] >= TF_obs[v]})/(B+1)
  p_corr <- fwer_p_from_maxnull(TF_obs, max_null, mask)

  list(TF=TF_obs, thr=thr, sig_mask=sig_mask, p_corr=p_corr)
}
```

---

## 3) Spatially Aware Cluster-Based FDR

### Why Not Voxelwise FDR?

Chumbley & Friston argue that voxelwise FDR is not the same as controlling the FDR of activations (topological features like clusters/peaks), especially for smooth distributed signals.

### Recommended: Topological Cluster-FDR (Cluster-Volume FDR)

Their "simple solution": define a list of clusters above an ad hoc height threshold, compute uncorrected cluster-volume p-values under RFT, then apply an FDR procedure to that finite list of clusters.

```r
cluster_fdr_rft <- function(stat, mask=NULL,
                            z0=3.0, q=0.05,
                            fwhm_mm,
                            connect="26-connect",
                            tail="pos") {

  if (is.null(mask)) mask <- as.mask(is.finite(stat))
  Z <- stat
  if (tail == "neg") Z <- -Z
  if (tail == "two") stop("wrap two-sided externally")

  # 1) Cluster at cluster-forming threshold z0
  comps <- neuroim2::conn_comp(Z * mask, threshold=z0,
                               cluster_table=TRUE, local_maxima=TRUE,
                               connect=connect)
  tab <- comps$cluster_table
  if (nrow(tab) == 0) {
    return(list(sig_mask=as.mask(FALSE * mask), table=tab, q_threshold=NA_real_))
  }

  # 2) Compute *uncorrected* cluster-extent p-values under RFT
  voxvol <- voxel_volume_mm3(Z)
  s_mm3 <- tab$N * voxvol
  s_resel <- s_mm3 / prod(fwhm_mm)

  c_const <- rft_c_const(fwhm_mm, z0, D=3)
  p_unc <- exp(-z0 * (s_resel / c_const)^(2/3))
  tab$p_unc <- p_unc

  # 3) BH (Benjamini-Hochberg) across clusters
  m <- length(p_unc)
  o <- order(p_unc)
  p_sorted <- p_unc[o]
  thresh_line <- (seq_len(m) / m) * q
  k <- max(which(p_sorted <= thresh_line), 0L)
  p_star <- if (k == 0L) NA_real_ else p_sorted[k]

  sig <- rep(FALSE, m)
  if (k > 0L) sig[o[seq_len(k)]] <- TRUE
  tab$cluster_fdr_sig <- sig
  tab$q_level <- q

  # 4) Convert significant clusters to voxel mask
  sig_ids <- tab$index[sig]
  sig_mask <- mask_from_cluster_ids(comps$index, sig_ids)

  list(sig_mask=sig_mask, table=tab, q_threshold=p_star, comps=comps)
}
```

---

## Method Registry Interface

```r
threshold_map <- function(stat, method=c("rft_peak","rft_cluster","tfce_fwer","cluster_fdr"),
                          ...) {
  method <- match.arg(method)
  switch(method,
    rft_peak    = rft_peak_fwer(stat, ...),
    rft_cluster = rft_cluster_fwer(stat, ...),
    tfce_fwer   = tfce_fwer(stat, ...),
    cluster_fdr = cluster_fdr_rft(stat, ...)
  )
}
```

---

# Canonicalization Layer: Support Z / t(df) / -log10(p)

## Adapter Rules

- If input is Z: Zeq = Z
- If input is t(df): use probability integral transform `u = F_{t,ν}(t), Zeq = qnorm(u)` — marginally exact N(0,1) under the null
- If input is -log10(p): recover p, invert to |Z|, then multiply by a sign source. Without a sign map you cannot do directional two-sided localization.

```r
.clamp01 <- function(u) pmin(pmax(u, .Machine$double.xmin), 1 - .Machine$double.eps)

log10p_to_p <- function(log10p) {
  # p = 10^(-log10p) = exp(-log(10) * log10p)
  p <- exp(-log(10) * log10p)
  .clamp01(p)
}

canonicalize_stat <- function(stat_vol,
                              stat_type = c("Z","t","log10p"),
                              df = NULL,
                              tail = c("pos","neg","two"),
                              p_side = c("one","two"),
                              sign_vol = NULL) {

  stat_type <- match.arg(stat_type)
  tail <- match.arg(tail)
  p_side <- match.arg(p_side)

  mask <- is.finite(stat_vol)
  s <- as.numeric(stat_vol[mask])

  Zeq <- rep(NA_real_, length(s))
  p_one <- rep(NA_real_, length(s))
  p_two <- rep(NA_real_, length(s))

  if (stat_type == "Z") {
    z <- s
    ppos <- 1 - pnorm(z)
    pneg <- pnorm(z)

    if (tail == "pos") {
      p_one <- .clamp01(ppos)
      Zeq   <- z
    } else if (tail == "neg") {
      p_one <- .clamp01(pneg)
      Zeq   <- z
    } else {
      p_two <- .clamp01(2 * (1 - pnorm(abs(z))))
      Zeq <- z
    }

    if (tail != "two") {
      p_two <- .clamp01(2 * pmin(ppos, pneg))
    }

  } else if (stat_type == "t") {
    if (is.null(df)) stop("df is required for stat_type='t'")

    tval <- s
    p_two <- .clamp01(2 * (1 - pt(abs(tval), df=df)))

    ppos <- 1 - pt(tval, df=df)
    pneg <- pt(tval, df=df)
    p_one <- if (tail == "neg") .clamp01(pneg) else .clamp01(ppos)

    u <- .clamp01(pt(tval, df=df))
    Zeq <- qnorm(u)

  } else if (stat_type == "log10p") {
    lp <- s
    p_raw <- log10p_to_p(lp)

    if (p_side == "one") {
      p_one <- p_raw
      p_two <- .clamp01(2 * p_raw)
      zabs  <- qnorm(1 - p_raw)
    } else {
      p_two <- p_raw
      p_one <- .clamp01(p_raw / 2)
      zabs  <- qnorm(1 - p_raw/2)
    }

    if (!is.null(sign_vol)) {
      sg <- sign(as.numeric(sign_vol[mask]))
      sg[sg == 0] <- 1
      Zeq <- sg * zabs
    } else {
      if (tail == "pos") Zeq <- zabs
      else if (tail == "neg") Zeq <- -zabs
      else Zeq <- zabs
    }
  }

  Zeq_vol  <- stat_vol; Zeq_vol[mask]  <- Zeq
  pone_vol <- stat_vol; pone_vol[mask] <- p_one
  ptwo_vol <- stat_vol; ptwo_vol[mask] <- p_two

  list(Zeq = Zeq_vol, p_one = pone_vol, p_two = ptwo_vol, mask = mask)
}
```

---

## Two-Sided Wrappers

### For Analytic RFT/GRF Methods: Split α/2 Per Tail

```r
two_sided_split_fwer <- function(run_one_sided, stat, alpha=0.05, ...) {
  res_pos <- run_one_sided(stat, tail="pos", alpha=alpha/2, ...)
  res_neg <- run_one_sided(stat, tail="neg", alpha=alpha/2, ...)

  sig <- (res_pos$sig_mask | res_neg$sig_mask)
  list(sig_mask=sig, pos=res_pos, neg=res_neg)
}
```

### For Permutation/Max-Stat Methods: Combined Max Null (Less Conservative)

For methods that produce an "enhanced" voxelwise map (TFCE, scan statistic):
- Compute two enhanced maps: Epos = enhance(Z), Eneg = enhance(-Z)
- Define two-sided evidence map: E2 = pmax(Epos, Eneg)
- For each permutation, compute E2_b and record max(E2_b)
- Threshold E2 by the (1-α) quantile of that null max distribution

---

## TFCE with Consistent Two-Sided Inference

```r
tfce_fwer_two_sided <- function(stat, perms,
                                stat_type=c("Z","t","log10p"), df=NULL,
                                p_side="two", sign_vol=NULL,
                                E=0.5, H=2, dh=0.1,
                                connect="26-connect",
                                alpha=0.05) {

  canon <- adapt_stat(stat, stat_type, df=df, p_side=p_side, sign_vol=sign_vol)
  Z <- canon$Zeq
  mask <- canon$mask

  # observed
  TF_pos <- tfce_transform(Z,  mask=mask, E=E, H=H, dh=dh, connect=connect, tail="pos")
  TF_neg <- tfce_transform(-Z, mask=mask, E=E, H=H, dh=dh, connect=connect, tail="pos")
  TF_2   <- pmax(TF_pos, TF_neg)

  # null max distribution of the two-sided enhanced map
  max_null <- numeric(length(perms))
  for (b in seq_along(perms)) {
    stat_b <- apply_permutation(stat, perms[[b]])
    canon_b <- adapt_stat(stat_b, stat_type, df=df, p_side=p_side, sign_vol=sign_vol)
    Zb <- canon_b$Zeq

    TFp <- tfce_transform(Zb,  mask=mask, E=E, H=H, dh=dh, connect=connect, tail="pos")
    TFn <- tfce_transform(-Zb, mask=mask, E=E, H=H, dh=dh, connect=connect, tail="pos")
    max_null[b] <- max(pmax(TFp, TFn)[mask], na.rm=TRUE)
  }

  thr <- unname(quantile(max_null, probs=1-alpha, na.rm=TRUE))

  sig_pos <- (TF_pos >= thr) & (TF_pos >= TF_neg) & mask
  sig_neg <- (TF_neg >= thr) & (TF_neg >  TF_pos) & mask
  sig_any <- (TF_2   >= thr) & mask

  list(TF_pos=TF_pos, TF_neg=TF_neg, TF_2=TF_2,
       thr=thr, sig_mask=sig_any, sig_pos=sig_pos, sig_neg=sig_neg,
       max_null=max_null)
}
```

---

## Cluster-FDR with Parametric-or-Permutation Switch

### Cluster Statistic Extractor

```r
cluster_stats <- function(Z, mask, u0, connect="26-connect", stat=c("extent","mass")) {
  stat <- match.arg(stat)
  comps <- neuroim2::conn_comp(Z * mask, threshold=u0,
                               cluster_table=TRUE,
                               local_maxima=TRUE,
                               connect=connect)

  tab <- comps$cluster_table
  if (nrow(tab) == 0) return(list(comps=comps, tab=tab, s=numeric(0)))

  if (stat == "extent") {
    voxvol <- voxel_volume_mm3(Z)
    s <- tab$N * voxvol
  } else {
    # cluster mass: sum(Z-u0) inside each cluster
    s <- vapply(tab$index, function(cid) {
      idx <- cluster_indices(comps$index, cid)
      sum((Z[idx] - u0), na.rm=TRUE)
    }, numeric(1))
  }
  list(comps=comps, tab=tab, s=s)
}
```

### Parametric RFT Cluster P-Values

```r
cluster_pvals_rft <- function(s, u0, fwhm_mm, D=3) {
  # Approximate c from Worsley eqns (3)-(5) (Gaussian case)
  # c = FWHM^D (2π/u0)^{D/2} (4 log 2)^{-D/2} / Γ(D/2+1)
  c0 <- prod(fwhm_mm)^1 * (2*pi/u0)^(D/2) * (4*log(2))^(-D/2) / gamma(D/2 + 1)

  # p_unc(s) ≈ exp(-u0 * (s/c)^(2/D))  (single cluster)
  p_unc <- exp(-u0 * (s / c0)^(2/D))
  pmin(pmax(p_unc, 0), 1)
}
```

### Nonparametric Permutation Cluster P-Values

```r
cluster_pvals_perm <- function(stat, perms,
                               stat_type, df=NULL, p_side="two", sign_vol=NULL,
                               u0, connect="26-connect",
                               stat_fun=c("extent","mass"),
                               tail=c("pos","neg")) {

  tail <- match.arg(tail)
  stat_fun <- match.arg(stat_fun)

  canon <- adapt_stat(stat, stat_type, df=df, p_side=p_side, sign_vol=sign_vol)
  Z <- canon$Zeq
  mask <- canon$mask
  if (tail == "neg") Z <- -Z

  obs <- cluster_stats(Z, mask, u0, connect=connect, stat=stat_fun)
  s_obs <- obs$s
  if (length(s_obs) == 0) return(list(p_unc=numeric(0), obs=obs, null_stats=numeric(0)))

  null_stats <- numeric(0)

  for (b in seq_along(perms)) {
    stat_b <- apply_permutation(stat, perms[[b]])
    canon_b <- adapt_stat(stat_b, stat_type, df=df, p_side=p_side, sign_vol=sign_vol)
    Zb <- canon_b$Zeq
    if (tail == "neg") Zb <- -Zb

    tmp <- cluster_stats(Zb, mask, u0, connect=connect, stat=stat_fun)
    if (length(tmp$s) > 0) null_stats <- c(null_stats, tmp$s)
  }

  if (length(null_stats) == 0) {
    p_unc <- rep(1/(length(perms)+1), length(s_obs))
  } else {
    p_unc <- vapply(s_obs, function(s) {
      (1 + sum(null_stats >= s)) / (1 + length(null_stats))
    }, numeric(1))
  }

  list(p_unc=p_unc, obs=obs, null_stats=null_stats)
}
```

### Cluster-FDR Driver with method = "rft" | "perm"

```r
bh_reject <- function(p, q) {
  m <- length(p)
  if (m == 0) return(rep(FALSE, 0))
  o <- order(p)
  ps <- p[o]
  thr <- (seq_len(m)/m) * q
  k <- max(which(ps <= thr), 0L)
  rej <- rep(FALSE, m)
  if (k > 0) rej[o[seq_len(k)]] <- TRUE
  rej
}

cluster_fdr <- function(stat, perms=NULL,
                        stat_type=c("Z","t","log10p"), df=NULL,
                        p_side="two", sign_vol=NULL,
                        p0=0.001, u0=NULL,
                        q=0.05,
                        connect="26-connect",
                        cluster_stat=c("extent","mass"),
                        p_method=c("rft","perm"),
                        fwhm_mm=NULL,
                        two_sided=TRUE,
                        two_sided_policy=c("BH_all","split_q")) {

  stat_type <- match.arg(stat_type)
  p_method <- match.arg(p_method)
  cluster_stat <- match.arg(cluster_stat)
  two_sided_policy <- match.arg(two_sided_policy)

  canon <- adapt_stat(stat, stat_type, df=df, p_side=p_side, sign_vol=sign_vol)
  Z <- canon$Zeq
  mask <- canon$mask
  if (is.null(u0)) u0 <- qnorm(1 - p0)

  run_tail <- function(tail) {
    if (p_method == "rft") {
      if (is.null(fwhm_mm)) stop("fwhm_mm required for p_method='rft'")
      Zt <- if (tail=="pos") Z else -Z
      obs <- cluster_stats(Zt, mask, u0, connect, cluster_stat)
      p_unc <- cluster_pvals_rft(obs$s, u0=u0, fwhm_mm=fwhm_mm, D=3)
      list(obs=obs, p_unc=p_unc)
    } else {
      if (is.null(perms)) stop("perms required for p_method='perm'")
      out <- cluster_pvals_perm(stat, perms, stat_type, df=df, p_side=p_side, sign_vol=sign_vol,
                               u0=u0, connect=connect, stat_fun=cluster_stat, tail=tail)
      list(obs=out$obs, p_unc=out$p_unc, null_stats=out$null_stats)
    }
  }

  if (!two_sided) {
    out <- run_tail("pos")
    rej <- bh_reject(out$p_unc, q)
    sig_ids <- out$obs$tab$index[rej]
    sig_mask <- mask_from_cluster_ids(out$obs$comps$index, sig_ids)
    out$obs$tab$p_unc <- out$p_unc
    out$obs$tab$reject <- rej
    return(list(sig_mask=sig_mask, table=out$obs$tab, comps=out$obs$comps))
  }

  pos <- run_tail("pos")
  neg <- run_tail("neg")

  if (two_sided_policy == "split_q") {
    rej_pos <- bh_reject(pos$p_unc, q/2)
    rej_neg <- bh_reject(neg$p_unc, q/2)
  } else {
    p_all <- c(pos$p_unc, neg$p_unc)
    rej_all <- bh_reject(p_all, q)
    rej_pos <- rej_all[seq_along(pos$p_unc)]
    rej_neg <- rej_all[length(pos$p_unc) + seq_along(neg$p_unc)]
  }

  sig_pos_ids <- pos$obs$tab$index[rej_pos]
  sig_neg_ids <- neg$obs$tab$index[rej_neg]

  sig_pos <- mask_from_cluster_ids(pos$obs$comps$index, sig_pos_ids)
  sig_neg <- mask_from_cluster_ids(neg$obs$comps$index, sig_neg_ids)
  sig_any <- sig_pos | sig_neg

  pos$obs$tab$p_unc <- pos$p_unc; pos$obs$tab$reject <- rej_pos; pos$obs$tab$tail <- "pos"
  neg$obs$tab$p_unc <- neg$p_unc; neg$obs$tab$reject <- rej_neg; neg$obs$tab$tail <- "neg"
  tab <- rbind(pos$obs$tab, neg$obs$tab)

  list(sig_mask=sig_any, sig_pos=sig_pos, sig_neg=sig_neg,
       table=tab,
       comps_pos=pos$obs$comps, comps_neg=neg$obs$comps,
       u0=u0, q=q, two_sided_policy=two_sided_policy, p_method=p_method)
}
```

---

## Summary: What You Now Have

- **One adapter** to support input as Z / t(df) / -log10(p), yielding a signed Zeq map
- **Consistent two-sided inference:**
  - Analytic RFT: α/2 split (standard)
  - Permutation enhanced-map methods (TFCE): combined max null (better power)
- **Cluster-FDR** that is truly "spatially aware":
  - Parametric cluster p-values (RFT; Worsley)
  - Or permutation cluster p-values (pooled null clusters)
  - BH across clusters (feature-level inference; aligns with topological FDR motivation)
- **TFCE implementation** anchored to Smith & Nichols and FSL's practice
