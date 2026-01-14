# Part 4: One-Pass Sibling Scoring and Production Rcpp Implementation

## Performance Upgrade: Score All Sibling Children in One Pass

The key idea: avoid slicing/extracting `Z[idx_child]` over and over for each permutation.

This is the same basic idea as why RFT peak inference is fast: it avoids repeated scanning by working with global objects (max/EC formulas). We can't avoid permutations here, but we can avoid a ton of R-level allocations and loops.

---

## Why the Current Sibling Scoring Is Slow

In the earlier pseudocode, at a node split you did something like:
- for each permutation b
- for each child c ∈ {1,...,m} (≤8)
- extract `Zb[idx_child]`, `pi[idx_child]`
- compute `U0(child) + max_kappa Softmax(child,kappa)`

This incurs:
- Repeated vector allocations (`Z[idx_child]`)
- Repeated weight extraction (`pi[idx_child]`)
- Repeated loops at the R interpreter level

Even though total arithmetic is proportional to the parent's voxel count, the constant factors are brutal.

---

## The Key Idea: Split Once, Then Accumulate Sufficient Statistics by Child Label

At split time, create a child label vector aligned with the parent's voxel list:
- Parent has mask-space indices `idx_parent` (length L)
- Define `child_id[t] ∈ {1..m}` telling you which child voxel `idx_parent[t]` belongs to

Then for any permutation map Z_b, you can compute all children's scores by scanning the parent voxels once (per κ) and updating accumulators for the appropriate child.

### What You Need to Precompute/Carry

For a given parent node:
- `idx_parent`: mask-space indices into Z_vec / pi_vec
- `w_parent`: `pi_vec[idx_parent]` (cached once)
- `child_id`: integer vector length L
- For each child j (1..m):
  - `log_den1[j] = log(sum(pi in child j))`
  - `inv_sqrt_den2[j] = 1 / sqrt(sum(pi^2 in child j))`

These last two are constants for that node (already computed from pi_mass and pi2_mass).

---

## One-Pass Sibling Scoring: κ=0 and κ>0

### κ=0 (Diffuse) Branch: Standardized Mean U_0

You want:
```
U_0(R) = Σ_{v∈R} π(v)·Z(v) / √(Σ_{v∈R} π(v)²)
```

For children, just accumulate:
- `num0[j] += w * z` for child j

Then:
- `U0[j] = num0[j] * inv_sqrt_den2[j]`

### κ>0 (Focal) Branch: Weighted Log-Sum-Exp S_κ

You want:
```
S_κ(R) = (1/κ)·[log(Σ_{v∈R} π(v)·e^{κ·Z(v)}) - log(Σ_{v∈R} π(v))]
```

For numerical stability and streaming (no storing vectors), maintain per child and per κ:
- `maxa[j,k]` (current max of a=κ*z)
- `sumexp[j,k]` (current sum of w*exp(a-maxa))

**Streaming update for one observation (w, a):**
```
if a > maxa:
  sumexp = sumexp * exp(maxa - a) + w
  maxa = a
else:
  sumexp += w * exp(a - maxa)
```

**End:**
```
logsum = log(sumexp) + maxa
S = (logsum - log_den1[j]) / κ
```

Then:
- `best_soft[j] = max_k S(j,k)`

**Final child score:**
- `T_child[j] = max(U0[j], best_soft[j])`

---

## Where to Plug This Into the Hierarchy

In `descend_node()`, instead of computing `t_perm` by looping children and calling `score_node()` repeatedly, do:
- Pre-split parent → store `child_id`, `w_parent`, and child constants
- For permutation b: compute `t_perm[b, ] = score_children_onepass(Zb_vec, parent_split_info, kappa_grid)`

Because m ≤ 8, this is extremely cheap (and step-down over siblings is trivial).

---

## R-Flavored Pseudocode

### 1) Build Split Info Once (Augment the Parent After Splitting)

```r
build_split_info <- function(parent, pi_vec, min_pi_mass=1e-10) {
  # parent$idx is mask-space indices; parent$grid are ijk for those idx.
  idx <- parent$idx
  g   <- parent$grid
  w_parent <- pi_vec[idx]              # cache weights

  # split plane
  x0 <- parent$bbox[1]; x1 <- parent$bbox[2]
  y0 <- parent$bbox[3]; y1 <- parent$bbox[4]
  z0 <- parent$bbox[5]; z1 <- parent$bbox[6]
  xm <- floor((x0 + x1)/2); ym <- floor((y0 + y1)/2); zm <- floor((z0 + z1)/2)

  # child_id in {1..8} based on which side of each midpoint
  # (this encodes the same dyadic octree boxes)
  child_id <- 1L +
    (g[,1] > xm) * 1L +
    (g[,2] > ym) * 2L +
    (g[,3] > zm) * 4L

  # compute child den1, den2 once (from w_parent and child_id)
  # (in real code: use fast tabulation; shown plainly here)
  den1 <- numeric(8); den2 <- numeric(8)
  for (t in seq_along(idx)) {
    j <- child_id[t]
    wt <- w_parent[t]
    den1[j] <- den1[j] + wt
    den2[j] <- den2[j] + wt*wt
  }

  # keep only non-empty kids (optional)
  keep <- which(den1 >= min_pi_mass)
  # remap child ids to 1..m for compactness
  map <- rep(0L, 8); map[keep] <- seq_along(keep)
  child_id2 <- map[child_id]
  m <- length(keep)

  list(
    idx_parent = idx,
    w_parent   = w_parent,
    child_id   = child_id2,      # now in {1..m}
    log_den1   = log(den1[keep]),
    inv_sqrt_den2 = 1/sqrt(den2[keep]),
    m = m
  )
}
```

### 2) Score All Children for One Permutation, One Pass Per κ

```r
score_children_onepass <- function(Z_vec, split, kappa_grid_pos) {
  # Z_vec: numeric vector over mask voxels (already abs() if two-sided)
  idx <- split$idx_parent
  w   <- split$w_parent
  cid <- split$child_id
  m <- split$m

  # κ=0 accumulators
  num0 <- numeric(m)

  # κ>0 accumulators: matrices m x K
  K <- length(kappa_grid_pos)
  maxa   <- matrix(-Inf, nrow=m, ncol=K)
  sumexp <- matrix(0.0, nrow=m, ncol=K)

  for (t in seq_along(idx)) {
    j <- cid[t]
    if (j == 0L) next
    z <- Z_vec[idx[t]]
    wt <- w[t]

    # κ=0
    num0[j] <- num0[j] + wt * z

    # κ>0
    for (kk in 1:K) {
      kappa <- kappa_grid_pos[kk]
      a <- kappa * z
      if (a > maxa[j,kk]) {
        sumexp[j,kk] <- sumexp[j,kk] * exp(maxa[j,kk] - a) + wt
        maxa[j,kk] <- a
      } else {
        sumexp[j,kk] <- sumexp[j,kk] + wt * exp(a - maxa[j,kk])
      }
    }
  }

  # compute child scores
  U0 <- num0 * split$inv_sqrt_den2

  best_soft <- rep(-Inf, m)
  for (kk in 1:K) {
    kappa <- kappa_grid_pos[kk]
    logsum <- log(sumexp[,kk]) + maxa[,kk]
    S <- (logsum - split$log_den1) / kappa
    best_soft <- pmax(best_soft, S)
  }

  pmax(U0, best_soft)
}
```

**Important note:** This R version is still loop-heavy; it's meant as clear pseudocode. In production, you'll want this in Rcpp (below) to eliminate per-voxel R overhead.

---

## Rcpp-Style Pseudocode (What You Actually Ship)

This is the same computation, but performed in C++ so you don't allocate `Z[idx]` sub-vectors.

```cpp
// PSEUDOCODE (Rcpp style)
NumericVector score_children_onepass_cpp(
    NumericVector Z_vec,          // length Nmask
    IntegerVector idx_parent,      // length L (mask-space indices)
    NumericVector w_parent,        // length L, pi weights for those indices
    IntegerVector child_id,        // length L, values 0..m-1
    NumericVector log_den1,        // length m
    NumericVector inv_sqrt_den2,   // length m
    NumericVector kappa_grid_pos   // length K (>0)
) {
  int L = idx_parent.size();
  int m = log_den1.size();
  int K = kappa_grid_pos.size();

  // κ=0
  std::vector<double> num0(m, 0.0);

  // κ>0 accumulators
  std::vector<double> maxa(m*K, -INFINITY);
  std::vector<double> sumexp(m*K, 0.0);

  for (int t=0; t<L; t++) {
    int j = child_id[t];           // 0..m-1, or -1 for "drop"
    if (j < 0) continue;

    int ix = idx_parent[t] - 1;    // R->C++ index shift
    double z = Z_vec[ix];
    double w = w_parent[t];

    // κ=0
    num0[j] += w * z;

    // κ>0
    for (int kk=0; kk<K; kk++) {
      double kappa = kappa_grid_pos[kk];
      double a = kappa * z;
      int pos = j + m*kk;

      double mx = maxa[pos];
      if (a > mx) {
        sumexp[pos] = sumexp[pos] * std::exp(mx - a) + w;
        maxa[pos] = a;
      } else {
        sumexp[pos] += w * std::exp(a - mx);
      }
    }
  }

  NumericVector out(m);

  for (int j=0; j<m; j++) {
    double U0 = num0[j] * inv_sqrt_den2[j];

    double best_soft = -INFINITY;
    for (int kk=0; kk<K; kk++) {
      double kappa = kappa_grid_pos[kk];
      int pos = j + m*kk;
      double logsum = std::log(sumexp[pos]) + maxa[pos];
      double S = (logsum - log_den1[j]) / kappa;
      if (S > best_soft) best_soft = S;
    }

    out[j] = std::max(U0, best_soft);
  }

  return out;
}
```

This turns the sibling scoring cost into essentially:
```
O(L · K)
```
with tiny constants (because m ≤ 8).

---

## How This Changes the Runtime Profile

With this change:
- **Root parcel stage** is still dominated by rowsum/grouping across all voxels and permutations (your #1 bottleneck)
- **Recursive octree stage** becomes much cheaper because:
  - You don't do 8× extraction/allocations per permutation per node
  - You stream once through parent voxels and update all child stats

In practice this is often the difference between "recursion is negligible" vs "recursion doubles the runtime."

---

## Practical Implementation Tip: Have perm_fun(b) Return a Mask-Space Vector

If `perm_fun(b)` currently returns a `DenseNeuroVol`, you then do `Zb_vec <- abs(Zb[mask_idx])` per permutation, which allocates. It's much faster if:
- You compute the permuted group statistic already as a numeric vector over mask voxels, or
- You provide a `perm_fun_vec(b)` wrapper that extracts only mask voxels once

This pairs perfectly with the Rcpp `score_children_onepass_cpp()` signature above.

---

## Combined Production Blueprint

- **Root level (parcels):** vectorized group scoring (one rowsum / sparse-matmul per κ per permutation)
- **Recursive level (octree children):** one-pass sibling scoring (single scan through parent voxels; implemented in Rcpp)

This is the closest you can get to "RFT-speed" while still doing permutation-based FWER.

---

## 1) Unify the Data Layout: Work in Mask-Space Vectors

Everything gets dramatically faster if you do all scoring in a compact vector space of masked voxels:

```r
mask_idx <- which(mask)  # full-volume linear indices inside mask
N <- length(mask_idx)
```

Convert DenseNeuroVol maps to mask-space numeric vectors:

```r
Z_obs_vec <- as.numeric(Z_obs[mask_idx])   # Z_obs is DenseNeuroVol
pi_vec    <- as.numeric(pi[mask_idx])      # prior map DenseNeuroVol
pi2_vec   <- pi_vec * pi_vec
```

Also precompute mask-space coordinates once (for octree splitting without building per-node grids):

```r
grid <- index_to_grid(Z_obs, mask_idx)     # N x 3 integer coords (i,j,k)  [neuroim2]
x <- grid[,1]; y <- grid[,2]; z <- grid[,3]
```

Now every node will store indices into 1..N (mask-space), and every permutation returns a Zb_vec of length N.

---

## 2) Root Parcels: Vectorized Scoring Over All Parcels (Fast Path)

### 2.1 Precompute Parcel Membership in Mask-Space

```r
lab_vec <- as.integer(atlas_labels[mask_idx])  # same mask ordering as Z_obs_vec
# Optionally drop label 0 (background) or keep as "remainder parcel"
```

Make sure parcel IDs are compact 1..P (relabel if needed).

### 2.2 Precompute Parcel Constants Once (Caching π and π²)

```r
den1 <- rowsum(pi_vec,  group=lab_vec, reorder=FALSE)   # P x 1
den2 <- rowsum(pi2_vec, group=lab_vec, reorder=FALSE)   # P x 1

log_den1      <- log(den1)
inv_sqrt_den2 <- 1 / sqrt(den2)
```

### 2.3 Vectorized Parcel Score Function (κ=0 Stabilized + κ>0 Softmax)

```r
score_parcels_vec <- function(Z_vec, pi_vec, lab_vec,
                              log_den1, inv_sqrt_den2,
                              kappa_grid_pos) {
  # κ=0: U0(p) = sum(pi*Z)/sqrt(sum(pi^2))
  num0 <- rowsum(pi_vec * Z_vec, group=lab_vec, reorder=FALSE)
  U0   <- num0 * inv_sqrt_den2

  best_soft <- rep(-Inf, length(U0))

  for (kappa in kappa_grid_pos) {
    # stable log-sum-exp with global max
    a <- kappa * Z_vec
    amax <- max(a)
    numk <- rowsum(pi_vec * exp(a - amax), group=lab_vec, reorder=FALSE)
    Sk   <- (amax + log(numk) - log_den1) / kappa
    best_soft <- pmax(best_soft, Sk)
  }

  pmax(U0, best_soft)
}
```

### Parcel Permutation Matrix (Feasible)

Because P is small (100-400), you can store `t_perm_par` (B×P) and do WY step-down at parcel level.

---

## 3) Recursive Octree: One-Pass Sibling Scoring (Fast Path)

### 3.1 Build a Split Descriptor Once Per Parent Node

Compute `child_id` without building a per-node grid; use the global x,y,z arrays:

```r
build_split_info <- function(idx_parent, bbox, pi_vec, x, y, z, min_pi_mass=1e-10) {
  # bbox = c(x0,x1,y0,y1,z0,z1) inclusive voxel coords
  x0 <- bbox[1]; x1 <- bbox[2]
  y0 <- bbox[3]; y1 <- bbox[4]
  z0 <- bbox[5]; z1 <- bbox[6]
  xm <- floor((x0 + x1)/2)
  ym <- floor((y0 + y1)/2)
  zm <- floor((z0 + z1)/2)

  # child in 1..8 using bit pattern
  cid8 <- 1L +
    (x[idx_parent] > xm) * 1L +
    (y[idx_parent] > ym) * 2L +
    (z[idx_parent] > zm) * 4L

  w_parent <- pi_vec[idx_parent]

  # accumulate den1 and den2 for the 8 children in one scan
  den1_8 <- numeric(8); den2_8 <- numeric(8)
  for (t in seq_along(idx_parent)) {
    j <- cid8[t]
    wt <- w_parent[t]
    den1_8[j] <- den1_8[j] + wt
    den2_8[j] <- den2_8[j] + wt*wt
  }

  keep <- which(den1_8 >= min_pi_mass)
  if (length(keep) == 0) return(NULL)

  # remap 1..8 -> 0..m-1 for compact child ids
  map <- rep(-1L, 8); map[keep] <- 0:(length(keep)-1L)
  child_id <- map[cid8]   # length L, values -1..m-1

  list(
    idx_parent = idx_parent,
    w_parent   = w_parent,
    child_id   = child_id,
    log_den1   = log(den1_8[keep]),
    inv_sqrt_den2 = 1 / sqrt(den2_8[keep]),
    keep8 = keep
  )
}
```

---

## 4) Full Algorithm Skeleton (Combined)

```r
run_hier_scan_fast <- function(Z_obs, mask, pi_raw, atlas_labels,
                               perm_fun_vec, B,
                               alpha=0.05,
                               kappa_grid_pos=c(0.5, 1, 2),
                               eta=0.9,
                               gamma_root=0.5,
                               gamma=0.5,
                               min_alpha=1e-6) {

  # ---- (A) Prep mask-space vectors ----
  mask_idx <- which(mask)
  Z_obs_vec <- abs(as.numeric(Z_obs[mask_idx]))
  pi <- prep_prior(pi_raw, mask, eta)              # as earlier
  pi_vec <- as.numeric(pi[mask_idx])
  pi2_vec <- pi_vec * pi_vec

  # coords for octree splitting (mask-space)
  grid <- index_to_grid(Z_obs, mask_idx)
  x <- grid[,1]; y <- grid[,2]; z <- grid[,3]

  # ---- (B) Parcel root scoring: vectorized ----
  lab_vec <- as.integer(atlas_labels[mask_idx])
  # relabel to 1..P if needed

  den1 <- rowsum(pi_vec,  group=lab_vec, reorder=FALSE)
  den2 <- rowsum(pi2_vec, group=lab_vec, reorder=FALSE)
  log_den1 <- log(den1)
  inv_sqrt_den2 <- 1 / sqrt(den2)

  t_obs_par <- score_parcels_vec(Z_obs_vec, pi_vec, lab_vec,
                                 log_den1, inv_sqrt_den2,
                                 kappa_grid_pos)

  # permutations for parcels (store BxP matrix)
  P <- length(t_obs_par)
  t_perm_par <- matrix(NA_real_, nrow=B, ncol=P)
  for (b in 1:B) {
    Zb_vec <- abs(perm_fun_vec(b))                 # already mask-space length N
    t_perm_par[b,] <- score_parcels_vec(Zb_vec, pi_vec, lab_vec,
                                        log_den1, inv_sqrt_den2,
                                        kappa_grid_pos)
  }

  # WY step-down at parcels
  sd_par <- wy_stepdown_maxT(t_obs_par, t_perm_par, alpha=gamma_root * alpha)
  sig_par <- which(sd_par$rejected)
  if (length(sig_par) == 0) return(list(hits=list()))

  # allocate remaining alpha among significant parcels by prior mass
  alpha_desc_total <- alpha - gamma_root * alpha
  parcel_pi_mass <- as.numeric(den1)               # same order as parcels in rowsum output
  w_par <- parcel_pi_mass[sig_par] / sum(parcel_pi_mass[sig_par])

  hits <- list()  # store discoveries

  # ---- (C) Recursive descent inside each significant parcel ----
  for (j in seq_along(sig_par)) {
    p <- sig_par[j]
    alpha_budget <- alpha_desc_total * w_par[j]

    # parent idx = all voxels in parcel p, in mask-space
    idx_parent <- which(lab_vec == p)

    # parent bbox from coords x,y,z restricted to idx_parent
    bbox <- c(min(x[idx_parent]), max(x[idx_parent]),
              min(y[idx_parent]), max(y[idx_parent]),
              min(z[idx_parent]), max(z[idx_parent]))

    hits <- descend_fast(idx_parent, bbox, alpha_budget,
                         Z_obs_vec, pi_vec, x,y,z,
                         perm_fun_vec, B, kappa_grid_pos,
                         gamma=gamma, min_alpha=min_alpha,
                         hits=hits)
  }

  list(hits=hits, parcel_stepdown=sd_par)
}
```

### Recursive Function: One-Pass Sibling Scoring + WY Step-Down

```r
descend_fast <- function(idx_parent, bbox, alpha_budget,
                         Z_obs_vec, pi_vec, x,y,z,
                         perm_fun_vec, B, kappa_grid_pos,
                         gamma=0.5, min_alpha=1e-6,
                         hits) {

  if (alpha_budget < min_alpha) return(hits)

  split <- build_split_info(idx_parent, bbox, pi_vec, x,y,z)
  if (is.null(split)) return(hits)

  m <- length(split$log_den1)      # number of non-empty children
  alpha_test <- gamma * alpha_budget
  alpha_desc <- alpha_budget - alpha_test

  # observed child scores (one-pass, ideally via the same C++ routine)
  t_obs_child <- score_children_onepass_cpp(
                   Z_obs_vec,
                   split$idx_parent,
                   split$w_parent,
                   split$child_id,
                   split$log_den1,
                   split$inv_sqrt_den2,
                   kappa_grid_pos)

  # permutation child scores: B x m matrix (m <= 8 so small)
  t_perm_child <- matrix(NA_real_, nrow=B, ncol=m)
  for (b in 1:B) {
    Zb_vec <- abs(perm_fun_vec(b))
    t_perm_child[b,] <- score_children_onepass_cpp(
                          Zb_vec,
                          split$idx_parent,
                          split$w_parent,
                          split$child_id,
                          split$log_den1,
                          split$inv_sqrt_den2,
                          kappa_grid_pos)
  }

  sd <- wy_stepdown_maxT(t_obs_child, t_perm_child, alpha=alpha_test)
  rej <- which(sd$rejected)
  if (length(rej) == 0) return(hits)

  # record rejected children as discoveries at this scale
  for (i in rej) {
    hits[[length(hits)+1]] <- list(
      type="cube",
      level="child_of_current",
      score=t_obs_child[i],
      p_adj_sibling=sd$p_adj[i],
      alpha_test=alpha_test
    )
  }

  # build child nodes for recursion (needs child idx lists + bbox)
  # (production: do this split in C++ to avoid R allocations)
  child_idx_list <- split_indices_cpp(split$idx_parent, split$child_id, m)

  # allocate descendant alpha proportional to prior mass in child
  child_mass <- vapply(child_idx_list, function(ii) sum(pi_vec[ii]), numeric(1))
  w <- child_mass[rej] / sum(child_mass[rej])

  for (k in seq_along(rej)) {
    ii <- rej[k]
    idx_child <- child_idx_list[[ii]]
    if (length(idx_child) == 0) next

    bbox_child <- c(min(x[idx_child]), max(x[idx_child]),
                    min(y[idx_child]), max(y[idx_child]),
                    min(z[idx_child]), max(z[idx_child]))

    hits <- descend_fast(idx_child, bbox_child,
                         alpha_budget = alpha_desc * w[k],
                         Z_obs_vec, pi_vec, x,y,z,
                         perm_fun_vec, B, kappa_grid_pos,
                         gamma=gamma, min_alpha=min_alpha,
                         hits=hits)
  }

  hits
}
```

---

## 5) What You Gain

### Root Stage: O(B·N·K) but with a Small Constant

- You're doing exactly what RFT tries to avoid (simulate the max tail), but:
  - Only at the parcel level (P small)
  - With vectorized operations (rowsum / sparse MV)

### Recursive Stage: O(B·L·K) per Visited Node, but One Pass

- For a parent with L voxels and m≤8 kids:
  - Naïve approach is ~m times slower due to repeated slicing/allocations
  - One-pass approach is essentially the theoretical minimum: stream once and update all children

This is the key engineering insight: you're using the same "scan once, compute many hypotheses" logic that makes max-statistic methods like RFT computationally attractive (even though the inference calibration differs).

---

## 6) Two Final Production Tips

1. **Make `perm_fun_vec(b)` return a numeric vector already in mask-space**
   - Avoid `DenseNeuroVol -> [mask_idx]` extraction inside the permutation loop

2. **Parallelize over permutations**
   - The `for (b in 1:B)` loops are embarrassingly parallel
   - Use `future.apply`, `parallel`, or `RcppParallel`

---

## Production Rcpp Package Implementation

### Package Wiring (Minimal)

**DESCRIPTION:**
```
Imports: Rcpp
LinkingTo: Rcpp
SystemRequirements: C++11
```

**NAMESPACE:**
```
useDynLib(yourpkg, .registration=TRUE)
importFrom(Rcpp, evalCpp)
```

**src/Makevars (optional but recommended):**
```
CXX_STD = CXX11
PKG_CXXFLAGS += -O3
```

**One-time setup in R:**
```r
Rcpp::compileAttributes()
```

---

## C++ File: src/hier_scan_rcpp.cpp

Three exported functions:
1. `octree_split_info_cpp()` - Builds child assignment and caches child constants
2. `score_children_onepass_cpp()` - Computes all children's omnibus scores in one pass
3. `split_indices_cpp()` - Converts (idx_parent, child_id) into list for recursion

```cpp
// src/hier_scan_rcpp.cpp
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>

// [[Rcpp::plugins(cpp11)]]

using Rcpp::IntegerVector;
using Rcpp::NumericVector;
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

  // κ=0 accumulators
  std::vector<double> num0(m, 0.0);

  // κ>0 streaming log-sum-exp accumulators
  const double NEG_INF = -std::numeric_limits<double>::infinity();
  std::vector<double> maxa(m * K, NEG_INF);
  std::vector<double> sumexp(m * K, 0.0);

  for (int t = 0; t < L; t++) {
    const int j = child_id[t];
    if (j < 0 || j >= m) continue;

    const int ix = idx_parent[t] - 1;
    double z = Z_vec[ix];
    if (Rcpp::NumericVector::is_na(z)) continue;
    if (do_abs) z = std::fabs(z);

    const double w = w_parent[t];

    // κ=0
    num0[j] += w * z;

    // κ>0
    for (int kk = 0; kk < K; kk++) {
      const double kappa = kappa_grid_pos[kk];
      const double a = kappa * z;
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
```

---

## R Wrappers (Minimal)

Create `R/rcpp_wrappers.R`:

```r
#' Build octree split info in mask-space
octree_split_info <- function(idx_parent, bbox, x, y, z, pi_vec, min_pi_mass=1e-10) {
  octree_split_info_cpp(idx_parent, bbox, x, y, z, pi_vec, min_pi_mass)
}

#' Score all children in one pass (one permutation)
score_children_onepass <- function(Z_vec, idx_parent, split, kappa_grid_pos, do_abs=FALSE) {
  score_children_onepass_cpp(
    Z_vec = Z_vec,
    idx_parent = idx_parent,
    w_parent = split$w_parent,
    child_id = split$child_id,
    log_den1 = split$log_den1,
    inv_sqrt_den2 = split$inv_sqrt_den2,
    kappa_grid_pos = kappa_grid_pos,
    do_abs = do_abs
  )
}

#' Split indices for recursion
split_child_indices <- function(idx_parent, split) {
  split_indices_cpp(idx_parent, split$child_id, split$m)
}
```

---

## Usage in Hierarchical Step-Down Recursion

Inside your `descend_fast()` loop, replace R scoring with:

```r
split <- octree_split_info(idx_parent, bbox, x, y, z, pi_vec)

# observed
t_obs_child <- score_children_onepass(Z_obs_vec, idx_parent, split, kappa_grid_pos)

# permutations
t_perm_child <- matrix(NA_real_, nrow=B, ncol=split$m)
for (b in 1:B) {
  Zb_vec <- perm_fun_vec(b)     # IMPORTANT: already mask-space numeric vector length N
  t_perm_child[b, ] <- score_children_onepass(Zb_vec, idx_parent, split, kappa_grid_pos)
}

sd <- wy_stepdown_maxT(t_obs_child, t_perm_child, alpha=alpha_test)

# for rejected children, get idx lists to recurse
idx_children <- split_child_indices(idx_parent, split)
```

That's it: no `Z[idx_child]` allocations, no loops over children in R, and the only loop left in R is over permutations (which you can parallelize later).

---

## Sanity Checklist

- **Indexing convention:** `idx_parent` must be mask-space indices into Z_vec/pi_vec/x/y/z, 1-based (R style). The C++ code subtracts 1 when reading from vectors.
- **`perm_fun_vec(b)`** should return mask-space vectors. Don't return DenseNeuroVol then subscript; that will allocate a lot.
- **`kappa_grid_pos`** must be strictly > 0. (κ=0 handled via the U0 branch inside the scorer.)
- **Two-sided:** either pass `abs(Z_vec)` from R, or use `do_abs=TRUE`.
