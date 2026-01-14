# Part 3: Hierarchical Top-Down Variant with Step-Down and Variance Stabilization

## Overview

A hierarchical (top-down) variant that is computationally practical and lets you report parcels, multiscale cubes, and (optionally) voxels, while keeping strong FWER control via an explicit α-budget ("α-spending") recursion.

This hierarchy is an **octree-style dyadic cube decomposition**, which makes the hierarchy nested (important for clean hierarchical error control). If you need sliding windows, the hierarchy/nesting guarantee is gone and you'd revert to global max-calibration or more complex closure.

---

## 1) Core Ingredients

### Set Score Used at Every Node (Parcel/Cube/Voxel)

For a region R with voxel index set idx_R, prior π(v) (sums to 1 over brain), and statistic map Z(v):
- Let `pi_R = pi[idx_R]`, `Z_R = Z[idx_R]`, `den = sum(pi_R)`

For a grid of κ values `kappa_grid = c(0, 0.5, 1, 2, ...)`:

**κ = 0 (weighted mean):**
```
S_0(R) = Σ_{v∈R} π(v)·Z(v) / Σ_{v∈R} π(v)
```

**κ > 0 (soft-max):**
```
S_κ(R) = (1/κ)·log((Σ_{v∈R} π(v)·e^{κ·Z(v)}) / (Σ_{v∈R} π(v)))
```

**Node's omnibus score:**
```
S(R) = max_{κ ∈ kappa_grid} S_κ(R)
```

### Permutation P-Value at a Node

Given permuted maps Z^{(b)}, define S^{(b)}(R) the same way and use:
```
p_R = (1 + #{b: S^{(b)}(R) ≥ S(R)}) / (B+1)
```

### Hierarchical FWER Control via α-Spending

Each node receives a budget `alpha_node`. We:
- Test the node at `alpha_test` (typically `alpha_test = gamma * alpha_node` for internal nodes; `alpha_test = alpha_node` for leaves)
- If significant, split remaining budget among children in proportion to prior mass:
```
α_child = (α_node - α_test) · w_child
w_child = π(child) / Σ π(children)
```

This yields a clean induction proof of strong FWER (union bound + nesting).

---

## 2) R-Style Pseudocode Using neuroim2 Volumes

### Node Structure

```r
new_node <- function(type, id, level, idx, grid, bbox, pi_mass) {
  list(
    type    = type,     # "parcel" | "cube" | "voxel"
    id      = id,
    level   = level,
    idx     = idx,      # integer vector of 1D voxel indices
    grid    = grid,     # matrix n x 3 of (i,j,k) for idx
    bbox    = bbox,     # c(x0,x1,y0,y1,z0,z1) integer voxel coords
    pi_mass = pi_mass,  # sum(pi[idx])
    children = NULL
  )
}
```

### Prior Prep: Normalize, Robust-Mix with Uniform, Renormalize

```r
prep_prior <- function(pi_raw, mask, eta=0.9) {
  # mask: LogicalNeuroVol
  # pi_raw: DenseNeuroVol (same space)
  idx_mask <- which(mask)        # 1D indices inside mask

  pi <- pi_raw
  pi[!mask] <- 0                 # zero out outside mask
  pi[pi < 0] <- 0

  s <- sum(pi[idx_mask])
  if (s <= 0) stop("Prior has zero mass in mask")

  pi <- pi / s                   # normalize to sum 1 in mask

  N <- length(idx_mask)
  pi[idx_mask] <- (1 - eta) * (1/N) + eta * pi[idx_mask]   # robust mix
  pi <- pi / sum(pi[idx_mask])

  return(pi)
}
```

### Score for One Set R, with Omnibus Over kappa_grid

```r
score_set <- function(Zvol, pi, idx_R, kappa_grid) {
  Zr  <- Zvol[idx_R]
  pir <- pi[idx_R]

  den <- sum(pir)
  if (den <= 0) return(list(score=-Inf, best_kappa=NA_real_))

  best <- -Inf
  best_k <- NA_real_

  for (kappa in kappa_grid) {
    if (kappa == 0) {
      sc <- sum(pir * Zr) / den
    } else {
      # stability trick in real code: subtract max(kappa*Zr) before exp
      num <- sum(pir * exp(kappa * Zr))
      sc  <- (1/kappa) * log(num / den)
    }
    if (sc > best) { best <- sc; best_k <- kappa }
  }

  list(score=best, best_kappa=best_k)
}
```

### Permutation P-Value for a Node

```r
pvalue_node <- function(node, Z_obs, pi, perm_fun, B, kappa_grid) {
  obs <- score_set(Z_obs, pi, node$idx, kappa_grid)
  s_obs <- obs$score

  ge <- 0
  for (b in 1:B) {
    Zb <- perm_fun(b)  # returns DenseNeuroVol
    sb <- score_set(Zb, pi, node$idx, kappa_grid)$score
    if (sb >= s_obs) ge <- ge + 1
  }

  p <- (ge + 1) / (B + 1)
  list(p=p, score=s_obs, best_kappa=obs$best_kappa)
}
```

### Build Parcel Roots

```r
build_parcel_roots <- function(Z_obs, mask, atlas_labels, pi) {
  idx_mask <- which(mask)

  # Restrict to mask
  lab <- atlas_labels[idx_mask]
  labs <- sort(unique(lab))

  roots <- list()
  for (L in labs) {
    idx_L <- idx_mask[lab == L]
    if (length(idx_L) == 0) next

    g <- index_to_grid(Z_obs, idx_L)
    bbox <- c(min(g[,1]), max(g[,1]),
              min(g[,2]), max(g[,2]),
              min(g[,3]), max(g[,3]))
    pi_mass <- sum(pi[idx_L])

    type <- if (L == 0) "parcel_remainder" else "parcel"
    roots[[length(roots)+1]] <- new_node(type, paste0("parcel_", L),
                                         level=1, idx=idx_L, grid=g,
                                         bbox=bbox, pi_mass=pi_mass)
  }

  roots
}
```

### Octree Split: Partition Node into Up to 8 Child Cubes (Dyadic Split)

Children are DISJOINT and UNION = parent (within the node's idx list), which makes nesting and alpha-spending clean.

```r
octree_split <- function(node, pi, min_voxels=30, min_edge=2, min_pi_mass=1e-10) {
  idx <- node$idx
  g   <- node$grid
  if (length(idx) < min_voxels) return(list())

  x0 <- node$bbox[1]; x1 <- node$bbox[2]
  y0 <- node$bbox[3]; y1 <- node$bbox[4]
  z0 <- node$bbox[5]; z1 <- node$bbox[6]

  ex <- x1 - x0 + 1
  ey <- y1 - y0 + 1
  ez <- z1 - z0 + 1
  if (max(ex,ey,ez) <= min_edge) return(list())

  xm <- floor((x0 + x1)/2)
  ym <- floor((y0 + y1)/2)
  zm <- floor((z0 + z1)/2)

  # define 8 sub-boxes: [x0..xm]x[y0..ym]x[z0..zm], etc.
  boxes <- list(
    c(x0,xm, y0,ym, z0,zm),
    c(xm+1,x1, y0,ym, z0,zm),
    c(x0,xm, ym+1,y1, z0,zm),
    c(xm+1,x1, ym+1,y1, z0,zm),
    c(x0,xm, y0,ym, zm+1,z1),
    c(xm+1,x1, y0,ym, zm+1,z1),
    c(x0,xm, ym+1,y1, zm+1,z1),
    c(xm+1,x1, ym+1,y1, zm+1,z1)
  )

  kids <- list()
  for (k in seq_along(boxes)) {
    bb <- boxes[[k]]

    # select indices that fall in the child bbox
    sel <- (g[,1] >= bb[1] & g[,1] <= bb[2] &
            g[,2] >= bb[3] & g[,2] <= bb[4] &
            g[,3] >= bb[5] & g[,3] <= bb[6])

    idx_k <- idx[sel]
    if (length(idx_k) == 0) next

    pi_mass_k <- sum(pi[idx_k])
    if (pi_mass_k < min_pi_mass) next

    gk <- g[sel, , drop=FALSE]
    kid <- new_node(type="cube",
                    id=paste0(node$id, "_c", k),
                    level=node$level + 1,
                    idx=idx_k, grid=gk, bbox=bb, pi_mass=pi_mass_k)
    kids[[length(kids)+1]] <- kid
  }

  kids
}
```

### Recursive Hierarchical Testing with Alpha-Spending

```r
test_tree <- function(node, Z_obs, pi, perm_fun, B, kappa_grid,
                      alpha_node,
                      gamma=0.5,
                      min_alpha=1e-6,
                      min_voxels=30,
                      min_edge=2,
                      min_pi_mass=1e-10,
                      hits_accum) {

  if (alpha_node < min_alpha) return(hits_accum)
  if (node$pi_mass < min_pi_mass) return(hits_accum)

  # Determine if splittable (children exist) WITHOUT building all descendants upfront
  children <- octree_split(node, pi,
                           min_voxels=min_voxels,
                           min_edge=min_edge,
                           min_pi_mass=min_pi_mass)

  is_leaf <- (length(children) == 0)

  alpha_test <- if (is_leaf) alpha_node else gamma * alpha_node

  # Compute p-value for this node
  tst <- pvalue_node(node, Z_obs, pi, perm_fun, B, kappa_grid)

  if (tst$p > alpha_test) {
    return(hits_accum)  # fail to reject => prune subtree
  }

  # Reject node => record it
  hits_accum[[length(hits_accum)+1]] <- list(
    id=node$id, type=node$type, level=node$level,
    bbox=node$bbox,
    n_vox=length(node$idx),
    pi_mass=node$pi_mass,
    p=tst$p, alpha_used=alpha_test,
    score=tst$score, best_kappa=tst$best_kappa
  )

  if (is_leaf) return(hits_accum)

  # Allocate remaining alpha to children, weighted by prior mass
  alpha_left <- alpha_node - alpha_test
  kid_masses <- vapply(children, function(ch) ch$pi_mass, numeric(1))
  w <- kid_masses / sum(kid_masses)

  for (i in seq_along(children)) {
    hits_accum <- test_tree(children[[i]], Z_obs, pi, perm_fun, B, kappa_grid,
                           alpha_node = alpha_left * w[i],
                           gamma=gamma,
                           min_alpha=min_alpha,
                           min_voxels=min_voxels,
                           min_edge=min_edge,
                           min_pi_mass=min_pi_mass,
                           hits_accum=hits_accum)
  }

  hits_accum
}
```

### Main Entry: Parcels at Top, Then Octree Inside Each Parcel

```r
hier_threshold <- function(Z_obs, mask, pi_raw, perm_fun, B,
                           atlas_labels,
                           alpha=0.05,
                           kappa_grid=c(0, 0.5, 1, 2),
                           eta=0.9,
                           gamma=0.5,
                           two_sided=TRUE,
                           min_edge=2,
                           min_voxels=30,
                           min_alpha=1e-6,
                           min_pi_mass=1e-10) {

  # Ensure Z is DenseNeuroVol; convert to |Z| if two-sided
  Z <- Z_obs
  if (two_sided) Z <- abs(Z)

  pi <- prep_prior(pi_raw, mask, eta)

  parcel_roots <- build_parcel_roots(Z, mask, atlas_labels, pi)

  # Allocate alpha among parcels proportional to prior mass
  masses <- vapply(parcel_roots, function(n) n$pi_mass, numeric(1))
  w_parc <- masses / sum(masses)

  hits <- list()

  for (j in seq_along(parcel_roots)) {
    hits <- test_tree(parcel_roots[[j]],
                      Z_obs=Z, pi=pi, perm_fun=perm_fun, B=B,
                      kappa_grid=kappa_grid,
                      alpha_node=alpha * w_parc[j],
                      gamma=gamma,
                      min_alpha=min_alpha,
                      min_voxels=min_voxels,
                      min_edge=min_edge,
                      min_pi_mass=min_pi_mass,
                      hits_accum=hits)
  }

  # Optional: build a voxelwise significance mask from leaf hits
  sig_idx <- integer(0)
  for (h in hits) {
    if (h$type == "voxel") sig_idx <- c(sig_idx, h$idx)
  }
  sig_mask <- as.mask(Z, sig_idx)

  list(hits=hits, sig_mask=sig_mask)
}
```

---

## 3) Flexibility Knobs

- **Parcels:** Get parcel calls directly as top-level nodes (parcel_17, etc.)
- **Multiscale cubes:** Get cubes at many scales because of octree splits
- **Voxel localization:** Set `min_edge = 1` and keep splitting until cubes are 1×1×1 (but only inside significant regions due to pruning)
- **Prior influences two places:**
  1. Inside the statistic S_κ(R) (evidence aggregation)
  2. Inside α allocation (more budget where prior mass is high)
- **Permutation engine (perm_fun):** Where your design enters (sign-flip for one-sample, Freedman-Lane, etc.). `perm_fun(b)` should return a DenseNeuroVol Z-map under the null.

---

## A) Statistic with n_eff Scaling (Variance Stabilization)

Let π(v) ≥ 0 with Σ_{v∈mask} π(v) = 1. For any set/region R:

**Define sums:**
- Π(R) = Σ_{v∈R} π(v)
- Π_2(R) = Σ_{v∈R} π(v)²

**Effective size:**
```
n_eff(R) = Π(R)² / Π_2(R)
```

### κ=0 "Diffuse" Score (Variance-Stabilized)

Start from prior-weighted mean in R: Z̄_π(R) = Σ_{v∈R} w_v·Z(v) with w_v = π(v)/Π(R).

Under an iid null, Var(Z̄_π) = Σ w_v² = Π_2(R)/Π(R)².

The standardized diffuse score is:
```
U_0(R) = (Σ_{v∈R} π(v)·Z(v)) / √(Σ_{v∈R} π(v)²)
       = Z̄_π(R) · √n_eff(R)
```

This is the "multiply by √n_eff" correction—fast to compute (just two sums).

### κ>0 "Focal" Score (Soft-Max)

```
S_κ(R) = (1/κ)·log((Σ_{v∈R} π(v)·e^{κ·Z(v)}) / (Σ_{v∈R} π(v)))
```

### Combined Node Score

```
T(R) = max(U_0(R), max_{κ∈K, κ>0} S_κ(R))
```

(Use |Z| for two-sided.)

---

## B) Westfall-Young Step-Down

For a family of hypotheses (parcels, or ≤8 children cubes of a node), compute:
- Observed scores `t_obs[i] = T(R_i)`
- Permutation score matrix `t_perm[b,i] = T^{(b)}(R_i)`

Then perform WY step-down maxT:

```r
wy_stepdown_maxT <- function(t_obs, t_perm, alpha) {
  # t_obs: numeric length m
  # t_perm: numeric matrix B x m
  m <- length(t_obs); B <- nrow(t_perm)

  ord <- order(t_obs, decreasing=TRUE)
  active <- ord
  rejected <- rep(FALSE, m)
  p_adj <- rep(NA_real_, m)

  p_running <- 0

  for (j in seq_len(m)) {
    i <- ord[j]
    if (!(i %in% active)) next

    # max over remaining active hypotheses for each permutation
    u <- apply(t_perm[, active, drop=FALSE], 1, max)

    # step p-value for this stage
    p_step <- (1 + sum(u >= t_obs[i])) / (B + 1)

    # step-down adjusted p-value must be monotone
    p_running <- max(p_running, p_step)
    p_adj[i] <- p_running

    if (p_adj[i] <= alpha) {
      rejected[i] <- TRUE
      active <- setdiff(active, i)
    } else {
      break
    }
  }

  list(rejected=rejected, p_adj=p_adj, order=ord)
}
```

This is the high-power upgrade over single-step maxT, and because we only apply it to (a) parcels and (b) sibling sets (≤8), it's feasible.

---

## C) Full Hierarchical Algorithm with Parcels + Octree Cubes + Step-Down

### Data Structure for a Node (with pi2_mass)

```r
new_node <- function(type, id, level, idx, grid, bbox, pi_mass, pi2_mass) {
  list(
    type=type, id=id, level=level,
    idx=idx,        # linear voxel indices
    grid=grid,      # N x 3 grid coords (i,j,k)
    bbox=bbox,      # c(x0,x1,y0,y1,z0,z1) inclusive
    pi_mass=pi_mass,
    pi2_mass=pi2_mass
  )
}
```

### Node Scoring with U_0 + Softmax (with Stability)

```r
score_node <- function(Z, pi, node, kappa_grid) {
  idx <- node$idx
  Zr <- Z[idx]
  pir <- pi[idx]

  # 1) U0: variance-stabilized diffuse score
  num0 <- sum(pir * Zr)
  den2 <- node$pi2_mass  # sum(pi^2) in node, precomputed
  U0 <- if (den2 > 0) num0 / sqrt(den2) else -Inf

  # 2) Softmax branch
  best_soft <- -Inf
  best_k <- NA_real_
  den1 <- node$pi_mass

  for (kappa in kappa_grid) {
    if (kappa <= 0) next
    # log-sum-exp stabilization
    a <- kappa * Zr
    amax <- max(a)
    num <- sum(pir * exp(a - amax))
    S <- (1/kappa) * (log(num / den1) + amax)
    if (S > best_soft) { best_soft <- S; best_k <- kappa }
  }

  list(score=max(U0, best_soft), U0=U0, soft=best_soft, best_kappa=best_k)
}
```

### Octree Split (with pi2_mass)

```r
octree_split <- function(node, pi, min_voxels=30, min_edge=2, min_pi_mass=1e-10) {

  idx <- node$idx
  g <- node$grid
  if (length(idx) < min_voxels) return(list())

  x0 <- node$bbox[1]; x1 <- node$bbox[2]
  y0 <- node$bbox[3]; y1 <- node$bbox[4]
  z0 <- node$bbox[5]; z1 <- node$bbox[6]

  ex <- x1 - x0 + 1
  ey <- y1 - y0 + 1
  ez <- z1 - z0 + 1
  if (max(ex,ey,ez) <= min_edge) return(list())

  xm <- floor((x0 + x1)/2)
  ym <- floor((y0 + y1)/2)
  zm <- floor((z0 + z1)/2)

  boxes <- list(
    c(x0,xm, y0,ym, z0,zm), c(xm+1,x1, y0,ym, z0,zm),
    c(x0,xm, ym+1,y1, z0,zm), c(xm+1,x1, ym+1,y1, z0,zm),
    c(x0,xm, y0,ym, zm+1,z1), c(xm+1,x1, y0,ym, zm+1,z1),
    c(x0,xm, ym+1,y1, zm+1,z1), c(xm+1,x1, ym+1,y1, zm+1,z1)
  )

  kids <- list()
  for (k in seq_along(boxes)) {
    bb <- boxes[[k]]
    sel <- (g[,1] >= bb[1] & g[,1] <= bb[2] &
            g[,2] >= bb[3] & g[,2] <= bb[4] &
            g[,3] >= bb[5] & g[,3] <= bb[6])
    idx_k <- idx[sel]
    if (length(idx_k) == 0) next

    pi_mass <- sum(pi[idx_k])
    if (pi_mass < min_pi_mass) next

    pi2_mass <- sum((pi[idx_k])^2)
    gk <- g[sel, , drop=FALSE]

    kids[[length(kids)+1]] <- new_node(
      type="cube",
      id=paste0(node$id, "_c", k),
      level=node$level + 1,
      idx=idx_k, grid=gk, bbox=bb,
      pi_mass=pi_mass, pi2_mass=pi2_mass
    )
  }

  kids
}
```

### Build Parcel Roots (with pi2_mass)

```r
build_parcel_roots <- function(Z, mask, atlas_labels, pi) {
  idx_mask <- which(mask)
  labs <- sort(unique(atlas_labels[idx_mask]))

  roots <- list()
  for (L in labs) {
    idx_L <- idx_mask[ atlas_labels[idx_mask] == L ]
    if (length(idx_L) == 0) next

    g <- index_to_grid(Z, idx_L)

    bbox <- c(min(g[,1]), max(g[,1]),
              min(g[,2]), max(g[,2]),
              min(g[,3]), max(g[,3]))

    roots[[length(roots)+1]] <- new_node(
      type="parcel",
      id=paste0("parcel_", L),
      level=1,
      idx=idx_L, grid=g, bbox=bbox,
      pi_mass=sum(pi[idx_L]),
      pi2_mass=sum((pi[idx_L])^2)
    )
  }
  roots
}
```

---

## D) Hierarchical Scan: Step-Down at Parcel Level + Step-Down Among Children

No giant permutation-by-hypothesis matrix beyond parcels.

### Helper: Score a Whole Family (Observed and Permuted)

```r
score_family_obs <- function(Z_obs, pi, nodes, kappa_grid) {
  vapply(nodes, function(nd) score_node(Z_obs, pi, nd, kappa_grid)$score, numeric(1))
}

score_family_perm <- function(perm_fun, B, pi, nodes, kappa_grid) {
  m <- length(nodes)
  t_perm <- matrix(NA_real_, nrow=B, ncol=m)
  for (b in 1:B) {
    Zb <- perm_fun(b)  # returns DenseNeuroVol Z map under null
    for (i in 1:m) {
      t_perm[b, i] <- score_node(Zb, pi, nodes[[i]], kappa_grid)$score
    }
  }
  t_perm
}
```

### Recursive Descent into Node's Subtree

```r
descend_node <- function(node, alpha_budget,
                         Z_obs, pi, perm_fun, B, kappa_grid,
                         gamma=0.5,
                         min_alpha=1e-6,
                         min_voxels=30, min_edge=2, min_pi_mass=1e-10,
                         hits) {

  if (alpha_budget < min_alpha) return(hits)

  kids <- octree_split(node, pi,
                       min_voxels=min_voxels, min_edge=min_edge,
                       min_pi_mass=min_pi_mass)
  if (length(kids) == 0) return(hits)

  # Spend part of budget to test the sibling family (step-down), reserve the rest for deeper descent
  alpha_test <- gamma * alpha_budget
  alpha_desc <- alpha_budget - alpha_test

  # Step-down across these children
  t_obs  <- score_family_obs(Z_obs, pi, kids, kappa_grid)
  t_perm <- score_family_perm(perm_fun, B, pi, kids, kappa_grid)

  sd <- wy_stepdown_maxT(t_obs, t_perm, alpha=alpha_test)
  rej <- which(sd$rejected)

  if (length(rej) == 0) return(hits)

  # Record rejected children (with WY step-down p_adj *within this sibling family*)
  for (i in rej) {
    hits[[length(hits)+1]] <- list(
      id=kids[[i]]$id, type=kids[[i]]$type, level=kids[[i]]$level,
      bbox=kids[[i]]$bbox,
      n_vox=length(kids[[i]]$idx),
      pi_mass=kids[[i]]$pi_mass,
      p_adj_sibling=sd$p_adj[i],
      alpha_test_sibling=alpha_test,
      score=t_obs[i]
    )
  }

  # Allocate remaining descendant budget among rejected kids by prior mass
  masses <- vapply(kids[rej], function(nd) nd$pi_mass, numeric(1))
  w <- masses / sum(masses)

  for (j in seq_along(rej)) {
    child <- kids[[rej[j]]]
    hits <- descend_node(child, alpha_budget = alpha_desc * w[j],
                         Z_obs=Z_obs, pi=pi, perm_fun=perm_fun, B=B,
                         kappa_grid=kappa_grid,
                         gamma=gamma,
                         min_alpha=min_alpha,
                         min_voxels=min_voxels, min_edge=min_edge, min_pi_mass=min_pi_mass,
                         hits=hits)
  }

  hits
}
```

### Main Entry: Parcel Step-Down + Recurse into Significant Parcels

```r
hier_scan <- function(Z_obs, mask, pi_raw, atlas_labels,
                      perm_fun, B,
                      alpha=0.05,
                      kappa_grid=c(0.5, 1, 2),   # κ=0 handled by U0 separately
                      two_sided=TRUE,
                      eta=0.9,
                      gamma_root=0.5,   # fraction of alpha spent on parcel testing
                      gamma=0.5,        # fraction of each subtree budget spent on each sibling test
                      min_edge=2, min_voxels=30,
                      min_alpha=1e-6, min_pi_mass=1e-10) {

  Z <- Z_obs
  if (two_sided) Z <- abs(Z)

  pi <- prep_prior(pi_raw, mask, eta)

  parcels <- build_parcel_roots(Z, mask, atlas_labels, pi)

  # ---- STEP-DOWN AT PARCEL ROOTS ----
  alpha_parcel_test <- gamma_root * alpha
  alpha_desc_total  <- alpha - alpha_parcel_test

  t_obs_par  <- score_family_obs(Z, pi, parcels, kappa_grid)
  t_perm_par <- score_family_perm(perm_fun, B, pi, parcels, kappa_grid)

  sd_par <- wy_stepdown_maxT(t_obs_par, t_perm_par, alpha=alpha_parcel_test)
  rej_par <- which(sd_par$rejected)

  hits <- list()

  # record rejected parcels
  for (i in rej_par) {
    hits[[length(hits)+1]] <- list(
      id=parcels[[i]]$id, type="parcel", level=parcels[[i]]$level,
      bbox=parcels[[i]]$bbox,
      n_vox=length(parcels[[i]]$idx),
      pi_mass=parcels[[i]]$pi_mass,
      p_adj_sibling=sd_par$p_adj[i],
      alpha_test_sibling=alpha_parcel_test,
      score=t_obs_par[i]
    )
  }

  if (length(rej_par) == 0) return(list(hits=hits))

  # Allocate descendant budget among rejected parcels by prior mass
  masses <- vapply(parcels[rej_par], function(nd) nd$pi_mass, numeric(1))
  w <- masses / sum(masses)

  # ---- RECURSE WITH STEP-DOWN AMONG OCTREE CHILDREN ----
  for (j in seq_along(rej_par)) {
    pnode <- parcels[[rej_par[j]]]
    hits <- descend_node(pnode,
                         alpha_budget = alpha_desc_total * w[j],
                         Z_obs=Z, pi=pi, perm_fun=perm_fun, B=B,
                         kappa_grid=kappa_grid,
                         gamma=gamma,
                         min_alpha=min_alpha,
                         min_voxels=min_voxels, min_edge=min_edge, min_pi_mass=min_pi_mass,
                         hits=hits)
  }

  list(hits=hits)
}
```

---

## Notes for Package Docs

- **Why "more sensitive than GRF peak-height":** GRF peak-height is fundamentally a max-statistic approach. The κ=0 standardized mean branch U_0 plus κ>0 softmax branch explicitly targets distributed + multifocal patterns that max-statistic procedures are weak for, while WY controls false positives via the joint permutation distribution.

- **Why U_0 fixes voxel-dominance:** Without scaling, voxel sets (high variance) dominate the max pool; U_0 is roughly N(0,1) for any set size under simple null models, so parcels/cubes get a fair shot.

- **Why step-down is feasible here:** You only step-down on (i) parcels (hundreds) and (ii) sibling groups (≤8), so you avoid the impossible global "store all scores for all cubes."

---

## Computational Profile: Where the Time Goes

### Notation
- N: # voxels in the analysis mask (often 100k-300k)
- P: # parcels (often 100-400)
- B: # permutations/sign-flips (often 1k-10k)
- K: # softmax κ values (>0) in the grid (often 2-6)
- Tree: visited nodes n_1,...,n_M with voxel counts |n_m|. Max depth D (typically 5-8)

### Stage-by-Stage Timing

**0) One-time preprocessing (usually not bottleneck)**
- Normalize π, compute π², build mask index map: O(N)
- Build parcel roots + voxel index lists: O(N)
- Compute and cache per-node constants: O(N) total across created nodes

**1) Parcel-root step-down testing (dominant in most runs)**

A good vectorized implementation costs:
- κ=0 branch: one pass through N voxels for parcel numerators (group-sum): O(N)
- for each κ>0: compute exp(κZ) (O(N)) plus another group-sum (O(N)) → O(KN)

Parcel permutation cost: **O(B · (K+1) · N)**

A bad loop implementation has the same big-O but much worse constant due to:
- Repeated `Z[idx]` vector allocations
- Repeated `pi[idx]` allocations
- R interpreter overhead for each parcel

The WY step-down itself is cheap: ~O(BP) for repeated max computations.

**2) Hierarchical descent (often smaller than parcel stage)**

At each visited node:
- Split into ≤8 children
- Compute child scores for all permutations
- Step-down among ≤8 kids

If you score each child separately by extracting `Z[idx_child]`, you pay repeated allocations. Total voxel count across children equals the parent, so arithmetic cost per level is:
```
O(B · (K+1) · |parent voxels|)
```

In practice, pruning makes this much smaller—only significant parcels spawn subtrees, most sibling families won't be significant.

**3) Memory profile**
- Parcel permutation matrix t_perm: size B × P (e.g., 5000 × 200 ≈ 8 MB)
- Sibling permutation matrices: B × ≤8 per family, tiny
- Node index lists: manageable if you prune

---

## Two Most Impactful Accelerations

### Acceleration 1: Vectorize Parcel Scoring

Goal: compute all parcel scores for a permutation in ~one pass over voxels per κ.

Work in compressed vector space of masked voxels:

```r
# Precompute once:
lab <- atlas_labels[mask_idx]              # integer labels length N
pi_vec  <- pi[mask_idx]
pi2_vec <- pi_vec * pi_vec

den1 <- rowsum(pi_vec,  group=lab, reorder=FALSE)   # P x 1
den2 <- rowsum(pi2_vec, group=lab, reorder=FALSE)   # P x 1
inv_sqrt_den2 <- 1 / sqrt(den2)
log_den1 <- log(den1)

# For each permutation b:
Zb_vec <- Zb[mask_idx]    # already two-sided if using abs()

# (1) κ=0 standardized diffuse statistic U0 for all parcels in one shot:
num0 <- rowsum(pi_vec * Zb_vec, group=lab, reorder=FALSE)
U0   <- num0 * inv_sqrt_den2

# (2) κ>0 softmax statistics for all parcels:
best_soft <- rep(-Inf, nrow(den1))

for (kappa in kappa_grid[kappa_grid > 0]) {
  a <- kappa * Zb_vec
  amax <- max(a)                     # global stabilization
  numk <- rowsum(pi_vec * exp(a - amax), group=lab, reorder=FALSE)
  Sk   <- (amax + log(numk) - log_den1) / kappa
  best_soft <- pmax(best_soft, Sk)
}

T_parcel <- pmax(U0, best_soft)      # omnibus score per parcel
```

This removes inner "for each parcel" loops—typically **10×-100× improvement**.

**Optional:** Use sparse membership matrix M (P × N indicator), then:
```r
num0 = M %*% (pi_vec * Zb_vec)
numk = M %*% (pi_vec * exp(κ * Zb_vec))
```
Extremely fast as single sparse MV multiply per κ.

### Acceleration 2: Cache π, π², and Node Index/Weight Lookups

Represent everything in masked-voxel vector index space 1..N:

```r
make_node_cached <- function(idx_maskspace, pi_vec) {
  pir <- pi_vec[idx_maskspace]
  den1 <- sum(pir)
  den2 <- sum(pir * pir)
  list(
    idx = idx_maskspace,
    pir = pir,                      # cache weights
    log_den1 = log(den1),
    inv_sqrt_den2 = 1 / sqrt(den2)
  )
}

score_node_cached <- function(Z_vec, node, kappa_grid) {
  Zr <- Z_vec[node$idx]             # still allocates in R
  pir <- node$pir

  # κ=0 standardized mean
  U0 <- sum(pir * Zr) * node$inv_sqrt_den2

  best_soft <- -Inf
  for (kappa in kappa_grid[kappa_grid > 0]) {
    a <- kappa * Zr
    amax <- max(a)
    soft <- (log(sum(pir * exp(a - amax))) + amax - node$log_den1) / kappa
    best_soft <- max(best_soft, soft)
  }

  max(U0, best_soft)
}
```

**Production-grade move:** Implement inner scoring loop in Rcpp so it iterates over idx and reads Z_vec by pointer without allocating Zr.

```cpp
double score_node_cpp(const NumericVector& Z_vec,
                      const IntegerVector& idx,
                      const NumericVector& pir,
                      double log_den1,
                      double inv_sqrt_den2,
                      const NumericVector& kappa_grid);
```

This single change is often another **5×-20× improvement** in the recursive stage.

---

## Summary

RFT peak inference is fast because it avoids resampling by using expected EC to approximate max tail in smooth fields. This method pays computational price to (i) be flexible across parcels/cubes and κ-grid and (ii) be empirically calibrated. The right engineering target: make "score many sets under many permutations" cheap, which is exactly what these two accelerations do.
