# Algorithm Specifications

## Core Algorithm: LR-MFT

**Full details:** `docs/discussion/part-01-grf-vs-lrmft.md`

### Prior-Weighted Soft Evidence Statistic

```
T_κ(R) = log Σ_{v∈R} π(v) · exp(κ · Z(v))
```

Where:
- R = region (set of voxel indices)
- π(v) = prior weight at voxel v (default: uniform)
- Z(v) = Z-score at voxel v
- κ = temperature parameter (default: 1.0)

### Variance-Stabilized Score

```
U_0(R) = Σ_{v∈R} π(v)·Z(v) / sqrt(Σ_{v∈R} π(v)²)
n_eff = (Σπ)² / Σπ²
```

Used for fair comparison across regions of different sizes.

---

## Hierarchical Inference

**Full details:** `docs/discussion/part-03-hierarchical-stepdown.md`

### α-Spending (Top-Down)
1. Test parent region at level α
2. If rejected, test children at α/n_children each
3. Recurse until leaf nodes or non-rejection

### Westfall-Young Step-Down
1. Sort observed scores descending
2. For each permutation, compute successive maxima
3. Adjusted p-value = proportion of perms where max ≥ observed
4. Enforce monotonicity: p_adj[i] = max(p_adj[i], p_adj[i-1])

---

## Octree Partitioning

**Full details:** `docs/discussion/part-04-onepass-rcpp.md`

### Split Logic
1. Compute bounding box of indices
2. Find midpoint in each dimension
3. Assign each voxel to one of 8 octants
4. Recurse until region < min_voxels

### One-Pass Scoring
Score all 8 children in single O(n) pass:
1. Compute octant assignment for each voxel
2. Accumulate π·exp(κ·Z) into 8 accumulators
3. Return log of each accumulator

---

## Baseline Methods

**Full details:** `docs/discussion/part-05-bbox-helpers-and-baselines.md`

### RFT (Random Field Theory)
- Peak-level: expected Euler characteristic approximation
- Cluster-level: cluster extent distribution from GRF theory

### TFCE (Threshold-Free Cluster Enhancement)
```
TFCE(v) = ∫₀^{Z(v)} e(h)^E · h^H dh
```
Default: E=0.5, H=2.0, dh=0.1

### Cluster-FDR
- Form clusters at primary threshold
- Compute cluster p-values (RFT or permutation)
- Apply Benjamini-Hochberg to cluster p-values

---

## Canonicalization

Convert input stat maps to Z-scores:

| Input Type | Conversion |
|------------|------------|
| Z | Identity |
| t(df) | `qnorm(pt(t, df))` |
| -log10(p) | `qnorm(1 - 10^(-x))` |
