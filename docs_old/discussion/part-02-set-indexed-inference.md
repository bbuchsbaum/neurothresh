# Part 2: Set-Indexed Inference with Prior-Weighted Evidence Aggregation

## Problem Statement

Given any probability map π(v) ≥ 0 over voxels/resels with Σ_{v∈V} π(v) = 1, we want an inference scheme that:
1. Uses π in a principled way (not just as a display overlay)
2. Lets you legitimately declare either a voxel significant or a region/set R ⊂ V significant (localization)

Key framing:
- **RFT peak inference** = hard-max on a smooth field (calibrated by EC)
- **What we need** = a family of tests indexed by regions/sets, with prior-weighted evidence aggregation, and simultaneous error control

---

## 1) RFT Baseline: Max (or Extent) Detector on a Smooth Field

Worsley (2003) describes two classic RFT approaches: (i) maximum statistic and (ii) cluster extent. The max approach chooses height z so Pr(max Z > z) ≈ 0.05 via EC/resels.

The "GRF peak-height baseline" is:
```
T_max = max_{v∈V} Z(v)
reject if T_max > u_α
```
with u_α calibrated by EC/RFT (or by permutation in modern practice).

This gives voxel localization by thresholding Z(v), but it's a hard-max rule.

---

## 2) Set-Indexed Hypotheses

Let each voxel/resel v ∈ V have a null hypothesis:
```
H_v: "no activation at v"
```

For any region/set R ⊆ V, define the set null:
```
H_R: ∩_{v∈R} H_v  ("no activation anywhere in R")
```

Requirements:
- "Call any voxel significant" = test H_{v}
- "Call any region/set R significant" = test H_R

We need:
1. A test statistic T(R) for each R that can incorporate π
2. Simultaneous error control across all the R's you might report

---

## 3) Prior-Weighted Statistic for Voxels AND Regions

Assume you have a smooth/whitened statistic field Z(v) (t→Z, or Z directly).

Define, for a tuning parameter κ > 0 (think "target effect size / temperature"):
```
T_κ(R) := log(Σ_{v∈R} π(v)·e^{κ·Z(v)})
```

**Key properties:**

**Voxel case R = {v}:**
```
T_κ({v}) = log π(v) + κ·Z(v)
```
A voxelwise score with an explicit prior term (large π(v) gives a boost).

**Region case:** Aggregates evidence across R (soft-max / log-sum-exp).

**As κ → ∞:**
```
(1/κ)·T_κ(R) → max_{v∈R} Z(v)
```
So RFT peak inference is a limiting special case (uniform π, κ large, region = whole brain).

This directly addresses the "any probability map π" requirement: π only needs to be nonnegative and sum to 1.

---

## 4) Simultaneous FWER Control via Permutation

### Define the Set Family F

Pick a family of sets F you want the right to report:
- All singletons {v} (voxels)
- Atlas ROIs (parcels, networks)
- Prior-defined sets like "top 10% prior mass", "top 1% prior mass"
- Multiscale tiles/windows (cubes/supervoxels at multiple sizes)

### Max Over Reported Objects Statistic

```
M_κ := max_{R∈F} T_κ(R)
```

### Permutation (or Sign-Flip) Calibration

Using standard neuroimaging permutation/sign-flip machinery (across subjects), generate null replicates of the whole field Z^{(b)}(·) and compute:
```
M_κ^{(b)} = max_{R∈F} T_κ^{(b)}(R)
```

Set the threshold c_α to the empirical (1-α) quantile of {M_κ^{(b)}}.

### Decision Rule (Simultaneous)

Declare any region/voxel R ∈ F significant if:
```
T_κ(R) ≥ c_α
```

### Theorem (Strong FWER Control over F)

Under the global null (all H_v true) and valid permutations/sign-flips:
```
Pr_0(∃ R∈F: we reject H_R) = Pr_0(M_κ ≥ c_α) ≤ α
```

**Proof sketch:** c_α is chosen from the permutation distribution of M_κ, so by construction the probability (under the null) that the observed M_κ exceeds c_α is α. Any false rejection implies M_κ ≥ c_α. ∎

This gives localization with honest familywise control across:
- Voxels (singletons)
- ROIs
- "Sets" induced by the prior map

No EC approximation is needed.

---

## 5) Why This Is More Sensitive Than GRF Peak-Height

RFT peak-height is essentially using max Z(v). That detector is optimal when the alternative is "one dominant peak."

But many fMRI alternatives are **multi-focal or distributed** (networks, bilateral activation, etc.). For those, a statistic that adds evidence is better.

### 5.1 Clean Alternative Where T_κ(R) Is Point-Optimal LR

Work in a region R. Suppose under the alternative there are K ≥ 2 activated foci inside R, with locations i.i.d. from normalized prior π_R(v) = π(v)/π(R). Assume matched-filtered field values at those foci are approximately independent.

For fixed "target" effect size μ > 0, the likelihood ratio for H_R vs that multi-focal alternative has the form:
```
Λ(y;R) ∝ (Σ_{v∈R} π_R(v)·e^{μ·Z(v)})^K
```

So rejecting for large Σ_{v∈R} π_R(v)·e^{μ·Z(v)} (equivalently large T_μ(R)) is the **Neyman-Pearson most powerful test** for that (μ,K).

RFT max is not equivalent unless you take μ→∞ (log-sum-exp → max) or K=1.

**For multi-focal alternatives, the region statistic is provably more powerful than the max statistic at the same size.**

### 5.2 Transparent Power Separation: Distributed Weak Signal

Assume (for intuition) independent Z(v) ~ N(0,1) under H_0. Let R contain m voxels, and under the alternative every voxel in R has a small mean shift δ > 0:
```
Z(v) ~ N(δ, 1) for v ∈ R
```

**Hard-max (RFT peak):** Requires at least one voxel exceed a high threshold u_α (on order of √(2 log N)). If δ ≪ u_α, power is poor.

**Aggregate test:** Consider sum statistic S_R = Σ_{v∈R} Z(v). Then:
```
S_R ~ N(m·δ, m)
```
The standardized sum S_R/√m ~ N(δ·√m, 1) gains a factor √m in mean. Even when each voxel is weak, the region-level aggregate can be strongly significant if m is large enough.

T_κ(R) behaves like a soft aggregate; for small/moderate κ it is closer to averaging/summing than to max.

**The regime of superiority is distributed/multi-focal, which is biologically plausible.**

---

## 6) Handling "Any Set R" vs "Any Set in a Family F"

**Trade-off:**
- If you want to claim significance for **any set R** you dream up after looking at the data (all 2^N subsets), you need closed testing / all-resolution inference machinery (very strong but complex).
- If you want to claim significance for **any set you predefine in advance** (voxels + atlas ROIs + multiscale tiling + prior-defined regions), then the permutation max over F approach is rigorous and straightforward.

In practice, defining F as:
- All voxels
- Atlas ROIs
- Multiscale tiling (cubes/supervoxels at multiple sizes)

already gives you "any reasonable area/set."

---

## 7) Localization Products (Outputs)

Once you have adjusted p-values / thresholds for all R ∈ F:

1. **Voxel map:** All v such that T_κ({v}) ≥ c_α
2. **ROI table:** All ROIs R such that T_κ(R) ≥ c_α
3. **Prior-aware "most likely locus"** inside a significant region (descriptive):
   ```
   w(v|R) ∝ π(v)·e^{κ·Z(v)}, v ∈ R
   ```
   Principled way to rank voxels inside a significant region, consistent with the prior.

---

## 8) Choosing κ Without Cheating

κ encodes the target effect size. Three honest options:

1. **Design choice:** Choose κ to target the smallest effect you care about
2. **Bayes mixture:** Integrate over κ with a prior g(κ):
   ```
   T(R) = log ∫ Σ_{v∈R} π(v)·e^{κ·Z(v)}·g(κ) dκ
   ```
   (a single statistic again)
3. **Grid search over κ** and include it in the max calibration (like "searching scale space")

Option (3) is the most common pragmatic solution.

---

## Pseudocode: Prior-Weighted Multiscale Set Scan with Permutation Max (FWER)

### Key Statistic

Let Z(v) be the (smoothed/standardized) statistic field inside brain mask. Let π(v) be the prior map on voxels (masked + normalized).

Define conditional prior mass in a set:
```
π(R) = Σ_{v∈R} π(v)
π_R(v) = π(v) / π(R)
```

Define soft-max family parameterized by κ ≥ 0:

**For κ > 0:**
```
S_κ(R) = (1/κ)·log(Σ_{v∈R} π_R(v)·e^{κ·Z(v)})
       = (1/κ)·log((Σ_{v∈R} π(v)·e^{κ·Z(v)}) / (Σ_{v∈R} π(v)))
```

**For κ = 0 (the limit):**
```
S_0(R) = Σ_{v∈R} π_R(v)·Z(v)
       = (Σ_{v∈R} π(v)·Z(v)) / (Σ_{v∈R} π(v))
```

**Important properties for localization:**
- For singleton voxel R = {v}: S_κ({v}) = Z(v) for all κ. Voxels remain "callable" in usual Z-units (no log π(v) penalty).
- Small κ behaves like prior-weighted average; large κ behaves like max over set.
- Scanning over small grid of κ gives flexible detector for distributed vs focal patterns.

### Algorithm

```
ALGORITHM PriorWeightedMultiscaleScanFWER

INPUTS:
  - data: subject-level contrast images or GLM-ready data
  - design: design matrix X and contrast c (or one-sample sign-flip setup)
  - mask: 3D boolean brain mask (V voxels)
  - pi_raw: nonnegative prior map over voxels (same grid)
  - atlas_labels: optional 3D integer label map for parcels (0 = none)
  - scales: list of cube edge lengths in voxels, e.g. [4, 6, 8, 12]
  - strides: dict scale -> stride (or default stride = floor(scale/2))
  - kappa_grid: list of κ values, include 0 plus a few >0, e.g. [0, 0.5, 1.0, 2.0]
  - alpha: desired FWER level (e.g. 0.05)
  - B: number of permutations / sign-flips (e.g. 5000+)
  - two_sided: boolean (if true, use |Z|; else use Z as-is)
  - eta: robust prior mixing in [0,1]; recommended eta ~ 0.9
         (eta=1 uses pure prior; eta<1 mixes with uniform to avoid zeros)

OUTPUTS:
  - significant_sets: list of (set_id, set_type, p_adj, S_obs, best_kappa)
  - significant_voxels: list of (voxel_index, p_adj, Z_obs_at_voxel)

PROCEDURE:

  # ----------------------------
  # 0) Index voxels and normalize prior
  # ----------------------------
  V = { v : mask[v] == 1 }
  N = |V|
  pi = pi_raw * mask
  if sum(pi) == 0: ERROR("Prior map has zero mass inside mask")

  pi = pi / sum(pi)                             # normalize to sum 1 over brain
  pi = (1-eta) * (1/N) + eta * pi               # robust mixture (prevents pi=0)
  pi = pi / sum(pi)                             # renormalize (numerical)

  # Precompute a prefix sum of pi for fast cube pi(R)
  PI_PREFIX = IntegralImage3D(pi)               # 3D summed-area table

  # ----------------------------
  # 1) Build candidate set family F (flexible: parcels + multiscale cubes + voxels)
  # ----------------------------
  F = empty list of sets
  # Each "set" record stores:
  #   - id
  #   - type in {"parcel","cube","voxel"}
  #   - representation:
  #       parcel: list_of_voxel_indices
  #       cube: bounding box (x0,x1,y0,y1,z0,z1) half-open
  #       voxel: single voxel index
  #   - pi_mass = pi(R)  (computed once)

  # 1a) Add parcels
  if atlas_labels provided:
    for each parcel label L > 0:
      idx = [v in V where atlas_labels[v] == L]
      if |idx| > 0:
        set = NewSet(type="parcel", rep=idx)
        set.pi_mass = sum(pi[v] for v in idx)
        if set.pi_mass > 0: append(F, set)

  # 1b) Add multiscale cubes
  # For cubes, use integral image to compute pi_mass quickly
  for each scale s in scales:
    step = strides.get(s, floor(s/2))
    for x0 in range(0, X_dim - s + 1, step):
      for y0 in range(0, Y_dim - s + 1, step):
        for z0 in range(0, Z_dim - s + 1, step):
          box = (x0, x0+s, y0, y0+s, z0, z0+s)
          # pi mass includes only masked voxels since pi outside mask = 0
          pi_mass = BoxSum3D(PI_PREFIX, box)
          if pi_mass >= MIN_PI_MASS:             # e.g. 1e-6 or tuned
            set = NewSet(type="cube", rep=box)
            set.pi_mass = pi_mass
            append(F, set)

  # 1c) Optionally add all voxels as singleton sets
  # NOTE: This can be very large; consider hierarchical testing (see below).
  for each voxel v in V:
    set = NewSet(type="voxel", rep=v)
    set.pi_mass = pi[v]                          # >0 due to robust mixing
    append(F, set)

  # ----------------------------
  # 2) Compute observed statistic field Z_obs
  # ----------------------------
  Z_obs = ComputeGroupZMap(data, design, mask)     # your GLM -> voxelwise Z
  if two_sided: Z_obs = abs(Z_obs)

  # ----------------------------
  # 3) Compute observed set scores for each κ, and the omnibus score per set
  # ----------------------------
  # We'll store:
  #   S_obs[R]  = max over κ in kappa_grid of S_κ(R)
  #   best_kappa[R] = argmax κ
  for each set R in F:
    S_obs[R] = -INF
    best_kappa[R] = NONE

  for each κ in kappa_grid:
    if κ == 0:
      # Need numerator field B0(v) = pi(v) * Z(v)
      B0 = pi * Z_obs
      B0_PREFIX = IntegralImage3D(B0)             # for cube sums
      # parcels and voxels computed by indexing
      for each set R in F:
        num = SumFieldOverSet(B0, B0_PREFIX, R)   # sum_{v in R} pi(v) Z(v)
        den = R.pi_mass                           # sum_{v in R} pi(v)
        score = num / den                         # S_0(R)
        if score > S_obs[R]:
          S_obs[R] = score
          best_kappa[R] = 0
    else:
      # Numerator field Aκ(v) = pi(v) * exp(κ Z(v))
      A = pi * exp(κ * Z_obs)
      A_PREFIX = IntegralImage3D(A)               # for cube sums
      for each set R in F:
        num = SumFieldOverSet(A, A_PREFIX, R)     # sum_{v in R} pi(v) exp(κ Z(v))
        den = R.pi_mass                           # sum_{v in R} pi(v)
        score = (1/κ) * log(num / den)            # S_κ(R)
        if score > S_obs[R]:
          S_obs[R] = score
          best_kappa[R] = κ

  # Observed global max across all candidate sets:
  M_obs = max_{R in F} S_obs[R]

  # ----------------------------
  # 4) Permutation calibration: build null distribution of max statistic
  # ----------------------------
  M_perm = array length B

  for b in 1..B:
    Z_b = ComputePermutedGroupZMap(data, design, mask, permutation_id=b)
    if two_sided: Z_b = abs(Z_b)

    # Compute the max over sets for this permutation (no need to store all set scores)
    max_score_b = -INF

    for each κ in kappa_grid:
      if κ == 0:
        B0b = pi * Z_b
        B0b_PREFIX = IntegralImage3D(B0b)
        for each set R in F:
          num = SumFieldOverSet(B0b, B0b_PREFIX, R)
          den = R.pi_mass
          score = num / den
          if score > max_score_b: max_score_b = score
      else:
        Ab = pi * exp(κ * Z_b)
        Ab_PREFIX = IntegralImage3D(Ab)
        for each set R in F:
          num = SumFieldOverSet(Ab, Ab_PREFIX, R)
          den = R.pi_mass
          score = (1/κ) * log(num / den)
          if score > max_score_b: max_score_b = score

    M_perm[b] = max_score_b

  # ----------------------------
  # 5) FWER-adjusted p-values for every set (voxels + parcels + cubes)
  # ----------------------------
  # Westfall-Young single-step maxT adjusted p-value:
  for each set R in F:
    p_adj[R] = (1 + count_b( M_perm[b] >= S_obs[R] )) / (B + 1)

  # Determine significance:
  significant_sets = []
  significant_voxels = []
  for each set R in F:
    if p_adj[R] <= alpha:
      if R.type == "voxel":
        append(significant_voxels, (R.rep, p_adj[R], Z_obs[R.rep]))
      else:
        append(significant_sets, (R.id, R.type, p_adj[R], S_obs[R], best_kappa[R]))

  RETURN significant_sets, significant_voxels
END
```

### Helper Routines

#### Summed-Area Table (Integral Image) for Fast Cube Sums

```
FUNCTION IntegralImage3D(A):
  # A is a 3D array; outside-mask entries are already 0
  P = zeros_like(A)
  for x in 0..X-1:
    for y in 0..Y-1:
      for z in 0..Z-1:
        P[x,y,z] = A[x,y,z]
                   + P[x-1,y,z] + P[x,y-1,z] + P[x,y,z-1]
                   - P[x-1,y-1,z] - P[x-1,y,z-1] - P[x,y-1,z-1]
                   + P[x-1,y-1,z-1]
        # Treat any P[...] with negative index as 0
  return P

FUNCTION BoxSum3D(P, box):
  # box = (x0,x1,y0,y1,z0,z1) half-open
  # Return sum_{x0<=x<x1, y0<=y<y1, z0<=z<z1} A[x,y,z]
  # Using inclusion-exclusion on prefix sums P.
  return P[x1-1,y1-1,z1-1]
         - P[x0-1,y1-1,z1-1] - P[x1-1,y0-1,z1-1] - P[x1-1,y1-1,z0-1]
         + P[x0-1,y0-1,z1-1] + P[x0-1,y1-1,z0-1] + P[x1-1,y0-1,z0-1]
         - P[x0-1,y0-1,z0-1]
  # Again treat P[negative index] as 0.
```

#### Summing a Field Over a Set (Parcel / Cube / Voxel)

```
FUNCTION SumFieldOverSet(fieldA, prefixA, set R):
  if R.type == "voxel":
    v = R.rep
    return fieldA[v]
  if R.type == "parcel":
    idx = R.rep                  # list of voxels
    return sum(fieldA[v] for v in idx)
  if R.type == "cube":
    box = R.rep
    return BoxSum3D(prefixA, box)
```

---

## Practical Notes

1. **Why permutation instead of EC/RFT calibration?**
   RFT uses EC/resels approximations for Pr(max Z > u) in smooth fields. Here we're taking maxima over a richer family (parcels + cubes + voxels + κ-grid), so permutation is the cleanest way to keep exact FWER.

2. **Computational load**
   If you include all voxels + all cubes at many scales, |F| can be huge. Common practical accelerations:
   - Reduce cube stride (e.g., stride = scale) for fewer windows
   - Drop voxel singletons initially and use hierarchical refinement:
     - Stage 1: parcels + big cubes (FWER α)
     - Stage 2: within significant sets, test voxels/small cubes (gatekeeper or closed testing)

3. **Two-sided**
   Using abs(Z) is simplest; for directional inference, run separately on Z and -Z and Bonferroni-combine (or include both in κ-grid with sign).

4. **No "double dipping"**
   The prior map π and the set family F should be fixed independently of the noise realization for the contrast you're testing (external priors, anatomical maps, independent localizer, etc.).

---

## Summary

You can get exactly what you asked for—any prior probability map π, plus the ability to declare either voxels or regions/sets significant—by:

1. Using a prior-weighted soft evidence statistic:
   ```
   T_κ(R) = log Σ_{v∈R} π(v)·e^{κ·Z(v)}
   ```
   which includes voxelwise inference as a special case

2. Calibrating simultaneously over voxels and your set family F using a permutation/max procedure to control FWER

This approach is rigorously more sensitive than GRF/RFT peak-height inference for biologically plausible multi-focal/distributed alternatives because the region statistic is (point-)NP-optimal for that alternative, while the RFT peak method is a hard-max detector.
