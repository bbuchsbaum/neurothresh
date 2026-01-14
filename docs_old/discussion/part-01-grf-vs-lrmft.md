# Part 1: GRF/RFT vs Likelihood-Ratio Matched-Filter Thresholding

## What GRF/RFT Actually Does

In the Gaussian random field / random field theory (RFT) framework (SPM-style inference), you model the statistic image (Z-, T-, or F-field over the brain) as a smooth Gaussian random field under the null. Then set a threshold so that the probability of any false positive anywhere in the search region is approximately α.

Key point: for high thresholds, the tail probability of the global maximum is approximated by the expected Euler characteristic (EC) of the excursion set above that threshold.

Worsley (2003) writes the max-tail approximation as:

```
Pr(max Z > z) ≈ Σ_{d=0}^D Resels_d · EC_d(z)
```

with Resels and EC densities depending on dimension, smoothness (FWHM), and search region geometry.

The canonical GRF/RFT peak-height approach:
1. Compute a statistic field Z(x)
2. Pick z_α so that Pr_0(sup_x Z(x) > z_α) ≈ α using EC/RFT
3. Declare voxels/peaks significant if Z(x) > z_α

This is a **"max detector"**: tuned for spiky/peak-like alternatives.

---

## A Different Thresholding Approach: LR-MFT

### Biological Assumptions

1. **Null field is approximately Gaussian** (after GLM + whitening) - the working assumption behind RFT peak and extent inference
2. **True activations are spatially extended "blobs"** at the resolution used for inference - BOLD effects are spatially smooth (hemodynamics + vasculature + point spread + explicit smoothing)
3. **Unknown activation location** (whole-brain search) - treating each eligible location as a priori equally likely

These are precisely the conditions where a "peak detector" is not information-theoretically optimal.

---

## The Proposed Method: Likelihood-Ratio Matched-Filter Thresholding (LR-MFT)

### Model Setup

Discretize the search region into N locations (voxels or resels). Let y ∈ ℝ^N be the (standardized) statistic image.

**Null (noise only):**
```
H_0: y ~ N(0, Σ)
```
where Σ captures the spatial covariance (smoothness).

**Signal model (one activation blob somewhere):**

Choose a biologically plausible spatial template s ∈ ℝ^N (e.g., isotropic Gaussian blob with FWHM comparable to expected spatial point-spread). Let s_g denote the template placed at location g.

Assume a single blob at unknown location g, with amplitude μ > 0:
```
H_{1,g}(μ): y ~ N(μ·s_g, Σ)
```

To encode "unknown location" without cheating, use the uniform-location mixture:
```
H_1(μ): y ~ (1/N) Σ_{g=1}^N N(μ·s_g, Σ)
```

### Derivation: Exact Neyman-Pearson Optimal Test

Because H_1(μ) is a single (mixture) distribution for each fixed μ, we can apply the Neyman-Pearson lemma directly. The UMP level-α test for H_0 vs H_1(μ) is the likelihood-ratio test (LRT).

**Step 1: Compute the likelihood ratio**

The Gaussian density ratio for a fixed mean shift μ·s_g is:
```
p(y|H_{1,g}(μ)) / p(y|H_0) = exp(μ·s_g^T·Σ^{-1}·y - (1/2)·μ²·s_g^T·Σ^{-1}·s_g)
```

For the mixture alternative H_1(μ):
```
Λ(y) = p(y|H_1(μ)) / p(y|H_0) = (1/N) Σ_{g=1}^N exp(μ·s_g^T·Σ^{-1}·y - (1/2)·μ²·s_g^T·Σ^{-1}·s_g)
```

If the field is approximately stationary, then s_g^T·Σ^{-1}·s_g is approximately constant in g. Call it A = ||s||²_{Σ^{-1}}. Then:
```
Λ(y) ∝ Σ_{g=1}^N exp(μ·s_g^T·Σ^{-1}·y)
```

Define the whitened matched-filter score:
```
M_g(y) := (s_g^T·Σ^{-1}·y) / √A
```

so that under H_0, each M_g is approximately standard normal. Then:
```
Λ(y) ∝ Σ_{g=1}^N exp(κ·M_g(y)),  where κ := μ·√A
```

**Step 2: The optimal test statistic is a "soft-max"**

Taking logs, the LRT is equivalent to thresholding:
```
T_LR(y) = log Σ_{g=1}^N exp(κ·M_g(y))
```

- When κ is large, T_LR is essentially κ·max_g M_g (because log-sum-exp → max)
- When κ is moderate, it aggregates evidence more efficiently than a hard max

**Step 3: Set the threshold to control familywise error**

Pick c_α such that under H_0:
```
Pr_0(T_LR(y) ≥ c_α) = α
```

This can be obtained by:
- Permutation/sign-flipping over subjects (exact under exchangeability)
- Parametric bootstrap using the fitted Σ
- Or a GRF approximation

---

## Proof of Superiority vs GRF Max-Threshold Approach

### What We Compare Against

The canonical GRF peak-height/RFT approach is a max-based whole-brain detector:
```
T_GRF(y) = max_{i=1,...,N} y_i
```
with threshold chosen so that Pr_0(T_GRF ≥ t_α) = α via the EC heuristic.

### Theorem (Uniform Dominance for the Blob-Alternative Mixture)

Fix any μ > 0. Consider testing H_0: y ~ N(0,Σ) versus the whole-brain "one blob somewhere" alternative H_1(μ): y ~ (1/N)Σ_{g=1}^N N(μ·s_g, Σ), where s_g is a spatially extended activation template and the mixture over g is uniform.

Let φ_LR be the level-α LRT based on T_LR, and let φ_GRF be any other level-α test (in particular, the GRF/RFT max test).

Then, for every μ > 0:
```
Pr_{H_1(μ)}(φ_LR = 1) ≥ Pr_{H_1(μ)}(φ_GRF = 1)
```
with strict inequality unless φ_GRF is almost surely a monotone function of Λ(y).

**Proof:** For each fixed μ, H_1(μ) is a single (mixture) distribution and H_0 is a single distribution. By the Neyman-Pearson lemma, among all tests with size α, the likelihood ratio test achieves maximum power against H_1(μ). Therefore, no other level-α test—including the GRF/RFT max-threshold test—can have higher power. QED.

### Why This Superiority is About fMRI Realism

- The GRF max test corresponds to a "delta-like" alternative: "one voxel gets big"
- LR-MFT corresponds to the biologically plausible alternative: "a spatially extended blob appears somewhere"

RFT is solving the wrong detection problem if you believe activations are blobs.

---

## Practical LR-MFT Workflow

1. Compute the usual voxelwise statistic map y (Z or t converted to Z)
2. Estimate Σ (or directly estimate a spatial whitening operator) from GLM residuals
3. Choose a template s (e.g., Gaussian blob with FWHM ≈ expected activation scale)
4. Compute the matched-filter score map:
   ```
   M_g(y) = (s_g^T·Σ^{-1}·y) / √(s_g^T·Σ^{-1}·s_g)
   ```
   Computationally, this is a convolution/correlation of a whitened map with a kernel.
5. Compute the LR statistic:
   ```
   T_LR = log Σ_g exp(κ·M_g)
   ```
6. Calibrate c_α using subject-level permutation/sign-flipping
7. For localization/thresholding, output either:
   - The posterior-like weights w_g ∝ exp(κ·M_g) and threshold them, or
   - A hard-threshold of M_g using a max-corrected threshold

---

## Simulation Verification

### Setup
- Grid: 64×64 (toy cortex/slice)
- Null: Gaussian-smoothed white noise (stationary smooth GRF) rescaled to unit variance
- Signal: add a constant mean shift μ in a 5×5 patch at a random location
- GRF baseline: threshold by global maximum of the raw field (peak-height test)
- LR-MFT: matched-filter scan with covariance-aware weights for a 5×5 patch

Both thresholds chosen so that under null, Pr(any false positive) ≈ 0.05.

### Empirical Null Control (FWER) - 5000 null simulations
- Peak-max test: empirical α ≈ 0.046
- LR-MFT scan: empirical α ≈ 0.055

### Empirical Power - 3000 simulations per μ

| mean shift μ | GRF peak-max power | LR-MFT power |
|--------------|-------------------|--------------|
| 0.0          | 0.045             | 0.054        |
| 0.5          | 0.042             | 0.481        |
| 1.0          | 0.054             | 0.999        |
| 2.0          | 0.176             | 1.000        |
| 3.0          | 0.632             | 1.000        |

**Interpretation:** At equal FWER, the GRF peak approach "almost never detects" until the blob mean is huge, while LR-MFT detects reliably at much smaller blob-level effects because it integrates spatial evidence optimally.

---

## Extension: Prior-Weighted LRT

### Prior Maps

A prior map is a function π(v) ≥ 0 over voxels with Σ_v π(v) = 1, representing a priori probability that the effect center is at v.

Plausible π(v) sources:
- **Anatomical plausibility:** π(v) high in gray matter, low in CSF/white matter
- **Task prior / meta-analysis prior:** π(v) proportional to meta-analytic activation likelihood map
- **Subject-specific functional localizer prior:** π(v) derived from separate run/contrast
- **Network prior:** π(v) high on known network mask

Robust version:
```
π_η(v) = (1-η)·(1/N) + η·π_prior(v),  0 ≤ η ≤ 1
```

### Template Bank for Shape/Scale

Let s_{v,ℓ} ∈ ℝ^N be a template centered at voxel v with "type" ℓ (scale/FWHM, anisotropy, etc.). Put a prior π(v,ℓ) with Σ_{v,ℓ} π(v,ℓ) = 1.

This covers:
- Bank of Gaussian blobs at multiple FWHM (scale prior)
- Cortical-surface diffusion/heat kernels (graph-based physics prior)
- Oriented blobs along sulci (geometry prior)
- Multi-blob network templates (systems prior)

### Prior-Weighted LRT Statistic

**Alternative with priors:**
```
H_1(μ): y ~ Σ_{v,ℓ} π(v,ℓ)·N(μ·s_{v,ℓ}, Σ)
```

**The LRT statistic:**
```
Λ(y) = Σ_{v,ℓ} π(v,ℓ)·exp(μ·s_{v,ℓ}^T·Σ^{-1}·y - (1/2)·μ²·s_{v,ℓ}^T·Σ^{-1}·s_{v,ℓ})
```

Equivalently:
```
T_PW-LR(y) = log Σ_{v,ℓ} π(v,ℓ)·exp(κ_ℓ·M_{v,ℓ}(y))
```

where:
```
M_{v,ℓ}(y) = (s_{v,ℓ}^T·Σ^{-1}·y) / √(s_{v,ℓ}^T·Σ^{-1}·s_{v,ℓ})
κ_ℓ = μ·√(s_{v,ℓ}^T·Σ^{-1}·s_{v,ℓ})
```

**Special cases:**
- Uniform prior π(v,ℓ) = 1/(N|L|) → noncommittal version
- Hard ROI prior: π(v) = 0 outside ROI → ROI-only search
- Spike template s_v = δ_v → prior-weighted max-like detector

### Prior Simulation Results

Power (reject probability) at α ≈ 0.05:

| blob amplitude μ | GRF max (peak) | uniform LR-MF | prior-weighted LR-MF |
|-----------------|----------------|---------------|---------------------|
| 0.2             | 0.051          | 0.063         | 0.081               |
| 0.3             | 0.053          | 0.071         | 0.101               |
| 0.4             | 0.059          | 0.083         | 0.132               |
| 0.5             | 0.073          | 0.100         | 0.168               |

---

## Corrected/Refined Analysis

### Why "Uniform Dominance" Cannot Mean "For All Signals, All Effect Sizes"

For composite alternatives (unknown effect size, unknown number of foci), no single statistic is UMP for all μ. The Neyman-Pearson lemma gives optimality for simple vs simple; once μ is unknown, the optimal rejection region can change with μ.

**Defensible statements:**
- **Point-optimal** for a specified target effect size μ* and signal model
- **Bayes-optimal** after defining a prior over μ
- **Adaptive/omnibus:** maximize over a small grid of μ values

### Multi-Focal Alternative Model

**Biological assumption:** Task fMRI effects are frequently multi-focal (multiple distinct regions respond).

**Model:**

Let y ∈ ℝ^N be a whitened smooth statistic field.

Null: H_0: y ~ N(0, I)

Assume K ≥ 2 non-overlapping foci with common amplitude μ > 0. Centers V_1,...,V_K drawn i.i.d. from prior map π(v):

```
H_1(μ,K): y ~ N(μ·Σ_{k=1}^K s_{V_k}, I),  V_k ~iid π
```

Assuming foci are sufficiently separated (templates approximately orthogonal):

```
Λ_{μ,K}(y) = e^{-(1/2)·μ²·K·||s||²} · (Σ_{v=1}^N π(v)·e^{μ·t_v})^K
```

The NP-optimal test statistic:
```
T_soft(y) := Σ_{v=1}^N π(v)·e^{μ·t_v}
```

**Relationship to RFT hard-max:**
As μ → ∞:
```
(1/μ)·log T_soft(y) → max_v t_v
```
So hard-max (RFT peak detector) is the large-effect / "single dominant peak" limit.

### Theorem (Most Powerful Against Multi-Focal Alternative)

Fix μ > 0 and K ≥ 2. Consider testing H_0 vs the multi-focal mixture alternative H_1(μ,K).

Among all tests with size α, the test that rejects for large Λ_{μ,K}(y) (equivalently large T_soft(y)) has the greatest power. In particular, it has power at least as large as the RFT peak-height hard-max test.

**Proof:** For fixed μ,K, both H_0 and H_1(μ,K) define single distributions for y. By the Neyman-Pearson lemma, the likelihood ratio test is uniformly most powerful among size-α tests for that pair of distributions. ∎

### Fair Simulation (Both Methods Use Same Smoothing)

**Setup:**
- Domain: 64×64
- Null: white noise smoothed by Gaussian kernel (smooth GRF)
- Alternative: add K Gaussian blobs to white noise, then apply same smoothing
- Baseline: hard-max T_max = max y
- Proposed: soft aggregation T_soft = log((1/N)Σ exp(y))
- Both calibrated to FWER ≈ 0.05 by Monte Carlo

**Empirical size (5000 null sims):**
- hard-max: 0.0486
- soft: 0.0496

**Power results (2000 sims per cell):**

| # blobs K | amplitude | hard-max power | soft power |
|-----------|-----------|----------------|------------|
| 1         | 1.2       | 0.282          | 0.251      |
| 1         | 1.6       | 0.635          | 0.562      |
| 2         | 1.0       | 0.245          | 0.268      |
| 2         | 1.4       | 0.652          | 0.672      |
| 3         | 1.0       | 0.319          | 0.371      |
| 3         | 1.4       | 0.756          | 0.836      |

**Interpretation:**
- For single-focal signals, hard-max is often better (tuned to "one dominant peak")
- For multi-focal signals (2-3 blobs), soft aggregation is consistently more powerful

**Corrected claim:** Soft evidence aggregation beats RFT peak-height hard-max under multi-focal alternatives, at equal smoothing and equal FWER. It is not uniformly better for single-focal signals.

---

## Effect-Size Uncertainty: Three Rigorous Options

**Option A: Point-optimal design**
Pick μ* = "smallest effect worth detecting" (biologically motivated). Use T_soft with that μ*. Then NP-optimal at μ*.

**Option B: Bayes-optimal (single statistic)**
Put a prior on μ (e.g., half-normal on standardized effect sizes):
```
T_Bayes(y) = Σ_v π(v) ∫ e^{μ·t_v - (1/2)·μ²·||s||²} dπ(μ)
```

**Option C: Adaptive omnibus over a small grid**
Compute T_soft(μ_j) for a small set {μ_j} and take max_j T_soft(μ_j). Calibrate by permutation/Monte Carlo.

---

## Detection vs Localization

The likelihood ratio is an **omnibus detector** ("is there activation somewhere?"), not a voxelwise "where" test.

### Localization Options

**1. Exploratory posterior localization (most coherent with model)**
```
w(v|y) ∝ π(v)·e^{μ·t_v}
```
Guarantee: coherent model-based localization, but not voxelwise FWER.

**2. Gatekeeper + standard corrected local threshold**
If omnibus test rejects, then localize using conventional corrected peak/cluster method on matched-filtered field t_v.
Guarantee: voxelwise/clusterwise error control.

**3. Full inferential localization**
Requires closed testing / hierarchical multiple testing / selective inference framework.

---

## Summary

**What can be proved cleanly:**
For a biologically plausible multi-focal activation model (multiple bumps; optionally guided by a prior map), the prior-weighted soft evidence aggregation statistic is the most powerful size-α detector (NP lemma). Therefore it is more sensitive than the RFT hard-max test for that model.

**What cannot be honestly claimed:**
A single test statistic is "uniformly more powerful for all μ" against the whole class of unknown-strength alternatives.

**What fair simulations show:**
With equal smoothing and equal FWER, soft aggregation beats hard-max when there are multiple moderate foci; hard-max can win for single-focal signals.
