# Package index

## All functions

- [`as.data.frame(`*`<neurothresh_result>`*`)`](https://bbuchsbaum.github.io/neurothresh/reference/as.data.frame.neurothresh_result.md)
  : Coerce a neurothresh result's clusters to a data frame
- [`as_subject_matrix()`](https://bbuchsbaum.github.io/neurothresh/reference/as_subject_matrix.md)
  : Coerce subject-level images to a canonical (subject by voxel) matrix
- [`bandpass_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/bandpass_energy.md)
  : Band-Pass Energy Statistic Channel
- [`bbox_from_indices()`](https://bbuchsbaum.github.io/neurothresh/reference/bbox_from_indices.md)
  : Compute bounding box from mask-space indices
- [`canonicalize()`](https://bbuchsbaum.github.io/neurothresh/reference/canonicalize.md)
  : Canonicalize a statistic map to a Z-equivalent volume
- [`canonicalize_stat()`](https://bbuchsbaum.github.io/neurothresh/reference/canonicalize_stat.md)
  : Canonicalize statistical maps to Z-scores
- [`cluster-fdr`](https://bbuchsbaum.github.io/neurothresh/reference/cluster-fdr.md)
  : Cluster-Level FDR Control
- [`cluster_fdr()`](https://bbuchsbaum.github.io/neurothresh/reference/cluster_fdr.md)
  : Cluster-FDR with RFT or Permutation P-values
- [`energy_mass_map()`](https://bbuchsbaum.github.io/neurothresh/reference/energy_mass_map.md)
  : Above-Null Energy Mass Map
- [`estimate_fwhm_vox()`](https://bbuchsbaum.github.io/neurothresh/reference/estimate_fwhm_vox.md)
  : Estimate Gaussian ACF FWHM from a single map (map-only)
- [`generate_null_scores()`](https://bbuchsbaum.github.io/neurothresh/reference/generate_null_scores.md)
  : Generate null distribution via sign-flipping
- [`generate_null_scores_parcels()`](https://bbuchsbaum.github.io/neurothresh/reference/generate_null_scores_parcels.md)
  : Generate null scores for parcel-level analysis
- [`global_band_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/global_band_energy.md)
  : Global Band-Pass Energy Statistic Channel
- [`global_lowpass_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/global_lowpass_energy.md)
  : Global Low-Pass Energy Statistic Channel
- [`hier_scan()`](https://bbuchsbaum.github.io/neurothresh/reference/hier_scan.md)
  : Hierarchical LR-MFT Analysis
- [`hierarchical`](https://bbuchsbaum.github.io/neurothresh/reference/hierarchical.md)
  : Hierarchical Testing with Alpha-Spending
- [`is_lisa_filter()`](https://bbuchsbaum.github.io/neurothresh/reference/is_lisa_filter.md)
  : Test whether an object is a LISA filter specification
- [`is_null_spec()`](https://bbuchsbaum.github.io/neurothresh/reference/is_null_spec.md)
  : Test whether an object is a neurothresh_null_spec
- [`is_smoothness()`](https://bbuchsbaum.github.io/neurothresh/reference/is_smoothness.md)
  : Test whether an object is a smoothness specification
- [`lisa()`](https://bbuchsbaum.github.io/neurothresh/reference/lisa.md)
  : LISA local adaptive inference
- [`lisa_apply_filter()`](https://bbuchsbaum.github.io/neurothresh/reference/lisa_apply_filter.md)
  : Apply a LISA filter
- [`lisa_bilateral()`](https://bbuchsbaum.github.io/neurothresh/reference/lisa_bilateral.md)
  : LISA bilateral filter specification
- [`lisa_fdr_map()`](https://bbuchsbaum.github.io/neurothresh/reference/lisa_fdr_map.md)
  : Empirical FDR for LISA filtered maps
- [`lisa_filtered_stat()`](https://bbuchsbaum.github.io/neurothresh/reference/lisa_filtered_stat.md)
  : Extract the LISA filtered statistic map
- [`lisa_glm()`](https://bbuchsbaum.github.io/neurothresh/reference/lisa_glm.md)
  : LISA GLM wrapper
- [`lisa_onesample()`](https://bbuchsbaum.github.io/neurothresh/reference/lisa_onesample.md)
  : LISA one-sample wrapper
- [`lisa_twosample()`](https://bbuchsbaum.github.io/neurothresh/reference/lisa_twosample.md)
  : LISA two-sample wrapper
- [`local_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/local_energy.md)
  : Local Energy Statistic Channel
- [`lowpass_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/lowpass_energy.md)
  : Low-Pass Energy Statistic Channel
- [`make_null_fun_glm()`](https://bbuchsbaum.github.io/neurothresh/reference/make_null_fun_glm.md)
  : GLM Freedman-Lane Null Generator
- [`make_null_fun_onesample()`](https://bbuchsbaum.github.io/neurothresh/reference/make_null_fun_onesample.md)
  : Subject-Level One-Sample Sign-Flip Null Generator
- [`make_null_fun_subject_signflip()`](https://bbuchsbaum.github.io/neurothresh/reference/make_null_fun_subject_signflip.md)
  : Subject-level sign-flip null generator (one-sample)
- [`make_null_fun_twosample()`](https://bbuchsbaum.github.io/neurothresh/reference/make_null_fun_twosample.md)
  : Two-Group Label-Permutation Null Generator
- [`make_perm_fun()`](https://bbuchsbaum.github.io/neurothresh/reference/make_perm_fun.md)
  : Create permutation function for hierarchical analysis
- [`make_within_subject_contrast_images()`](https://bbuchsbaum.github.io/neurothresh/reference/make_within_subject_contrast_images.md)
  : Construct Within-Subject Contrast Images
- [`maxT_singlestep()`](https://bbuchsbaum.github.io/neurothresh/reference/maxT_singlestep.md)
  : Single-step MaxT Procedure
- [`null_from_function()`](https://bbuchsbaum.github.io/neurothresh/reference/null_from_function.md)
  : Wrap a bare draw function in a null_spec
- [`null_from_maps()`](https://bbuchsbaum.github.io/neurothresh/reference/null_from_maps.md)
  : External observed/null statistic maps
- [`null_glm()`](https://bbuchsbaum.github.io/neurothresh/reference/null_glm.md)
  : GLM Freedman-Lane null
- [`null_mc_acf()`](https://bbuchsbaum.github.io/neurothresh/reference/null_mc_acf.md)
  : Monte-Carlo ACF null (AFNI-style autocorrelation)
- [`null_mc_fwhm()`](https://bbuchsbaum.github.io/neurothresh/reference/null_mc_fwhm.md)
  : Monte-Carlo FWHM null (Gaussian-smoothed white noise)
- [`null_one_sample()`](https://bbuchsbaum.github.io/neurothresh/reference/null_one_sample.md)
  : One-sample subject sign-flip null
- [`null_parcel_perm()`](https://bbuchsbaum.github.io/neurothresh/reference/null_parcel_perm.md)
  : Parcel-level permutation null
- [`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md)
  : Null specification S3 class
- [`null_subject_signflip()`](https://bbuchsbaum.github.io/neurothresh/reference/null_subject_signflip.md)
  : Subject sign-flip null (alias of one-sample)
- [`null_two_sample()`](https://bbuchsbaum.github.io/neurothresh/reference/null_two_sample.md)
  : Two-sample label-permutation null
- [`null_voxel_signflip()`](https://bbuchsbaum.github.io/neurothresh/reference/null_voxel_signflip.md)
  : Voxel-level sign-flip null
- [`octree-wrappers`](https://bbuchsbaum.github.io/neurothresh/reference/octree-wrappers.md)
  : R wrappers for octree operations
- [`octree_scan_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/octree_scan_fwer.md)
  : Octree Max-Scan with permutation FWER (U0-only)
- [`octree_scan_stepdown()`](https://bbuchsbaum.github.io/neurothresh/reference/octree_scan_stepdown.md)
  : Octree scan with hierarchical step-down inference (U0-only)
- [`octree_split()`](https://bbuchsbaum.github.io/neurothresh/reference/octree_split.md)
  : Recursive octree split
- [`octree_split_info()`](https://bbuchsbaum.github.io/neurothresh/reference/octree_split_info.md)
  : Build octree split info for a parent node
- [`permutation_glm()`](https://bbuchsbaum.github.io/neurothresh/reference/permutation_glm.md)
  : GLM Permutation Specification
- [`plot_resolution_curve()`](https://bbuchsbaum.github.io/neurothresh/reference/plot_resolution_curve.md)
  : Resolution Evidence Curve
- [`plot_scale_space_map()`](https://bbuchsbaum.github.io/neurothresh/reference/plot_scale_space_map.md)
  : Extract a Scale-Space Corrected Map
- [`positive_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/positive_energy.md)
  : Positive-Direction Local Energy Statistic Channel
- [`print(`*`<neurothresh_null_spec>`*`)`](https://bbuchsbaum.github.io/neurothresh/reference/print.neurothresh_null_spec.md)
  : Print a neurothresh null specification
- [`print(`*`<neurothresh_result>`*`)`](https://bbuchsbaum.github.io/neurothresh/reference/print.neurothresh_result.md)
  : Print a neurothresh result
- [`print(`*`<neurothresh_smoothness>`*`)`](https://bbuchsbaum.github.io/neurothresh/reference/print.neurothresh_smoothness.md)
  : Print a smoothness specification
- [`rft-methods`](https://bbuchsbaum.github.io/neurothresh/reference/rft-methods.md)
  : Random Field Theory (RFT) Thresholding Methods
- [`rft_cluster_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/rft_cluster_fwer.md)
  : RFT Cluster-Level FWER
- [`rft_peak_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/rft_peak_fwer.md)
  : RFT Peak-Level FWER
- [`run_resampling_fwer_study()`](https://bbuchsbaum.github.io/neurothresh/reference/run_resampling_fwer_study.md)
  : Simulation Study: Resampling Choice and Null Rejection Rates
- [`scale_space()`](https://bbuchsbaum.github.io/neurothresh/reference/scale_space.md)
  : Define a Gaussian Scale-Space Grid
- [`scale_space_diagnostics()`](https://bbuchsbaum.github.io/neurothresh/reference/scale_space_diagnostics.md)
  : Diagnose Scale-Space Evidence Curves
- [`scale_space_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/scale_space_energy.md)
  : Scale-Space Localized Energy Inference
- [`scale_space_onesample()`](https://bbuchsbaum.github.io/neurothresh/reference/scale_space_onesample.md)
  : Scale-Space One-Sample Wrapper
- [`scale_space_twosample()`](https://bbuchsbaum.github.io/neurothresh/reference/scale_space_twosample.md)
  : Scale-Space Two-Sample Wrapper
- [`score_children_onepass()`](https://bbuchsbaum.github.io/neurothresh/reference/score_children_onepass.md)
  : Score all children in one pass through parent voxels
- [`score_set()`](https://bbuchsbaum.github.io/neurothresh/reference/score_set.md)
  : Compute prior-weighted soft evidence statistic T_kappa(R)
- [`score_set_omnibus()`](https://bbuchsbaum.github.io/neurothresh/reference/score_set_omnibus.md)
  : Compute omnibus score combining diffuse and focal detection
- [`score_set_stabilized()`](https://bbuchsbaum.github.io/neurothresh/reference/score_set_stabilized.md)
  : Compute variance-stabilized score U_0(R)
- [`signed_coherence()`](https://bbuchsbaum.github.io/neurothresh/reference/signed_coherence.md)
  : Signed Coherence Statistic Channel
- [`smoothness_acf()`](https://bbuchsbaum.github.io/neurothresh/reference/smoothness_acf.md)
  : AFNI-style ACF smoothness specification
- [`smoothness_fwhm()`](https://bbuchsbaum.github.io/neurothresh/reference/smoothness_fwhm.md)
  : FWHM-based smoothness specification
- [`split_child_indices()`](https://bbuchsbaum.github.io/neurothresh/reference/split_child_indices.md)
  : Split parent indices into child index lists
- [`split_indices_bbox()`](https://bbuchsbaum.github.io/neurothresh/reference/split_indices_bbox.md)
  : Split indices with bounding boxes
- [`split_indices_bbox_mass()`](https://bbuchsbaum.github.io/neurothresh/reference/split_indices_bbox_mass.md)
  : Split indices with bounding boxes and prior masses
- [`summary(`*`<neurothresh_result>`*`)`](https://bbuchsbaum.github.io/neurothresh/reference/summary.neurothresh_result.md)
  : Summarize a neurothresh result
- [`tfce-methods`](https://bbuchsbaum.github.io/neurothresh/reference/tfce-methods.md)
  : Threshold-Free Cluster Enhancement (TFCE)
- [`tfce_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/tfce_fwer.md)
  : TFCE with Permutation FWER
- [`tfce_transform()`](https://bbuchsbaum.github.io/neurothresh/reference/tfce_transform.md)
  : TFCE Transform
- [`validate_null_spec()`](https://bbuchsbaum.github.io/neurothresh/reference/validate_null_spec.md)
  : Validate a neurothresh_null_spec
- [`wy_stepdown()`](https://bbuchsbaum.github.io/neurothresh/reference/wy_stepdown.md)
  : Westfall-Young Step-Down Procedure
