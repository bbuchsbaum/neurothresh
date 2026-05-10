# Invariant / sanity tests for user-facing thresholding methods

make_blob_vol <- function(dims, blob_corner = c(2, 2, 2), blob_size = c(2, 2, 2),
                          amp = 6, noise_sd = 0) {
  stopifnot(length(dims) == 3)
  z <- array(0, dim = dims)
  if (noise_sd > 0) {
    z[] <- rnorm(prod(dims), sd = noise_sd)
  }

  x <- blob_corner[1]:(blob_corner[1] + blob_size[1] - 1)
  y <- blob_corner[2]:(blob_corner[2] + blob_size[2] - 1)
  k <- blob_corner[3]:(blob_corner[3] + blob_size[3] - 1)
  z[x, y, k] <- z[x, y, k] + amp
  z
}

make_null_fun_iid <- function(n_mask, seed) {
  force(n_mask)
  force(seed)
  function(b) {
    set.seed(seed + b)
    rnorm(n_mask)
  }
}

test_that("alpha monotonicity: stricter alpha does not increase discoveries", {
  set.seed(1)
  z <- array(rnorm(6 * 6 * 6), dim = c(6, 6, 6))

  # octree_scan_fwer threshold should increase as alpha decreases
  res_loose <- octree_scan_fwer(z, alpha = 0.2, n_perm = 20, null = "mc_fwhm", fwhm_vox = 1, seed = 1)
  res_strict <- octree_scan_fwer(z, alpha = 0.05, n_perm = 20, null = "mc_fwhm", fwhm_vox = 1, seed = 1)
  expect_gte(res_strict$u, res_loose$u)

  # stepdown: selected hits should not increase as alpha decreases (same null draws via seed)
  sd_loose <- octree_scan_stepdown(z, alpha = 0.2, n_perm = 20, null = "mc_fwhm", fwhm_vox = 1, min_voxels = 2, seed = 1)
  sd_strict <- octree_scan_stepdown(z, alpha = 0.05, n_perm = 20, null = "mc_fwhm", fwhm_vox = 1, min_voxels = 2, seed = 1)
  expect_lte(nrow(sd_strict$selected), nrow(sd_loose$selected))

  # RFT peak threshold monotonicity
  rft_loose <- rft_peak_fwer(z, smoothness = smoothness_fwhm(8), alpha = 0.2)
  rft_strict <- rft_peak_fwer(z, smoothness = smoothness_fwhm(8), alpha = 0.05)
  expect_gte(rft_strict$threshold, rft_loose$threshold)
})

test_that("seed reproducibility for stochastic methods", {
  z <- array(rnorm(6 * 6 * 6), dim = c(6, 6, 6))

  a <- octree_scan_fwer(z, n_perm = 15, null = "mc_fwhm", fwhm_vox = 1, seed = 123)
  b <- octree_scan_fwer(z, n_perm = 15, null = "mc_fwhm", fwhm_vox = 1, seed = 123)
  expect_identical(a$u, b$u)
  expect_identical(a$M_obs, b$M_obs)

  s1 <- octree_scan_stepdown(z, n_perm = 15, null = "mc_fwhm", fwhm_vox = 1, min_voxels = 2, seed = 123)
  s2 <- octree_scan_stepdown(z, n_perm = 15, null = "mc_fwhm", fwhm_vox = 1, min_voxels = 2, seed = 123)
  expect_identical(nrow(s1$selected), nrow(s2$selected))

  h1 <- hier_scan(z, n_perm = 15, kappa = 1.0, min_voxels = 4, seed = 123)
  h2 <- hier_scan(z, n_perm = 15, kappa = 1.0, min_voxels = 4, seed = 123)
  expect_identical(h1$n_significant, h2$n_significant)

  ns_t <- null_voxel_signflip(z, n_perm = 15, seed = 123)
  t1 <- tfce_fwer(z, null_spec = ns_t, seed = 123)
  t2 <- tfce_fwer(z, null_spec = ns_t, seed = 123)
  expect_identical(as.numeric(t1$threshold), as.numeric(t2$threshold))
})

test_that("two-sided equivalence where implemented as abs()", {
  set.seed(1)
  z <- array(rnorm(6 * 6 * 6), dim = c(6, 6, 6))

  # octree_scan_fwer: two_sided should match abs(z)
  a <- octree_scan_fwer(z, n_perm = 10, null = "mc_fwhm", fwhm_vox = 1, seed = 1, two_sided = TRUE)
  b <- octree_scan_fwer(abs(z), n_perm = 10, null = "mc_fwhm", fwhm_vox = 1, seed = 1, two_sided = FALSE)
  expect_equal(a$M_obs, b$M_obs, tolerance = 1e-12)
  expect_equal(a$u, b$u, tolerance = 1e-12)

  # octree_scan_stepdown: two_sided should match abs(z)
  s_a <- octree_scan_stepdown(z, n_perm = 10, null = "mc_fwhm", fwhm_vox = 1, seed = 1, min_voxels = 2, two_sided = TRUE)
  s_b <- octree_scan_stepdown(abs(z), n_perm = 10, null = "mc_fwhm", fwhm_vox = 1, seed = 1, min_voxels = 2, two_sided = FALSE)
  expect_identical(nrow(s_a$selected), nrow(s_b$selected))

  # hier_scan: two_sided should match abs(z)
  h_a <- hier_scan(z, n_perm = 10, kappa = 1.0, seed = 1, min_voxels = 4, two_sided = TRUE)
  h_b <- hier_scan(abs(z), n_perm = 10, kappa = 1.0, seed = 1, min_voxels = 4, two_sided = FALSE)
  expect_identical(h_a$n_significant, h_b$n_significant)
})

test_that("two-sided wrappers return unions of one-sided results", {
  set.seed(1)
  z <- array(rnorm(8 * 8 * 8), dim = c(8, 8, 8))

  # RFT two-sided: tail="two" matches union of one-sided tails
  s8 <- smoothness_fwhm(8)
  tw <- rft_peak_fwer(z, smoothness = s8, alpha = 0.1, tail = "two")
  pos <- rft_peak_fwer(z, smoothness = s8, alpha = 0.1, tail = "pos")
  neg <- rft_peak_fwer(z, smoothness = s8, alpha = 0.1, tail = "neg")
  expect_identical(as.logical(tw$reject), as.logical(pos$reject) | as.logical(neg$reject))

  # TFCE two-sided: reject equals reject_pos | reject_neg
  ns_tf <- null_voxel_signflip(z, n_perm = 10, seed = 1)
  tf <- tfce_fwer(z, null_spec = ns_tf, alpha = 0.2, seed = 1, tail = "two")
  expect_identical(as.logical(tf$reject),
                   as.logical(tf$extras$reject_pos) | as.logical(tf$extras$reject_neg))
})

test_that("TFCE transform increases inside a strong positive cluster", {
  z <- make_blob_vol(c(10, 10, 10), blob_corner = c(4, 4, 4), blob_size = c(3, 3, 3), amp = 5, noise_sd = 0)
  tf <- tfce_transform(z, tail = "pos", dh = 0.5)

  blob_idx <- which(make_blob_vol(c(10, 10, 10), blob_corner = c(4, 4, 4), blob_size = c(3, 3, 3), amp = 1) > 0)
  bg_idx <- setdiff(which(rep(TRUE, length(z))), blob_idx)

  tf_vec <- as.numeric(tf)
  expect_gt(mean(tf_vec[blob_idx]), mean(tf_vec[bg_idx]))
})

test_that("octree scan detects an obvious blob when given an independent null_fun", {
  set.seed(1)
  z <- make_blob_vol(c(8, 8, 8), blob_corner = c(3, 3, 3), blob_size = c(2, 2, 2), amp = 6, noise_sd = 0.5)
  mask <- is.finite(as.numeric(z))
  n_mask <- sum(mask)
  null_fun <- make_null_fun_iid(n_mask, seed = 9000)

  res <- octree_scan_fwer(z, mask = mask, n_perm = 50, alpha = 0.1, null_fun = null_fun, seed = 1)
  expect_true(res$M_obs >= res$u)

  sd <- octree_scan_stepdown(z, mask = mask, n_perm = 50, alpha = 0.1, null_fun = null_fun, seed = 1, min_voxels = 2)
  expect_gt(nrow(sd$selected), 0)
})

test_that("RFT/cluster-FDR detect an obvious blob", {
  # Large, noiseless blob so inference is deterministic and shouldn't be fragile.
  z <- make_blob_vol(c(20, 20, 20), blob_corner = c(8, 8, 8), blob_size = c(5, 5, 5), amp = 8, noise_sd = 0)

  # Use a modest smoothness so that the cluster is large in resel units.
  s4 <- smoothness_fwhm(4)
  pk <- rft_peak_fwer(z, smoothness = s4, alpha = 0.05, tail = "pos")
  expect_true(any(pk$reject))

  cl <- rft_cluster_fwer(z, smoothness = s4, cluster_thresh = 4.0, alpha = 0.2, tail = "pos")
  expect_true(any(cl$reject))

  fdr <- cluster_fdr(z, cluster_thresh = 4.0, q = 0.2, p_method = "rft", fwhm_mm = 4, tail = "pos")
  expect_true(any(fdr$sig_mask))
})

test_that("canonicalize_stat t->Z matches qnorm(pt())", {
  t_map <- array(c(-2, -1, 0, 1, 2), dim = c(5, 1, 1))
  df <- 20
  out <- canonicalize_stat(t_map, stat_type = "t", df = df, tail = "pos")

  expected <- qnorm(pt(as.numeric(t_map), df = df))
  got <- as.numeric(out$Zeq)
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("hier_scan parallel=TRUE matches serial when future is available", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  z <- array(rnorm(6 * 6 * 6), dim = c(6, 6, 6))

  serial <- hier_scan(z, n_perm = 20, kappa = 1.0, seed = 1, min_voxels = 4, parallel = FALSE)
  # future may emit harmless "built under R version" warnings depending on the local setup
  par <- suppressWarnings(hier_scan(z, n_perm = 20, kappa = 1.0, seed = 1, min_voxels = 4, parallel = TRUE, n_workers = 2))
  expect_identical(serial$n_significant, par$n_significant)
})
