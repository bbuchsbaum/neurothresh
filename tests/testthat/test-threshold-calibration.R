# Threshold calibration should match discrete permutation theory

test_that("fwer_threshold_from_maxnull uses expected discrete order statistic", {
  fwer_cut <- neurothresh:::fwer_threshold_from_maxnull
  max_null <- c(10, 9, 8, 7, 6)

  # B = 5, alpha = 0.2 -> k = floor(0.2 * 6) = 1 => largest
  expect_equal(fwer_cut(max_null, alpha = 0.2), 10)

  # alpha = 0.5 -> k = floor(0.5 * 6) = 3 => 3rd largest
  expect_equal(fwer_cut(max_null, alpha = 0.5), 8)
})

test_that("octree_scan_fwer threshold equals discrete cutoff from max null", {
  fwer_cut <- neurothresh:::fwer_threshold_from_maxnull
  set.seed(5)
  z <- array(rnorm(6 * 6 * 6), dim = c(6, 6, 6))
  alpha <- 0.2

  res <- octree_scan_fwer(
    z,
    alpha = alpha,
    n_perm = 25,
    null = "mc_fwhm",
    fwhm_vox = 1,
    seed = 10
  )

  expect_equal(res$u, fwer_cut(res$M_null, alpha), tolerance = 1e-12)
})

test_that("tfce_fwer threshold equals discrete cutoff and reject mask matches corrected p-values", {
  fwer_cut <- neurothresh:::fwer_threshold_from_maxnull
  set.seed(6)
  z <- array(rnorm(6 * 6 * 6), dim = c(6, 6, 6))
  alpha <- 0.2

  ns <- null_voxel_signflip(z, n_perm = 30, seed = 12)
  res <- tfce_fwer(
    z,
    alpha = alpha,
    seed = 12,
    tail = "pos",
    null_spec = ns
  )
  p_corr <- as.numeric(res$pvals)
  rej <- as.logical(res$reject)

  expect_equal(res$threshold, fwer_cut(res$null_max, alpha), tolerance = 1e-12)
  expect_identical(rej, p_corr <= alpha)
})

test_that("rft_peak_fwer computes threshold on 3D arrays without error", {
  z <- array(rnorm(10 * 10 * 10), dim = c(10, 10, 10))
  res <- rft_peak_fwer(z, smoothness = smoothness_fwhm(6), alpha = 0.1)

  expect_s3_class(res, "neurothresh_result")
  expect_true(is.finite(res$threshold))
  expect_true(is.logical(as.logical(res$reject)))
})
