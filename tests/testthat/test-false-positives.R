# Sanity checks for false positive control on null (random) data

estimate_any_fp <- function(n_sims, sim_fun) {
  any_fp <- logical(n_sims)
  for (s in seq_len(n_sims)) {
    set.seed(1000 + s)
    any_fp[s] <- isTRUE(sim_fun(s))
  }
  mean(any_fp)
}

make_smooth_null_array <- function(side, fwhm_vox, seed) {
  grid <- expand.grid(
    x0 = 0:(side - 1),
    y0 = 0:(side - 1),
    z0 = 0:(side - 1)
  )
  x0 <- as.integer(grid$x0)
  y0 <- as.integer(grid$y0)
  z0 <- as.integer(grid$z0)
  mask_full_idx0 <- neurothresh:::.full_index0(x0, y0, z0, side)
  null_fun <- neurothresh:::make_null_fun_mc_fwhm(x0, y0, z0, mask_full_idx0, side, fwhm_vox)

  set.seed(seed)
  z_vec <- null_fun(1)
  z_full <- numeric(side * side * side)
  z_full[mask_full_idx0] <- z_vec
  array(z_full, dim = c(side, side, side))
}

test_that("hier_scan FWER is not wildly anti-conservative under null", {
  alpha <- 0.2
  n_sims <- 15
  rate <- estimate_any_fp(n_sims, function(s) {
    z <- array(rnorm(5 * 5 * 5), dim = c(5, 5, 5))
    res <- hier_scan(
      z,
      alpha = alpha,
      n_perm = 30,
      method = "stepdown",
      kappa = 1.0,
      min_voxels = 4,
      seed = s
    )
    length(res$significant_regions) > 0
  })
  expect_lte(rate, alpha + 0.2)
})

test_that("octree_scan_fwer FWER is not wildly anti-conservative under null", {
  alpha <- 0.2
  n_sims <- 15
  rate <- estimate_any_fp(n_sims, function(s) {
    z <- array(rnorm(6 * 6 * 6), dim = c(6, 6, 6))
    res <- octree_scan_fwer(
      z,
      alpha = alpha,
      n_perm = 30,
      null = "mc_fwhm",
      fwhm_vox = 1,
      seed = s
    )
    res$M_obs >= res$u
  })
  expect_lte(rate, alpha + 0.2)
})

test_that("octree_scan_stepdown FWER is not wildly anti-conservative under null", {
  alpha <- 0.2
  n_sims <- 15
  rate <- estimate_any_fp(n_sims, function(s) {
    z <- array(rnorm(6 * 6 * 6), dim = c(6, 6, 6))
    res <- octree_scan_stepdown(
      z,
      alpha = alpha,
      n_perm = 30,
      min_voxels = 4,
      null = "mc_fwhm",
      fwhm_vox = 1,
      seed = s
    )
    nrow(res$selected) > 0
  })
  expect_lte(rate, alpha + 0.2)
})

test_that("tfce_fwer FWER is not wildly anti-conservative under null", {
  alpha <- 0.2
  n_sims <- 6
  rate <- estimate_any_fp(n_sims, function(s) {
    z <- array(rnorm(5 * 5 * 5), dim = c(5, 5, 5))
    ns <- null_voxel_signflip(z, n_perm = 20, seed = s)
    res <- tfce_fwer(z, null_spec = ns, alpha = alpha, seed = s)
    any(res$reject)
  })
  expect_lte(rate, alpha + 0.35)
})

test_that("rft_peak_fwer is not wildly anti-conservative under null", {
  alpha <- 0.2
  side <- 12
  fwhm_vox <- 3
  fwhm_mm <- fwhm_vox * 2
  n_sims <- 20
  rate <- estimate_any_fp(n_sims, function(s) {
    z <- make_smooth_null_array(side, fwhm_vox, seed = 2000 + s)
    res <- rft_peak_fwer(z, smoothness = smoothness_fwhm(fwhm_mm), alpha = alpha)
    any(res$reject)
  })
  expect_lte(rate, alpha + 0.2)
})

test_that("rft_cluster_fwer is not wildly anti-conservative under null", {
  alpha <- 0.2
  side <- 12
  fwhm_vox <- 3
  fwhm_mm <- fwhm_vox * 2
  n_sims <- 20
  rate <- estimate_any_fp(n_sims, function(s) {
    z <- make_smooth_null_array(side, fwhm_vox, seed = 3000 + s)
    res <- rft_cluster_fwer(z, smoothness = smoothness_fwhm(fwhm_mm),
                            cluster_thresh = 3.0, alpha = alpha)
    any(res$reject)
  })
  expect_lte(rate, alpha + 0.2)
})

test_that("cluster_fdr is not wildly anti-conservative under null", {
  q <- 0.2
  side <- 12
  fwhm_vox <- 3
  fwhm_mm <- fwhm_vox * 2
  n_sims <- 20
  rate <- estimate_any_fp(n_sims, function(s) {
    z <- make_smooth_null_array(side, fwhm_vox, seed = 4000 + s)
    res <- cluster_fdr(
      z,
      cluster_thresh = 3.0,
      q = q,
      p_method = "rft",
      fwhm_mm = fwhm_mm
    )
    any(res$sig_mask)
  })
  expect_lte(rate, q + 0.2)
})
