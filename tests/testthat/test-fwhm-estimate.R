# Tests for map-only FWHM estimator

test_that("estimate_fwhm_vox returns finite values on a smooth null field", {
  set.seed(1)
  dims <- c(20L, 20L, 20L)
  n <- prod(dims)
  eps <- rnorm(n)

  # Generate a smooth field by Gaussian smoothing; note generator uses
  # sigma_kernel = sigma_field / sqrt(2).
  fwhm_true <- rep(4, 3)
  sigma_field <- fwhm_true / sqrt(8 * log(2))
  sigma_kernel <- sigma_field / sqrt(2)

  sm <- neurothresh:::gaussian_smooth3d_cpp(eps, dims = dims, sigma_xyz = sigma_kernel)
  z <- array(sm, dim = dims)
  z <- (z - mean(z)) / sd(z)

  fwhm_hat <- estimate_fwhm_vox(z, isotropic = FALSE, z_trim = 1.5)
  expect_true(all(is.finite(fwhm_hat)))
  # Loose tolerance; single realization + robust sign estimator.
  expect_true(all(abs(fwhm_hat - fwhm_true) / fwhm_true < 0.6))
})
