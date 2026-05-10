# Lightweight contract checks for resampling FWER study utility

test_that("run_resampling_fwer_study returns expected structure", {
  st <- run_resampling_fwer_study(
    n_sims = 2,
    dims = c(6, 6, 6),
    n_subj = 6,
    n_perm = 12,
    alpha = 0.2,
    cluster_thresh = 2.0,
    fwhm_vox = 1.5,
    seed = 9
  )

  expect_true(is.list(st))
  expect_true(all(c("draws", "summary", "params") %in% names(st)))

  expect_s3_class(st$draws, "data.frame")
  expect_equal(nrow(st$draws), 4 * 2)
  expect_true(all(c("sim", "method", "null_model", "reject", "time_sec") %in% names(st$draws)))

  expect_s3_class(st$summary, "data.frame")
  expect_equal(nrow(st$summary), 4)
  expect_true(all(st$summary$reject_rate >= 0 & st$summary$reject_rate <= 1))
  expect_true(all(st$summary$mc_se >= 0))
})

test_that("run_resampling_fwer_study is reproducible with fixed seed", {
  a <- run_resampling_fwer_study(
    n_sims = 2,
    dims = c(5, 5, 5),
    n_subj = 6,
    n_perm = 10,
    alpha = 0.2,
    cluster_thresh = 2.0,
    fwhm_vox = 1.2,
    seed = 42
  )
  b <- run_resampling_fwer_study(
    n_sims = 2,
    dims = c(5, 5, 5),
    n_subj = 6,
    n_perm = 10,
    alpha = 0.2,
    cluster_thresh = 2.0,
    fwhm_vox = 1.2,
    seed = 42
  )

  expect_equal(
    a$draws[, c("sim", "method", "null_model", "reject")],
    b$draws[, c("sim", "method", "null_model", "reject")],
    tolerance = 0
  )
  expect_equal(
    a$summary[, c("method", "null_model", "n_sims", "reject_rate", "mc_se", "ci_low", "ci_high")],
    b$summary[, c("method", "null_model", "n_sims", "reject_rate", "mc_se", "ci_low", "ci_high")],
    tolerance = 0
  )
})
