# End-to-end exercise of the canonical API: every public inference method
# called with the new canonical idiom (stat_vol + tail + null_spec/smoothness +
# neurothresh_result return). Designed as a quick smoke + documentation test;
# numerical fixtures live in the per-method test files.

set.seed(2026)

make_canonical_z <- function() {
  z <- array(rnorm(8 * 8 * 8), dim = c(8, 8, 8))
  z[3:5, 3:5, 3:5] <- z[3:5, 3:5, 3:5] + 4
  z
}

test_that("canonicalize -> rft_peak_fwer returns a neurothresh_result", {
  z <- make_canonical_z()
  zvol <- canonicalize(z, stat_type = "Z")
  res <- rft_peak_fwer(zvol, alpha = 0.1, smoothness = smoothness_fwhm(4))
  expect_s3_class(res, "neurothresh_result")
  expect_identical(res$method, "rft_peak_fwer")
  expect_true(is.finite(res$threshold))
  expect_equal(dim(res$reject), dim(z))
})

test_that("rft_cluster_fwer returns clusters and reject volume", {
  z <- make_canonical_z()
  res <- rft_cluster_fwer(z, alpha = 0.2, tail = "pos",
                          cluster_thresh = 3.0,
                          smoothness = smoothness_fwhm(4))
  expect_s3_class(res, "neurothresh_result")
  expect_true(is.data.frame(res$clusters))
})

test_that("tfce_fwer with null_voxel_signflip returns neurothresh_result", {
  z <- make_canonical_z()
  ns <- null_voxel_signflip(z, n_perm = 30, seed = 7)
  res <- tfce_fwer(z, alpha = 0.2, tail = "pos", null_spec = ns, seed = 7)
  expect_s3_class(res, "neurothresh_result")
  expect_equal(length(res$null_max), 30)
  expect_true(is.finite(res$threshold))
})

test_that("cluster_fdr (perm) consumes a null_spec and returns clusters table", {
  z <- make_canonical_z()
  ns <- null_voxel_signflip(z, n_perm = 20, seed = 11)
  res <- cluster_fdr(z, q = 0.2, tail = "pos", cluster_thresh = 3.0,
                    p_method = "perm", null_spec = ns, seed = 11)
  expect_s3_class(res, "neurothresh_result")
  expect_true(is.data.frame(res$clusters))
})

test_that("cluster_fdr (rft) consumes a smoothness object", {
  z <- make_canonical_z()
  res <- cluster_fdr(z, q = 0.2, tail = "pos", cluster_thresh = 3.0,
                    p_method = "rft", smoothness = smoothness_fwhm(4))
  expect_s3_class(res, "neurothresh_result")
  expect_true(is.data.frame(res$clusters))
})

test_that("octree_scan_fwer accepts canonical signature with null_spec", {
  z <- make_canonical_z()
  ns <- null_voxel_signflip(z, n_perm = 30, seed = 13)
  res <- octree_scan_fwer(z, alpha = 0.2, tail = "pos",
                          null_spec = ns, n_perm = 30, seed = 13)
  expect_s3_class(res, "neurothresh_result")
  expect_true(is.finite(res$threshold))
})

test_that("octree_scan_stepdown accepts canonical signature with null_spec", {
  z <- make_canonical_z()
  ns <- null_voxel_signflip(z, n_perm = 30, seed = 17)
  res <- octree_scan_stepdown(z, alpha = 0.2, tail = "pos",
                              min_voxels = 4,
                              null_spec = ns, n_perm = 30, seed = 17)
  expect_s3_class(res, "neurothresh_result")
})

test_that("hier_scan accepts the canonical tail vocabulary", {
  z <- make_canonical_z()
  res <- hier_scan(z, alpha = 0.2, tail = "pos", kappa = 1.0,
                   min_voxels = 4, n_perm = 15, seed = 19)
  expect_s3_class(res, "neurothresh_result")
})

test_that("null_one_sample produces a null_spec consumable by tfce_fwer", {
  set.seed(31)
  subj <- replicate(6, array(rnorm(125), c(5, 5, 5)), simplify = FALSE)
  ns <- null_one_sample(subj, n_perm = 25, seed = 31)
  expect_true(is_null_spec(ns))
  res <- tfce_fwer(ns$stat_vol, mask = ns$mask,
                   alpha = 0.2, tail = "pos", null_spec = ns, seed = 31)
  expect_s3_class(res, "neurothresh_result")
})

test_that("lisa accepts the canonical stat_vol + null_spec API", {
  z <- make_canonical_z()
  ns <- null_voxel_signflip(z, n_perm = 6, seed = 41)
  filter <- lisa_bilateral(radius = 1, sigma_space = 1, sigma_range = 2,
                           n_iter = 1)
  res <- lisa(z, mask = ns$mask, q = 0.2, tail = "two",
              null_spec = ns, n_perm = 6, seed = 41,
              filter = filter)
  expect_s3_class(res, "neurothresh_result")
  expect_identical(res$method, "lisa")
  expect_equal(dim(res$pvals), dim(z))
  expect_identical(res$params$tail, "two")
  expect_equal(res$params$n_perm, 6L)
})
