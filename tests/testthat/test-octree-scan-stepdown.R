# Tests for octree_scan_stepdown (U0-only hierarchical step-down)

test_that("octree_scan_stepdown returns no hits for all-zero map", {
  z <- array(0, dim = c(4, 4, 4))
  res <- octree_scan_stepdown(z, n_perm = 30, seed = 1, min_voxels = 1,
                              null = "mc_fwhm", fwhm_vox = 1)
  expect_s3_class(res$hits, "data.frame")
  expect_equal(nrow(res$hits), 0)
  expect_equal(length(res$sig_indices), 0)
})

test_that("node selection modes behave sensibly on a nested hit table", {
  hits <- data.frame(
    node_id = c(1L, 2L, 3L, 4L),
    parent_id = c(NA_integer_, 1L, 1L, NA_integer_),
    depth = c(1L, 2L, 2L, 1L),
    score = c(5, 4, 3, 4.5),
    p_adj = c(0.01, 0.02, 0.03, 0.01),
    alpha_test = c(0.05, 0.05, 0.05, 0.05),
    n_voxels = c(8L, 1L, 1L, 8L),
    mass = c(0.2, 0.025, 0.025, 0.2),
    x0 = c(1L, 1L, 2L, 3L), x1 = c(2L, 1L, 2L, 4L),
    y0 = c(1L, 1L, 2L, 3L), y1 = c(2L, 1L, 2L, 4L),
    z0 = c(1L, 1L, 2L, 3L), z1 = c(2L, 1L, 2L, 4L),
    indices = I(list(1:8, 1L, 8L, 9:16))
  )

  co <- neurothresh:::.select_nodes(hits, report = "coarsest")
  expect_equal(sort(co$node_id), c(1L, 4L))

  fi <- neurothresh:::.select_nodes(hits, report = "finest")
  expect_equal(sort(fi$node_id), c(2L, 3L, 4L))

  gr <- neurothresh:::.select_nodes(hits, report = "greedy")
  expect_true(all(gr$node_id %in% c(1L, 4L)))
})

test_that("octree_scan_stepdown returns expected structure on random input", {
  z <- array(rnorm(4 * 4 * 4), dim = c(4, 4, 4))
  res <- octree_scan_stepdown(z, n_perm = 10, seed = 1, min_voxels = 2,
                              null = "mc_fwhm", fwhm_vox = 1)
  expect_true(all(c("hits", "selected", "sig_mask", "sig_indices", "params") %in% names(res)))
  expect_equal(dim(res$sig_mask), dim(z))
})

test_that("octree_scan_stepdown is not wildly anti-conservative under null (small Monte Carlo)", {
  set.seed(123)
  n_sims <- 10
  alpha <- 0.2
  n_perm <- 40

  any_fp <- logical(n_sims)
  for (s in seq_len(n_sims)) {
    z <- array(rnorm(4 * 4 * 4), dim = c(4, 4, 4))
    res <- octree_scan_stepdown(z, alpha = alpha, n_perm = n_perm, seed = s, min_voxels = 2,
                                null = "mc_fwhm", fwhm_vox = 1)
    any_fp[s] <- nrow(res$selected) > 0
  }

  # Loose guardrail: should not be grossly above alpha with this tiny setup.
  expect_lte(mean(any_fp), 0.6)
})
