# Tests for custom null_fun hooks (spatial-covariance-preserving resampling)

test_that("octree_scan_fwer uses user-supplied null_fun", {
  z <- array(rnorm(8), dim = c(2, 2, 2))
  cnt <- 0L

  # null draws: constant zeros (not a valid null, but good for verifying wiring)
  null_fun <- function(b) {
    cnt <<- cnt + 1L
    rep(0, 8)
  }

  res <- octree_scan_fwer(z, n_perm = 5, null_fun = null_fun, seed = 1)
  expect_equal(cnt, 5L)
  expect_equal(length(res$M_null), 5L)
})

test_that("octree_scan_fwer applies tail transforms to external null draws", {
  z <- array(rnorm(8), dim = c(2, 2, 2))
  null_fun <- function(b) rep(-100, length(z))

  res_neg <- octree_scan_fwer(z, tail = "neg", n_perm = 3, null_fun = null_fun)
  expect_true(all(res_neg$M_null > 0))

  res_two <- octree_scan_fwer(z, tail = "two", n_perm = 3, null_fun = null_fun)
  expect_true(all(res_two$M_null > 0))
})

test_that("octree_scan_stepdown uses user-supplied null_fun at root family", {
  z <- array(0, dim = c(4, 4, 4))
  cnt <- 0L
  n_perm <- 7

  null_fun <- function(b) {
    cnt <<- cnt + 1L
    rep(0, length(z))
  }

  res <- octree_scan_stepdown(z, n_perm = n_perm, null_fun = null_fun, seed = 1, min_voxels = 2)
  expect_equal(cnt, n_perm)
  expect_s3_class(res$hits, "data.frame")
})

test_that("octree_scan_stepdown applies tail transforms to external null draws", {
  z <- array(-10, dim = c(4, 4, 4))
  null_fun <- function(b) rep(-100, length(z))

  res <- octree_scan_stepdown(
    z,
    tail = "neg",
    n_perm = 5,
    null_fun = null_fun,
    alpha = 0.5,
    min_voxels = 2
  )

  expect_equal(nrow(res$hits), 0L)
})

test_that("octree null_spec paths default and bound n_perm", {
  z <- array(rnorm(8), dim = c(2, 2, 2))
  cnt <- 0L
  ns <- null_from_function(function(b) {
    cnt <<- cnt + 1L
    rep(0, length(z))
  }, stat_vol = z, n_perm = 3)

  res <- octree_scan_fwer(z, null_spec = ns)
  expect_equal(length(res$M_null), 3L)
  expect_equal(cnt, 3L)
  expect_error(
    octree_scan_fwer(z, null_spec = ns, n_perm = 4),
    "exceeds null_spec"
  )
  expect_error(
    octree_scan_stepdown(z, null_spec = ns, null_fun = function(b) rep(0, length(z))),
    "provide only one"
  )
})

test_that("tfce_fwer consumes draws from a user-supplied null_spec", {
  z <- array(rnorm(4 * 4 * 4), dim = c(4, 4, 4))
  cnt <- 0L
  n_perm <- 6

  null_fun <- function(b) {
    cnt <<- cnt + 1L
    rep(0, sum(is.finite(z)))
  }
  ns <- null_from_function(null_fun, stat_vol = z, n_perm = n_perm)

  res <- tfce_fwer(z, n_perm = n_perm, null_spec = ns, seed = 1, tail = "pos")
  expect_equal(cnt, n_perm)
  expect_equal(length(res$null_max), n_perm)
})

test_that("tfce_fwer requires a null source for calibrated inference", {
  z <- array(rnorm(4 * 4 * 4), dim = c(4, 4, 4))
  expect_error(
    tfce_fwer(z, tail = "pos"),
    "null_fun or null_spec"
  )
})

test_that("cluster_fdr permutation path uses user-supplied null_fun", {
  z <- array(1, dim = c(4, 4, 4))
  cnt <- 0L
  n_perm <- 8

  # Return mask-space draws to exercise coercion logic.
  null_fun <- function(b) {
    cnt <<- cnt + 1L
    rep(0, sum(is.finite(as.numeric(z))))
  }

  res <- cluster_fdr(
    z,
    cluster_thresh = 0.5,
    q = 0.2,
    p_method = "perm",
    n_perm = n_perm,
    null_fun = null_fun,
    tail = "pos",
    seed = 1
  )

  expect_equal(cnt, n_perm)
  expect_true(is.data.frame(res$table))
  expect_true(nrow(res$table) >= 1)
})

test_that("cluster_fdr permutation path requires null_fun unless unsafe fallback is explicitly allowed", {
  z <- array(rnorm(4 * 4 * 4), dim = c(4, 4, 4))
  expect_error(
    cluster_fdr(z, p_method = "perm", n_perm = 4, cluster_thresh = 2.0),
    "requires null_spec"
  )
})
