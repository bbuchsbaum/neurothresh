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

