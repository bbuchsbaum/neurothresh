# Tests for octree max-scan (U0 pyramid)

test_that("pyramid U0 reduces to Z at voxel leaves", {
  z <- 3.2
  pi <- 0.123
  x0 <- 0L; y0 <- 0L; z0 <- 0L

  den <- neurothresh:::build_den_pyramid_xptr_cpp(x0 = x0, y0 = y0, z0 = z0, pi_vec = pi, side = 1L)
  obs <- neurothresh:::pyramid_scores_u0_cpp(Z_vec = z, pi_vec = pi, x0 = x0, y0 = y0, z0 = z0, den_xptr = den)

  leaf <- obs$nodes[obs$nodes$level == 0 & obs$nodes$i == 1 & obs$nodes$j == 1 & obs$nodes$k == 1, ]
  expect_equal(nrow(leaf), 1)
  expect_equal(leaf$score, z, tolerance = 1e-12)
})

test_that("pyramid U0 aggregates correctly for a 2x2x2 block", {
  # 8 voxels at (0/1, 0/1, 0/1)
  coords <- expand.grid(x0 = 0:1, y0 = 0:1, z0 = 0:1)
  x0 <- as.integer(coords$x0)
  y0 <- as.integer(coords$y0)
  z0 <- as.integer(coords$z0)

  Z_vec <- rep(1.0, 8)
  pi_vec <- rep(1 / 8, 8)

  den <- neurothresh:::build_den_pyramid_xptr_cpp(x0, y0, z0, pi_vec, side = 2L)
  obs <- neurothresh:::pyramid_scores_u0_cpp(Z_vec, pi_vec, x0, y0, z0, den)

  # Level 1 root node for side=2
  root <- obs$nodes[obs$nodes$level == 1 & obs$nodes$i == 1 & obs$nodes$j == 1 & obs$nodes$k == 1, ]
  expect_equal(nrow(root), 1)
  expect_equal(root$score, sqrt(8), tolerance = 1e-12)
})
