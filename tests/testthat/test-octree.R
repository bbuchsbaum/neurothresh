# Tests for octree R wrappers

test_that("bbox_from_indices computes correct bounding box", {
  x <- c(1L, 3L, 5L, 2L, 4L)
  y <- c(2L, 4L, 1L, 3L, 5L)
  z <- c(1L, 2L, 3L, 4L, 5L)
  idx <- 1:5

  bb <- bbox_from_indices(idx, x, y, z)

  expect_length(bb, 6)
  expect_equal(bb[1], 1L)
  expect_equal(bb[2], 5L)
  expect_equal(bb[3], 1L)
  expect_equal(bb[4], 5L)
  expect_equal(bb[5], 1L)
  expect_equal(bb[6], 5L)
})

test_that("bbox_from_indices with single voxel", {
  x <- c(3L, 5L)
  y <- c(7L, 9L)
  z <- c(2L, 4L)

  bb <- bbox_from_indices(1L, x, y, z)
  expect_equal(bb[1], 3L)
  expect_equal(bb[2], 3L)
  expect_equal(bb[3], 7L)
  expect_equal(bb[4], 7L)
})

test_that("octree_split_info returns valid structure", {
  n <- 20
  x <- as.integer(rep(1:4, each = 5))
  y <- as.integer(rep(1:5, times = 4))
  z <- as.integer(rep(c(1L, 2L), length.out = n))
  pi_vec <- rep(1.0, n)
  idx_parent <- 1:n
  bb <- bbox_from_indices(idx_parent, x, y, z)

  info <- octree_split_info(idx_parent, bb, x, y, z, pi_vec)

  expect_true(is.list(info))
  expect_true("m" %in% names(info))
  expect_true("child_id" %in% names(info))
  expect_true("w_parent" %in% names(info))
  expect_true("log_den1" %in% names(info))
  expect_true("inv_sqrt_den2" %in% names(info))
  expect_true(info$m >= 1 && info$m <= 8)
  expect_length(info$child_id, length(idx_parent))
  expect_length(info$w_parent, length(idx_parent))
})

test_that("score_children_onepass returns scores for each child", {
  n <- 20
  x <- as.integer(rep(1:4, each = 5))
  y <- as.integer(rep(1:5, times = 4))
  z <- as.integer(rep(c(1L, 2L), length.out = n))
  pi_vec <- rep(1.0, n)
  idx_parent <- 1:n
  bb <- bbox_from_indices(idx_parent, x, y, z)

  info <- octree_split_info(idx_parent, bb, x, y, z, pi_vec)

  Z_vec <- rnorm(n)
  kappa_grid <- c(0.5, 1.0, 2.0)

  scores <- score_children_onepass(Z_vec, idx_parent, info, kappa_grid)

  expect_length(scores, info$m)
  expect_true(all(is.finite(scores)))
})

test_that("split_child_indices partitions correctly", {
  n <- 20
  x <- as.integer(rep(1:4, each = 5))
  y <- as.integer(rep(1:5, times = 4))
  z <- as.integer(rep(c(1L, 2L), length.out = n))
  pi_vec <- rep(1.0, n)
  idx_parent <- 1:n
  bb <- bbox_from_indices(idx_parent, x, y, z)

  info <- octree_split_info(idx_parent, bb, x, y, z, pi_vec)
  children <- split_child_indices(idx_parent, info)

  expect_length(children, info$m)
  all_idx <- sort(unlist(children))
  expect_equal(all_idx, sort(idx_parent))
})

test_that("split_indices_bbox returns indices and bboxes", {
  n <- 20
  x <- as.integer(rep(1:4, each = 5))
  y <- as.integer(rep(1:5, times = 4))
  z <- as.integer(rep(c(1L, 2L), length.out = n))
  pi_vec <- rep(1.0, n)
  idx_parent <- 1:n
  bb <- bbox_from_indices(idx_parent, x, y, z)

  info <- octree_split_info(idx_parent, bb, x, y, z, pi_vec)
  res <- split_indices_bbox(idx_parent, info, x, y, z)

  expect_true("idx" %in% names(res))
  expect_true("bbox" %in% names(res))
  expect_length(res$idx, info$m)
  expect_equal(nrow(res$bbox), info$m)
  expect_equal(ncol(res$bbox), 6)
})

test_that("split_indices_bbox_mass includes prior masses", {
  n <- 20
  x <- as.integer(rep(1:4, each = 5))
  y <- as.integer(rep(1:5, times = 4))
  z <- as.integer(rep(c(1L, 2L), length.out = n))
  pi_vec <- rep(1.0, n)
  idx_parent <- 1:n
  bb <- bbox_from_indices(idx_parent, x, y, z)

  info <- octree_split_info(idx_parent, bb, x, y, z, pi_vec)
  res <- split_indices_bbox_mass(idx_parent, info, x, y, z, pi_vec)

  expect_true("den1" %in% names(res))
  expect_true("den2" %in% names(res))
  expect_length(res$den1, info$m)
  expect_length(res$den2, info$m)
  expect_true(all(res$den1 > 0))
  expect_true(all(res$den2 > 0))
})

test_that("octree_split returns leaf regions", {
  n <- 64
  coords <- expand.grid(x = 1:4, y = 1:4, z = 1:4)
  x <- as.integer(coords$x)
  y <- as.integer(coords$y)
  z <- as.integer(coords$z)
  indices <- seq_len(n)

  leaves <- octree_split(indices, x, y, z, min_voxels = 8)

  expect_true(is.list(leaves))
  expect_true(length(leaves) >= 1)
  all_idx <- sort(unlist(leaves))
  expect_equal(all_idx, seq_len(n))
})

test_that("octree_split with small region returns single leaf", {
  x <- as.integer(c(1, 2, 3))
  y <- as.integer(c(1, 1, 1))
  z <- as.integer(c(1, 1, 1))

  leaves <- octree_split(1:3, x, y, z, min_voxels = 8)

  expect_length(leaves, 1)
  expect_equal(leaves[[1]], 1:3)
})

test_that("octree_split with single-point bbox returns single leaf", {
  x <- as.integer(rep(1L, 10))
  y <- as.integer(rep(1L, 10))
  z <- as.integer(rep(1L, 10))

  leaves <- octree_split(1:10, x, y, z, min_voxels = 2)

  expect_length(leaves, 1)
  expect_equal(leaves[[1]], 1:10)
})
