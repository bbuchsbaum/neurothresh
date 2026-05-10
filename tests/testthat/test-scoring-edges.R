# Tests for scoring function edge cases

test_that("score_set with kappa=0 returns weighted sum", {
  z_vec <- c(1.0, 2.0, 3.0)
  pi_vec <- c(0.5, 0.3, 0.2)
  indices <- 1:3

  score <- score_set(indices, z_vec, pi_vec, kappa = 0)
  expect_equal(score, sum(pi_vec * z_vec))
})

test_that("score_set handles NA values", {
  z_vec <- c(1.0, NA, 3.0)
  pi_vec <- c(1.0, 1.0, 1.0)

  score <- score_set(1:3, z_vec, pi_vec, kappa = 1.0)
  expect_true(is.finite(score))
})

test_that("score_set with all-NA returns -Inf", {
  z_vec <- c(NA, NA)
  pi_vec <- c(1.0, 1.0)

  score <- score_set(1:2, z_vec, pi_vec, kappa = 1.0)
  expect_equal(score, -Inf)
})

test_that("score_set with weighted priors", {
  z_vec <- c(5.0, 0.0, 0.0)
  pi_vec <- c(10.0, 0.001, 0.001)

  score_heavy <- score_set(1:3, z_vec, pi_vec, kappa = 2.0)
  score_uniform <- score_set(1:3, z_vec, rep(1, 3), kappa = 2.0)

  expect_true(score_heavy > score_uniform)
})

test_that("score_set_stabilized handles empty indices", {
  result <- score_set_stabilized(integer(0), rnorm(10), rep(1, 10))
  expect_equal(result$U0, -Inf)
  expect_equal(result$n_eff, 0)
})

test_that("score_set_stabilized handles all-NA", {
  result <- score_set_stabilized(1:3, c(NA, NA, NA), c(1, 1, 1))
  expect_equal(result$U0, -Inf)
  expect_equal(result$n_eff, 0)
})

test_that("score_set_omnibus handles empty indices", {
  result <- score_set_omnibus(integer(0), rnorm(10), rep(1, 10))
  expect_equal(result$T_omnibus, -Inf)
  expect_true(is.na(result$kappa_best))
})

test_that("score_set_omnibus handles all-NA", {
  result <- score_set_omnibus(1:3, c(NA, NA, NA), c(1, 1, 1))
  expect_equal(result$T_omnibus, -Inf)
})

test_that("score_set_omnibus kappa_best is from grid", {
  set.seed(99)
  z_vec <- rnorm(100)
  pi_vec <- rep(1, 100)
  grid <- c(0.5, 1.0, 2.0)

  result <- score_set_omnibus(1:50, z_vec, pi_vec, kappa_grid = grid)
  expect_true(result$kappa_best %in% grid)
})
