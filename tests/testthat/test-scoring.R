# Tests for scoring functions

test_that("score_set computes correct log-sum-exp", {
  z_vec <- c(1.0, 2.0, 3.0, 4.0, 5.0)
  pi_vec <- rep(1.0, 5)
  indices <- 1:5
  kappa <- 1.0

  score <- score_set(indices, z_vec, pi_vec, kappa)

  # Manual calculation: log(sum(exp(1:5)))
  expected <- log(sum(exp(1:5)))
  expect_equal(score, expected, tolerance = 1e-10)
})

test_that("score_set handles empty indices", {
  z_vec <- 1:10

pi_vec <- rep(1, 10)

  score <- score_set(integer(0), z_vec, pi_vec, kappa = 1.0)
  expect_equal(score, -Inf)
})

test_that("score_set_stabilized has unit variance under null", {
  # For uniform weights, variance should be ~1 under N(0,1) null
  set.seed(42)

  # Many regions of same size
  n_regions <- 1000
  region_size <- 100
  n_mask <- 10000

  scores <- numeric(n_regions)
  for (i in seq_len(n_regions)) {
    z_vec <- rnorm(n_mask)
    pi_vec <- rep(1, n_mask)
    indices <- sample(n_mask, region_size)
    result <- score_set_stabilized(indices, z_vec, pi_vec)
    scores[i] <- result$U0
  }

  # Variance should be close to 1
  expect_true(abs(var(scores) - 1) < 0.2)
})

test_that("score_set_stabilized n_eff is correct for uniform weights", {
  pi_vec <- rep(1, 100)
  indices <- 1:50

  result <- score_set_stabilized(indices, rnorm(100), pi_vec)

  # For uniform weights, n_eff = n
  expect_equal(result$n_eff, 50)
})

test_that("score_set_omnibus returns max of U0 and S_kappa", {
  set.seed(123)
  z_vec <- rnorm(100)
  pi_vec <- rep(1, 100)
  indices <- 1:50

  result <- score_set_omnibus(indices, z_vec, pi_vec,
                              kappa_grid = c(0.5, 1.0, 2.0))

  expect_true(result$T_omnibus >= result$U0)
  expect_true(result$T_omnibus >= result$S_best || is.infinite(result$S_best))
})
