# Tests for stepdown helper functions

test_that("fwer_p_from_maxnull computes correct p-values", {
  max_null <- c(3.0, 2.5, 4.0, 1.0, 3.5)
  observed <- c(3.0, 5.0)

  p <- neurothresh:::fwer_p_from_maxnull(observed, max_null)

  expect_length(p, 2)
  expect_equal(p[1], (1 + sum(max_null >= 3.0)) / (5 + 1))
  expect_equal(p[2], (1 + sum(max_null >= 5.0)) / (5 + 1))
})

test_that("fwer_threshold_from_maxnull computes correct cutoff", {
  max_null <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

  cutoff <- neurothresh:::fwer_threshold_from_maxnull(max_null, alpha = 0.05)

  expect_true(is.numeric(cutoff))
  expect_length(cutoff, 1)
})

test_that("fwer_threshold_from_maxnull validates inputs", {
  expect_error(
    neurothresh:::fwer_threshold_from_maxnull(numeric(0), 0.05),
    "non-empty numeric"
  )
  expect_error(
    neurothresh:::fwer_threshold_from_maxnull(1:10, 0),
    "finite scalar"
  )
  expect_error(
    neurothresh:::fwer_threshold_from_maxnull(1:10, 1),
    "finite scalar"
  )
})

test_that("fwer_threshold_from_maxnull returns Inf for very small alpha", {
  max_null <- 1:10
  cutoff <- neurothresh:::fwer_threshold_from_maxnull(max_null, alpha = 1e-10)
  expect_equal(cutoff, Inf)
})

test_that("maxT_singlestep handles empty input", {
  result <- maxT_singlestep(numeric(0), matrix(nrow = 100, ncol = 0), alpha = 0.05)
  expect_equal(nrow(result), 0)
})

test_that("maxT_singlestep rejects strong signals", {
  set.seed(42)
  observed <- c(10.0, 0.5, 0.1)
  null_matrix <- matrix(rnorm(3000), nrow = 1000, ncol = 3)

  result <- maxT_singlestep(observed, null_matrix, alpha = 0.05)

  expect_true(result$rejected[1])
  expect_true(result$p_adj[1] < 0.05)
})
