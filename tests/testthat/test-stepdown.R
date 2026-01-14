# Tests for Westfall-Young step-down procedure

test_that("wy_stepdown returns correct structure", {
  observed <- c(3.5, 2.1, 4.2, 1.8, 2.9)
  null_matrix <- matrix(rnorm(5000), nrow = 1000, ncol = 5)

  result <- wy_stepdown(observed, null_matrix, alpha = 0.05)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)
  expect_true(all(c("region", "score", "p_adj", "rejected") %in% names(result)))
})

test_that("wy_stepdown p-values are monotonic after sorting", {
  set.seed(42)
  observed <- c(5.0, 4.0, 3.0, 2.0, 1.0)  # Descending order
  null_matrix <- matrix(rnorm(5000), nrow = 1000, ncol = 5)

  result <- wy_stepdown(observed, null_matrix, alpha = 0.05)

  # P-values should be monotonically non-decreasing in original order
  # when sorted by score descending
  ord <- order(result$score, decreasing = TRUE)
  p_sorted <- result$p_adj[ord]

  for (i in 2:length(p_sorted)) {
    expect_true(p_sorted[i] >= p_sorted[i - 1])
  }
})

test_that("wy_stepdown controls FWER under global null", {
  set.seed(123)
  n_sims <- 100
  n_tests <- 10
  n_perm <- 500
  alpha <- 0.05

  false_positives <- 0

  for (sim in seq_len(n_sims)) {
    # All tests under null
    observed <- rnorm(n_tests)
    null_matrix <- matrix(rnorm(n_perm * n_tests), nrow = n_perm, ncol = n_tests)

    result <- wy_stepdown(observed, null_matrix, alpha = alpha)

    if (any(result$rejected)) {
      false_positives <- false_positives + 1
    }
  }

  # FWER should be close to alpha
  observed_fwer <- false_positives / n_sims
  expect_true(observed_fwer <= alpha + 0.03)  # Allow some Monte Carlo error
})

test_that("maxT_singlestep returns correct structure", {
  observed <- c(3.5, 2.1, 4.2)
  null_matrix <- matrix(rnorm(3000), nrow = 1000, ncol = 3)

  result <- maxT_singlestep(observed, null_matrix, alpha = 0.05)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_true(all(c("region", "score", "p_adj", "rejected") %in% names(result)))
})

test_that("wy_stepdown handles empty input", {
  result <- wy_stepdown(numeric(0), matrix(nrow = 100, ncol = 0), alpha = 0.05)
  expect_equal(nrow(result), 0)
})
