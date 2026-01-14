test_that("prep_prior normalizes and mixes with uniform", {
  pi_vec <- c(0, 1, 1, NA_real_)
  out <- neurothresh:::.prep_prior(pi_vec, eta = 0.5)

  expect_equal(sum(out), 1, tolerance = 1e-12)
  expect_true(all(out > 0))

  uniform <- neurothresh:::.prep_prior(pi_vec, eta = 0)
  expect_equal(uniform, rep(1 / length(pi_vec), length(pi_vec)))
})

test_that("array_to_vol preserves dims for array template", {
  template <- array(0, dim = c(2, 2, 2))
  vals <- seq_len(8)

  out <- neurothresh:::.array_to_vol(vals, dim(template), template)
  expect_equal(dim(out), dim(template))
  expect_equal(as.vector(out), vals)
})
