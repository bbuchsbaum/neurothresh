# Contract/differential/metamorphic tests for permutation utilities

test_that("generate_null_scores returns expected shape and finite values", {
  z <- c(1.0, -0.5, 0.2, 2.1, -1.3, 0.8)
  pi <- c(1, 2, 1, 1, 3, 2)
  regions <- list(1:2, 3:4, 5:6)

  out <- generate_null_scores(
    z_vec = z,
    pi_vec = pi,
    regions = regions,
    n_perm = 15,
    kappa_grid = c(0.5, 1.0, 2.0),
    seed = 11
  )

  expect_equal(dim(out), c(15, 3))
  expect_true(all(is.finite(out)))
})

test_that("parcel vectorized null scores match region-loop null scores", {
  z <- c(1.2, -0.3, 0.4, 2.5, -1.1, 0.9, -0.2, 1.7)
  pi <- c(1, 1, 2, 2, 1, 3, 1, 2)
  lab <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)
  regions <- lapply(sort(unique(lab)), function(p) which(lab == p))

  loop_scores <- generate_null_scores(
    z_vec = z,
    pi_vec = pi,
    regions = regions,
    n_perm = 25,
    kappa_grid = c(0.5, 1.0, 2.0),
    seed = 19
  )
  vec_scores <- generate_null_scores_parcels(
    z_vec = z,
    pi_vec = pi,
    lab_vec = lab,
    n_perm = 25,
    kappa_grid = c(0.5, 1.0, 2.0),
    seed = 19
  )

  expect_equal(vec_scores, loop_scores, tolerance = 1e-12)
})

test_that("region-order permutation only permutes output columns", {
  z <- rnorm(10)
  pi <- abs(rnorm(10)) + 0.2
  regions <- list(1:3, 4:7, 8:10)
  perm_idx <- c(3, 1, 2)

  base <- generate_null_scores(
    z_vec = z,
    pi_vec = pi,
    regions = regions,
    n_perm = 20,
    kappa_grid = c(0.5, 1.0),
    seed = 101
  )
  permuted <- generate_null_scores(
    z_vec = z,
    pi_vec = pi,
    regions = regions[perm_idx],
    n_perm = 20,
    kappa_grid = c(0.5, 1.0),
    seed = 101
  )

  expect_equal(permuted, base[, perm_idx, drop = FALSE], tolerance = 1e-12)
})

test_that("seed reproducibility holds for parcel null generation", {
  z <- rnorm(12)
  pi <- abs(rnorm(12)) + 0.1
  lab <- rep(1:3, each = 4)

  a <- generate_null_scores_parcels(z, pi, lab, n_perm = 20, seed = 77)
  b <- generate_null_scores_parcels(z, pi, lab, n_perm = 20, seed = 77)
  expect_identical(a, b)
})

test_that("make_perm_fun is deterministic per index and checks bounds", {
  z <- c(1, -2, 3, -4)
  pf <- make_perm_fun(z, n_perm = 5, seed = 7)

  expect_equal(pf(3), pf(3), tolerance = 0)
  expect_error(pf(0), "Permutation index out of range")
  expect_error(pf(6), "Permutation index out of range")
})
