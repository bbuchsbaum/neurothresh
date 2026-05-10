test_that("new_null_spec validates inputs", {
  expect_error(neurothresh:::new_null_spec(NULL, NULL, draw = "x", n_perm = 1),
               "draw must be a function")
  expect_error(neurothresh:::new_null_spec(NULL, NULL, draw = function(b) 0,
                                            n_perm = 0),
               "positive integer scalar")
  expect_error(neurothresh:::new_null_spec(NULL, NULL, draw = function(b) 0,
                                            n_perm = 1, kind = c("a", "b")),
               "length-1 character")
})

test_that("is_null_spec recognises spec objects", {
  expect_true(is_null_spec(null_voxel_signflip(array(rnorm(8), c(2, 2, 2)),
                                                n_perm = 5, seed = 1)))
  expect_false(is_null_spec(list()))
  expect_false(is_null_spec(NULL))
})

test_that("null_voxel_signflip draw is deterministic by index", {
  set.seed(1)
  z <- array(rnorm(27), c(3, 3, 3))
  ns <- null_voxel_signflip(z, n_perm = 10, seed = 7)
  expect_identical(ns$draw(3), ns$draw(3))
  expect_equal(length(ns$draw(1)), sum(is.finite(z)))
  expect_error(ns$draw(11), "out of range")
})

test_that("Monte Carlo null specs are deterministic by index and bounded", {
  z <- array(rnorm(27), c(3, 3, 3))
  mask <- array(TRUE, c(3, 3, 3))

  ns_f <- null_mc_fwhm(z, mask = mask, smoothness = 2, n_perm = 4, seed = 11)
  expect_identical(ns_f$draw(2), ns_f$draw(2))
  expect_equal(length(ns_f$draw(1)), sum(mask))
  expect_error(ns_f$draw(5), "Permutation index out of range")

  ns_a <- null_mc_acf(
    z, mask = mask,
    smoothness = smoothness_acf(a = 0.7, b = 3, c = 8),
    n_perm = 4, seed = 12
  )
  expect_identical(ns_a$draw(3), ns_a$draw(3))
  expect_equal(length(ns_a$draw(1)), sum(mask))
  expect_error(ns_a$draw(5), "Permutation index out of range")
})

test_that("null_one_sample produces null_spec from list and 4D array equally", {
  set.seed(1)
  arrs <- replicate(5, array(rnorm(27), c(3, 3, 3)), simplify = FALSE)
  arr4d <- array(NA_real_, c(3, 3, 3, 5))
  for (i in seq_along(arrs)) arr4d[, , , i] <- arrs[[i]]

  a <- null_one_sample(arrs, n_perm = 8, seed = 2)
  b <- null_one_sample(arr4d, n_perm = 8, seed = 2)

  expect_true(is_null_spec(a))
  expect_equal(a$df, 4)
  expect_identical(a$draw(3), b$draw(3))
  expect_equal(length(a$draw(1)), sum(a$mask))
})

test_that("null_subject_signflip is identical to null_one_sample", {
  arrs <- replicate(4, array(rnorm(8), c(2, 2, 2)), simplify = FALSE)
  a <- null_one_sample(arrs, n_perm = 6, seed = 3)
  b <- null_subject_signflip(arrs, n_perm = 6, seed = 3)
  expect_identical(a$draw(2), b$draw(2))
  expect_identical(a$kind, "one_sample")
})

test_that("null_two_sample uses level-2-minus-level-1 by default", {
  set.seed(7)
  arrs <- replicate(6, array(rnorm(8), c(2, 2, 2)), simplify = FALSE)
  ns <- null_two_sample(arrs, group = c("a", "a", "a", "b", "b", "b"),
                        n_perm = 5, seed = 1)
  expect_true(is_null_spec(ns))
  expect_identical(ns$kind, "two_sample")
  expect_match(ns$meta$contrast, "b minus a")
})

test_that("null_two_sample contrast= flips direction", {
  arrs <- replicate(6, array(rnorm(8), c(2, 2, 2)), simplify = FALSE)
  default <- null_two_sample(arrs, group = c("a", "a", "a", "b", "b", "b"),
                             n_perm = 3, seed = 1)
  flipped <- null_two_sample(arrs, group = c("a", "a", "a", "b", "b", "b"),
                             contrast = c("a", "b"),
                             n_perm = 3, seed = 1)
  expect_equal(as.numeric(default$stat_vol), -as.numeric(flipped$stat_vol),
               tolerance = 1e-12)
})

test_that("null_two_sample rejects unknown explicit group levels", {
  arrs <- replicate(4, array(rnorm(8), c(2, 2, 2)), simplify = FALSE)
  expect_error(
    null_two_sample(arrs, group = c("a", "a", "b", "c"), levels = c("a", "b")),
    "not present in levels"
  )
})

test_that("null_parcel_perm tags kind correctly", {
  z <- array(rnorm(8), c(2, 2, 2))
  ns <- null_parcel_perm(z, n_perm = 3, seed = 1)
  expect_identical(ns$kind, "parcel_perm")
})

test_that("validate_null_spec catches inconsistent draw outputs", {
  z <- array(rnorm(8), c(2, 2, 2))
  ns <- null_voxel_signflip(z, n_perm = 3, seed = 1)
  expect_silent(validate_null_spec(ns))

  bad <- ns
  bad$draw <- function(b) numeric(0)
  expect_error(validate_null_spec(bad), "expected")
})
