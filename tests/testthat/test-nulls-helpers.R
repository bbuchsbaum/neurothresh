# Tests for null distribution helper functions

test_that(".fwhm_mm_to_vox converts correctly", {
  result <- neurothresh:::.fwhm_mm_to_vox(6, 2)
  expect_equal(result, c(3, 3, 3))

  result2 <- neurothresh:::.fwhm_mm_to_vox(c(6, 8, 10), c(2, 2, 2))
  expect_equal(result2, c(3, 4, 5))
})

test_that(".fwhm_mm_to_vox validates lengths", {
  expect_error(neurothresh:::.fwhm_mm_to_vox(c(1, 2), 1), "length 1 or 3")
  expect_error(neurothresh:::.fwhm_mm_to_vox(1, c(1, 2)), "length 1 or 3")
})

test_that(".prep_fwhm_vox normalizes and validates", {
  expect_equal(neurothresh:::.prep_fwhm_vox(3), c(3, 3, 3))
  expect_equal(neurothresh:::.prep_fwhm_vox(c(1, 2, 3)), c(1, 2, 3))
  expect_error(neurothresh:::.prep_fwhm_vox(NULL), "must be provided")
  expect_error(neurothresh:::.prep_fwhm_vox(c(1, 2)), "length 1 or 3")
  expect_error(neurothresh:::.prep_fwhm_vox(c(1, -1, 2)), "finite and >= 0")
})

test_that(".standardize_mask produces zero mean", {
  x <- c(5, 10, 15, 20, 25)
  result <- neurothresh:::.standardize_mask(x)
  expect_equal(mean(result), 0, tolerance = 1e-10)
  expect_equal(sd(result), 1, tolerance = 1e-10)
})

test_that(".standardize_mask handles zero variance", {
  x <- rep(5, 10)
  result <- neurothresh:::.standardize_mask(x)
  expect_true(all(result == 0))
})

test_that(".coerce_null_draw handles full-length vector", {
  mask_idx <- c(2, 5, 8)
  draw <- rep(0.5, 10)

  result <- neurothresh:::.coerce_null_draw(draw, 10, mask_idx)
  expect_length(result, 10)
  expect_equal(result, draw)
})

test_that(".coerce_null_draw handles mask-length vector", {
  mask_idx <- c(2, 5, 8)
  draw <- c(1.0, 2.0, 3.0)

  result <- neurothresh:::.coerce_null_draw(draw, 10, mask_idx)
  expect_length(result, 10)
  expect_equal(result[2], 1.0)
  expect_equal(result[5], 2.0)
  expect_equal(result[8], 3.0)
  expect_equal(result[1], 0)
})

test_that(".coerce_null_draw rejects wrong-length vector", {
  mask_idx <- c(2, 5, 8)
  draw <- c(1.0, 2.0)

  expect_error(
    neurothresh:::.coerce_null_draw(draw, 10, mask_idx),
    "length n_mask or full volume length"
  )
})

test_that(".coerce_null_draw rejects non-finite mask values", {
  mask_idx <- c(2, 5, 8)
  draw <- c(1.0, NaN, 3.0)

  expect_error(
    neurothresh:::.coerce_null_draw(draw, 10, mask_idx),
    "non-finite values"
  )
})

test_that("make_null_fun_voxel_signflip produces correct length", {
  z_vec <- rnorm(100)
  fn <- neurothresh:::make_null_fun_voxel_signflip(z_vec)

  result <- fn(1)
  expect_length(result, 100)
  expect_equal(abs(result), abs(z_vec))
})
