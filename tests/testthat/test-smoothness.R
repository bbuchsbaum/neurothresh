test_that("smoothness_fwhm constructs and validates", {
  s <- smoothness_fwhm(6)
  expect_s3_class(s, "neurothresh_smoothness")
  expect_s3_class(s, "neurothresh_smoothness_fwhm")
  expect_equal(s$fwhm_mm, c(6, 6, 6))

  s3 <- smoothness_fwhm(c(6, 6, 8))
  expect_equal(s3$fwhm_mm, c(6, 6, 8))

  expect_error(smoothness_fwhm("a"), "must be numeric")
  expect_error(smoothness_fwhm(c(6, 6)), "length 1 or 3")
  expect_error(smoothness_fwhm(0), "strictly positive")
  expect_error(smoothness_fwhm(-1), "strictly positive")
  expect_error(smoothness_fwhm(NA_real_), "strictly positive")
})

test_that("smoothness_acf constructs and validates", {
  s <- smoothness_acf(a = 0.7, b = 3.0, c = 8.0)
  expect_s3_class(s, "neurothresh_smoothness")
  expect_s3_class(s, "neurothresh_smoothness_acf")
  expect_equal(s$a, 0.7)
  expect_equal(s$b, 3.0)
  expect_equal(s$c, 8.0)

  expect_error(smoothness_acf(0, 1, 1), "between 0 and 1")
  expect_error(smoothness_acf(1, 1, 1), "between 0 and 1")
  expect_error(smoothness_acf(0.5, 0, 1), "strictly positive")
  expect_error(smoothness_acf(0.5, 1, -1), "strictly positive")
  expect_error(smoothness_acf(c(0.5, 0.5), 1, 1), "length-1")
})

test_that("is_smoothness predicate", {
  expect_true(is_smoothness(smoothness_fwhm(6)))
  expect_true(is_smoothness(smoothness_acf(0.5, 1, 2)))
  expect_false(is_smoothness(list(fwhm_mm = 6)))
  expect_false(is_smoothness(6))
})

test_that("print methods work", {
  out_f <- utils::capture.output(print(smoothness_fwhm(c(6, 6, 8))))
  expect_true(any(grepl("smoothness_fwhm", out_f)))
  expect_true(any(grepl("6.+8", out_f)))

  out_a <- utils::capture.output(print(smoothness_acf(0.7, 3, 8)))
  expect_true(any(grepl("smoothness_acf", out_a)))
})
