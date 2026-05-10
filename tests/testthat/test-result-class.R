test_that("new_neurothresh_result builds an object with expected slots", {
  stat_vol <- array(rnorm(8), dim = c(2, 2, 2))
  mask <- array(TRUE, dim = c(2, 2, 2))
  reject <- array(FALSE, dim = c(2, 2, 2))
  pvals <- array(1, dim = c(2, 2, 2))

  r <- neurothresh:::new_neurothresh_result(
    method = "tfce_fwer",
    stat_vol = stat_vol,
    mask = mask,
    threshold = 3.5,
    reject = reject,
    pvals = pvals,
    params = list(alpha = 0.05, tail = "pos", n_perm = 1000),
    null_max = rnorm(1000),
    call = quote(tfce_fwer(stat_vol, alpha = 0.05))
  )

  expect_s3_class(r, "neurothresh_result")
  expect_identical(r$method, "tfce_fwer")
  expect_equal(r$threshold, 3.5)
  expect_identical(dim(r$stat_vol), c(2L, 2L, 2L))
  expect_identical(r$params$alpha, 0.05)
})

test_that("new_neurothresh_result rejects malformed inputs", {
  expect_error(
    neurothresh:::new_neurothresh_result(
      method = c("a", "b"), stat_vol = NULL, mask = NULL,
      threshold = 1, reject = NULL, pvals = NULL
    ),
    "length-1 character"
  )
  expect_error(
    neurothresh:::new_neurothresh_result(
      method = "x", stat_vol = NULL, mask = NULL,
      threshold = c(1, 2), reject = NULL, pvals = NULL
    ),
    "length-1 numeric"
  )
})

test_that("print and summary methods work without error", {
  r <- neurothresh:::new_neurothresh_result(
    method = "rft_peak_fwer",
    stat_vol = array(0, c(2, 2, 2)),
    mask = array(TRUE, c(2, 2, 2)),
    threshold = 4.2,
    reject = array(FALSE, c(2, 2, 2)),
    pvals = array(1, c(2, 2, 2)),
    clusters = data.frame(id = 1:2, size = c(10L, 5L), peak = c(4.5, 4.3)),
    params = list(alpha = 0.05, tail = "pos", smoothness = "fwhm_mm")
  )

  out <- utils::capture.output(print(r))
  expect_true(any(grepl("rft_peak_fwer", out)))
  expect_true(any(grepl("threshold", out)))

  out2 <- utils::capture.output(summary(r))
  expect_true(any(grepl("Clusters", out2)))
  expect_true(any(grepl("Parameters", out2)))
})

test_that("as.data.frame returns the cluster table", {
  cl <- data.frame(id = 1:3, size = c(20L, 12L, 7L), peak = c(5.1, 4.7, 4.3))
  r <- neurothresh:::new_neurothresh_result(
    method = "cluster_fdr",
    stat_vol = NULL, mask = NULL, threshold = 0,
    reject = NULL, pvals = NULL,
    clusters = cl
  )
  df <- as.data.frame(r)
  expect_equal(nrow(df), 3L)
  expect_equal(df$size, c(20L, 12L, 7L))

  r2 <- neurothresh:::new_neurothresh_result(
    method = "tfce_fwer",
    stat_vol = NULL, mask = NULL, threshold = 0,
    reject = NULL, pvals = NULL
  )
  expect_equal(nrow(as.data.frame(r2)), 0L)
})
