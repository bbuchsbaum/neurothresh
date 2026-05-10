test_that("canonicalize() returns just the Zeq stat_vol", {
  set.seed(1)
  z <- array(rnorm(8), c(2, 2, 2))

  out <- canonicalize(z, stat_type = "Z")
  expect_true(is.array(out) || inherits(out, "NeuroVol"))
  expect_equal(dim(out), c(2, 2, 2))
  expect_equal(as.numeric(out), as.numeric(z))
})

test_that("canonicalize() converts t-maps via probability integral transform", {
  set.seed(2)
  t_map <- array(stats::rt(27, df = 24), c(3, 3, 3))
  zvol <- canonicalize(t_map, stat_type = "t", df = 24)

  expect_equal(dim(zvol), c(3, 3, 3))
  ref <- canonicalize_stat(t_map, stat_type = "t", df = 24)$Zeq
  expect_equal(as.numeric(zvol), as.numeric(ref))
})

test_that("canonicalize() converts -log10(p) maps", {
  lp <- array(c(0.5, 1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.4), c(2, 2, 2))
  zvol <- canonicalize(lp, stat_type = "neglog10p", p_side = "one")
  ref  <- canonicalize_stat(lp, stat_type = "neglog10p", p_side = "one")$Zeq
  expect_equal(as.numeric(zvol), as.numeric(ref))
})

test_that("hier_scan no longer accepts stat_type or df", {
  set.seed(3)
  z <- array(rnorm(64), c(4, 4, 4))
  expect_error(
    hier_scan(z, n_perm = 5, kappa = 1.0, min_voxels = 4, stat_type = "t"),
    "stat_type"
  )
  expect_error(
    hier_scan(z, n_perm = 5, kappa = 1.0, min_voxels = 4, df = 10),
    "df"
  )
})

test_that("hier_scan still runs on already-canonicalized Z input", {
  set.seed(4)
  z <- array(rnorm(64), c(4, 4, 4))
  res <- hier_scan(z, n_perm = 10, kappa = 1.0, min_voxels = 4, seed = 1)
  expect_s3_class(res, "neurothresh_result")
})
