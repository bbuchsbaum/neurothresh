# Tests for subject-level sign-flip null generator

test_that("make_null_fun_subject_signflip returns correct shapes", {
  set.seed(1)
  subj <- replicate(5, array(rnorm(4 * 4 * 4), dim = c(4, 4, 4)), simplify = FALSE)
  out <- make_null_fun_subject_signflip(subj, seed = 1)
  expect_equal(dim(out$z_vol), c(4, 4, 4))
  expect_true(is.function(out$null_fun))
  z1 <- out$null_fun(1)
  expect_equal(length(z1), sum(out$mask))
})

test_that("subject-level null_fun can drive octree_scan_fwer", {
  set.seed(2)
  subj <- replicate(6, array(rnorm(6 * 6 * 6), dim = c(6, 6, 6)), simplify = FALSE)
  out <- make_null_fun_subject_signflip(subj, seed = 1)
  res <- octree_scan_fwer(out$z_vol, mask = out$mask, n_perm = 10, null_fun = out$null_fun)
  expect_true(is.numeric(res$u))
})

