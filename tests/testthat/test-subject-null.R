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

test_that("subject-level null_spec drives tfce_fwer", {
  set.seed(3)
  subj <- replicate(6, array(rnorm(6 * 6 * 6), dim = c(6, 6, 6)), simplify = FALSE)
  ns <- null_one_sample(subj, n_perm = 10, seed = 1)

  res <- tfce_fwer(
    ns$stat_vol,
    mask = ns$mask,
    alpha = 0.2,
    null_spec = ns,
    tail = "pos"
  )

  expect_s3_class(res, "neurothresh_result")
  expect_true(is.finite(as.numeric(res$threshold)))
  expect_equal(length(res$null_max), 10)
})

test_that("subject-level null_fun can drive cluster_fdr permutation path", {
  set.seed(4)
  subj <- replicate(6, array(rnorm(6 * 6 * 6), dim = c(6, 6, 6)), simplify = FALSE)
  out <- make_null_fun_subject_signflip(subj, seed = 1)

  res <- cluster_fdr(
    out$z_vol,
    mask = out$mask,
    cluster_thresh = 2.0,
    q = 0.2,
    p_method = "perm",
    n_perm = 10,
    null_fun = out$null_fun,
    tail = "pos"
  )

  expect_true(is.data.frame(res$table))
  expect_equal(res$n_clusters, sum(res$table$significant))
})
