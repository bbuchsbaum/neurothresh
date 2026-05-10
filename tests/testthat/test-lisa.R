test_that("null_from_maps wraps list, 4D array, and matrix null maps", {
  set.seed(1)
  z <- array(rnorm(27), c(3, 3, 3))
  mask <- array(TRUE, c(3, 3, 3))
  mask[1, 1, 1] <- FALSE
  null_list <- replicate(4, array(rnorm(27), c(3, 3, 3)), simplify = FALSE)

  ns_list <- null_from_maps(z, null_list, mask = mask)
  expect_true(is_null_spec(ns_list))
  expect_identical(ns_list$kind, "external_maps")
  expect_equal(ns_list$n_perm, 4L)
  expect_length(ns_list$draw(1), sum(mask))
  expect_identical(ns_list$draw(2), ns_list$draw(2))

  arr4 <- array(NA_real_, c(3, 3, 3, 4))
  for (i in seq_len(4)) arr4[, , , i] <- null_list[[i]]
  ns_4d <- null_from_maps(z, arr4, mask = mask)
  expect_equal(ns_list$draw(3), ns_4d$draw(3))

  mat <- t(vapply(null_list, function(x) as.numeric(x)[mask], numeric(sum(mask))))
  ns_mat <- null_from_maps(z, mat, mask = mask)
  expect_equal(ns_list$draw(4), ns_mat$draw(4))
  expect_error(ns_mat$draw(5), "out of range")
})

test_that("lisa_bilateral validates filter settings", {
  f <- lisa_bilateral(radius = 1, sigma_space = 1.5, sigma_range = 2, n_iter = 1)
  expect_true(is_lisa_filter(f))
  expect_identical(f$type, "bilateral")
  expect_error(lisa_bilateral(radius = 0), "positive integer")
  expect_error(lisa_bilateral(sigma_space = 0), "positive")
  expect_error(lisa_bilateral(sigma_range = -1), "positive")
})

test_that("reference LISA filter preserves constant maps inside the mask", {
  z <- array(3, c(3, 3, 3))
  mask <- array(TRUE, c(3, 3, 3))
  f <- lisa_bilateral(radius = 1, sigma_space = 1, sigma_range = 1,
                      n_iter = 1, scale = "none",
                      boundary = "mask_normalized")
  out <- lisa_apply_filter(z, mask = mask, filter = f, backend = "reference")
  expect_equal(as.numeric(out)[mask], rep(3, sum(mask)), tolerance = 1e-12)
})

test_that("reference LISA filter smooths an impulse locally", {
  z <- array(0, c(5, 5, 5))
  z[3, 3, 3] <- 10
  mask <- array(TRUE, c(5, 5, 5))
  f <- lisa_bilateral(radius = 1, sigma_space = 2, sigma_range = 100,
                      n_iter = 1, scale = "none",
                      boundary = "mask_normalized")
  out <- lisa_apply_filter(z, mask = mask, filter = f, backend = "reference")
  expect_lt(out[3, 3, 3], 10)
  expect_gt(out[3, 3, 2], 0)
})

test_that("neuroim2 backend returns a filtered map with input dimensions", {
  skip_if_not_installed("neuroim2")
  z <- array(rnorm(27), c(3, 3, 3))
  mask <- array(TRUE, c(3, 3, 3))
  f <- lisa_bilateral(radius = 1, sigma_space = 1, sigma_range = 1,
                      n_iter = 1, boundary = "neuroim2")
  out <- lisa_apply_filter(z, mask = mask, filter = f, backend = "neuroim2")
  expect_equal(dim(out), dim(z))
  expect_true(all(is.finite(as.numeric(out)[mask])))
})

test_that("lisa_fdr_map returns bounded monotone q-values and threshold", {
  obs <- c(4, 3, 2, 0)
  nulls <- rbind(c(0, 0.5, 1, 1.5), c(0.2, 0.4, 0.8, 1.0))
  fit <- lisa_fdr_map(obs, nulls, tail = "pos", q = 0.5)
  expect_true(all(fit$fdr >= 0 & fit$fdr <= 1))
  expect_true(all(diff(fit$curve$fdr) >= -1e-12))
  expect_true(is.finite(fit$threshold))
  expect_equal(fit$reject, fit$fdr <= 0.5)
})

test_that("lisa consumes null_spec and returns neurothresh_result", {
  set.seed(11)
  z <- array(rnorm(27), c(3, 3, 3))
  z[2, 2, 2] <- 4
  null_maps <- replicate(5, array(rnorm(27), c(3, 3, 3)), simplify = FALSE)
  ns <- null_from_maps(z, null_maps)
  f <- lisa_bilateral(radius = 1, sigma_space = 1, sigma_range = 2,
                      n_iter = 1, boundary = "mask_normalized")
  res <- lisa(stat_vol = z, mask = ns$mask, null_spec = ns,
              filter = f, backend = "reference", tail = "pos", q = 0.5)
  expect_s3_class(res, "neurothresh_result")
  expect_identical(res$method, "lisa")
  expect_equal(dim(res$pvals), dim(z))
  expect_equal(dim(res$reject), dim(z))
  expect_equal(dim(lisa_filtered_stat(res)), dim(z))
  expect_true(all(as.numeric(res$pvals) >= 0 & as.numeric(res$pvals) <= 1))
  expect_null(res$extras$filtered_nulls)
  expect_length(res$extras$null_tail_counts, sum(is.finite(z)))
})

test_that("lisa keeps null_spec positional convenience but canonical n_perm wins", {
  set.seed(15)
  z <- array(rnorm(8), c(2, 2, 2))
  ns <- null_from_maps(z, replicate(4, array(rnorm(8), c(2, 2, 2)), simplify = FALSE))
  f <- lisa_bilateral(radius = 1, sigma_space = 1, sigma_range = 2,
                      n_iter = 1, boundary = "mask_normalized")
  res <- lisa(ns, filter = f, backend = "reference", q = 0.5,
              n_perm = 2, store_filtered_nulls = TRUE)
  expect_equal(res$params$n_perm, 2L)
  expect_equal(dim(res$extras$filtered_nulls), c(2L, 8L))
  expect_error(
    lisa(stat_vol = z, null_spec = ns, n_perm = 5, filter = f,
         backend = "reference"),
    "n_perm may not exceed null_spec\\$n_perm"
  )
})

test_that("lisa can optionally store filtered null maps for diagnostics", {
  set.seed(12)
  z <- array(rnorm(8), c(2, 2, 2))
  ns <- null_from_maps(z, replicate(3, array(rnorm(8), c(2, 2, 2)), simplify = FALSE))
  f <- lisa_bilateral(radius = 1, sigma_space = 1, sigma_range = 2,
                      n_iter = 1, boundary = "mask_normalized")
  res <- lisa(stat_vol = z, mask = ns$mask, null_spec = ns,
              filter = f, backend = "reference", q = 0.5,
              store_filtered_nulls = TRUE)
  expect_equal(dim(res$extras$filtered_nulls), c(3L, 8L))
})

test_that("same null_spec can feed lisa and scale_space_energy", {
  set.seed(12)
  z <- array(rnorm(27), c(3, 3, 3))
  null_maps <- replicate(4, array(rnorm(27), c(3, 3, 3)), simplify = FALSE)
  ns <- null_from_maps(z, null_maps)
  f <- lisa_bilateral(radius = 1, sigma_space = 1, sigma_range = 2,
                      n_iter = 1, boundary = "mask_normalized")
  lisa_fit <- lisa(stat_vol = z, mask = ns$mask, null_spec = ns,
                   filter = f, backend = "reference", tail = "two", q = 0.5)
  ss_fit <- scale_space_energy(
    null_spec = ns,
    scales = scale_space(2, include_global = FALSE),
    statistics = list(signed_coherence(tail = "two")),
    n_perm = 4
  )
  expect_s3_class(lisa_fit, "neurothresh_result")
  expect_s3_class(ss_fit, "scale_space_energy_result")
})

test_that("lisa wrappers delegate through shared null factories", {
  set.seed(13)
  imgs <- replicate(5, array(rnorm(27), c(3, 3, 3)), simplify = FALSE)
  f <- lisa_bilateral(radius = 1, sigma_space = 1, sigma_range = 2,
                      n_iter = 1, boundary = "mask_normalized")
  one <- lisa_onesample(imgs, n_perm = 4, seed = 1, filter = f,
                        backend = "reference", q = 0.5)
  expect_s3_class(one, "neurothresh_result")

  two <- lisa_twosample(imgs, group = c("a", "a", "b", "b", "b"),
                        n_perm = 4, seed = 2, filter = f,
                        backend = "reference", q = 0.5)
  expect_s3_class(two, "neurothresh_result")
})

test_that("lisa_glm delegates through null_glm", {
  set.seed(14)
  imgs <- replicate(6, array(rnorm(27), c(3, 3, 3)), simplify = FALSE)
  dat <- data.frame(group = factor(rep(c("a", "b"), each = 3)))
  f <- lisa_bilateral(radius = 1, sigma_space = 1, sigma_range = 2,
                      n_iter = 1, boundary = "mask_normalized")
  res <- lisa_glm(
    imgs,
    design = ~ group,
    contrast = "groupb",
    data = dat,
    n_perm = 4,
    seed = 3,
    statistic = "z",
    filter = f,
    backend = "reference",
    q = 0.5
  )
  expect_s3_class(res, "neurothresh_result")
  expect_identical(res$params$null_kind, "glm")
})

test_that("unsupported LISA backends and lisa2018 boundary fail clearly", {
  z <- array(rnorm(8), c(2, 2, 2))
  ns <- null_from_maps(z, replicate(2, array(rnorm(8), c(2, 2, 2)), simplify = FALSE))
  expect_error(
    lisa(ns, filter = lisa_bilateral(boundary = "lisa2018"), backend = "reference"),
    "not implemented"
  )
  expect_error(
    lisa(ns, backend = "adjoin_sparse"),
    "planned but not implemented"
  )
})
