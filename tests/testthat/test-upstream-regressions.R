test_that("cluster permutation tails transform observed and null maps alike (#1)", {
  d <- c(3, 3, 3)
  mask <- array(TRUE, d)
  for (s in c(1, -1)) {
    obs <- array(s * 4, d)
    for (masked in c(FALSE, TRUE)) {
      mask[1] <- !masked
      ns <- null_from_function(function(b) rep(s * 4, sum(mask)), obs, 39,
                               mask = mask)
      for (tail in c(if (s == 1) "pos" else "neg", "two")) {
        fit <- cluster_fdr(obs, mask, tail = tail, p_method = "perm",
                           null_spec = ns, cluster_thresh = 3)
        expect_equal(fit$clusters$p_unc, 1)
        expect_false(any(fit$reject))
      }
    }
  }
})

test_that("two-sided cluster FDR paints significant clusters from both tails", {
  z <- array(0, c(5, 2, 2))
  z[1:2, , ] <- 4
  z[4:5, , ] <- -4
  ns <- null_from_function(function(b) rep(0, length(z)), z, 39)
  fit <- cluster_fdr(z, tail = "two", p_method = "perm", null_spec = ns)
  expect_true(all(fit$clusters$significant))
  expect_equal(fit$reject, abs(z) >= 3)
  expect_equal(fit$extras$reject_pos, z >= 3)
  expect_equal(fit$extras$reject_neg, z <= -3)
})

test_that("two-sided RFT peak spends half alpha per tail (#2)", {
  z <- array(0, c(1, 1, 1))
  fit <- rft_peak_fwer(z, alpha = .05, tail = "two", smoothness = 1e6)
  expect_equal(fit$threshold, qnorm(.975), tolerance = 1e-10)
  expect_equal(2 * pnorm(fit$threshold, lower.tail = FALSE), .05,
               tolerance = 1e-10)
})

test_that("two-sided RFT cluster p-values describe the combined family (#2)", {
  z <- array(0, c(8, 4, 4))
  z[1:3, , ] <- 4
  z[6:8, , ] <- -4
  for (a in c(.05, .5)) {
    pos <- rft_cluster_fwer(z, tail = "pos", alpha = a / 2, smoothness = 6)
    neg <- rft_cluster_fwer(z, tail = "neg", alpha = a / 2, smoothness = 6)
    two <- rft_cluster_fwer(z, tail = "two", alpha = a, smoothness = 6)
    expect_equal(two$clusters$p_fwer,
                 pmin(1, 2 * c(pos$clusters$p_fwer, neg$clusters$p_fwer)))
    expect_equal(two$clusters$significant, two$clusters$p_fwer <= a)
    expect_equal(two$reject, pos$reject | neg$reject)
  }
  boundary <- array(0, c(20, 20, 20))
  boundary[1:4, 1:4, 1:3] <- 4
  one <- rft_cluster_fwer(boundary, tail = "pos", smoothness = 8)
  two <- rft_cluster_fwer(boundary, tail = "two", smoothness = 8)
  expect_true(one$clusters$significant)
  expect_false(two$clusters$significant)
  expect_gt(sum(one$reject), 0)
  expect_false(any(two$reject))
})

test_that("hier_scan refuses unsupported supplied nulls before drawing (#3)", {
  z <- array(4, c(4, 4, 4))
  calls <- 0L
  ns <- null_from_function(function(b) { calls <<- calls + 1L; rep(1e4, 64) },
                           z, 39)
  for (method in c("stepdown", "hierarchical")) {
    for (parcels in list(NULL, array(rep(1:2, each = 32), dim(z)))) {
      expect_error(hier_scan(z, null_spec = ns, n_perm = 39,
                             method = method, parcels = parcels),
                   "null_spec.*not supported")
    }
  }
  expect_identical(calls, 0L)
})

test_that("signed t-to-Z conversion preserves extreme tails and zero (#4)", {
  t <- c(0, .1, 1, 10, 100, 1e4, 1e100)
  for (two in c(FALSE, TRUE)) {
    pos <- neurothresh:::.ss_t_to_z(t, 30, two_sided = two)
    neg <- neurothresh:::.ss_t_to_z(-t, 30, two_sided = two)
    expect_equal(neg, -pos, tolerance = 1e-12)
    expect_identical(pos[1], 0)
    expect_true(all(is.finite(pos)))
    expect_true(all(diff(pos) > 0))
    # Equality of log survival probabilities checks the quantile conversion.
    expect_equal(pnorm(pos, lower.tail = FALSE, log.p = TRUE),
                 pt(t, 30, lower.tail = FALSE, log.p = TRUE), tolerance = 1e-10)
  }
  x <- lapply(1:31, function(i) array(100 + i / 31, c(1, 1, 1)))
  p <- null_one_sample(x, n_perm = 3, seed = 1)
  n <- null_one_sample(lapply(x, function(v) -v), n_perm = 3, seed = 1)
  expect_equal(as.numeric(n$stat_vol), -as.numeric(p$stat_vol), tolerance = 1e-12)
  for (b in 1:3) expect_equal(n$draw(b), -p$draw(b), tolerance = 1e-12)
})

test_that("null validation rejects any nonfinite value in the first draw (#5)", {
  z <- array(0, c(2, 2, 2))
  for (value in c(NA_real_, NaN, Inf, -Inf)) {
    ns <- null_from_function(function(b) c(rep(0, 7), value), z, 2)
    expect_error(validate_null_spec(ns), "draw\\(1\\).*non-finite")
  }
  expect_silent(validate_null_spec(null_from_function(function(b) rep(0, 8), z, 2)))
})

test_that("LISA batch and streaming match brute-force tail counts with ties (#6)", {
  obs <- c(-3, -1, 0, 0, 1, 1, 3)
  nulls <- rbind(obs, rev(obs), c(-4, -1, -1, 0, 0, 1, 4))
  for (tail in c("pos", "neg", "two")) {
    transform <- switch(tail, pos = identity, neg = function(x) -x, two = abs)
    score <- transform(obs)
    ns <- transform(as.numeric(nulls))
    ord <- order(score, decreasing = TRUE)
    raw <- vapply(score[ord], function(s) {
      min(1, .8 * mean(ns >= s) / mean(score >= s))
    }, numeric(1))
    for (monotone in c(FALSE, TRUE)) {
      expected <- if (monotone) rev(cummin(rev(raw))) else raw
      fit <- lisa_fdr_map(obs, nulls, tail = tail, p0 = .8, q = .5,
                          monotone = monotone)
      counter <- neurothresh:::.lisa_fdr_counter_init(obs, tail)
      for (b in seq_len(nrow(nulls))) {
        counter <- neurothresh:::.lisa_fdr_counter_update(counter, nulls[b, ], tail)
      }
      stream <- neurothresh:::.lisa_fdr_counter_finalize(counter, .5, .8, monotone)
      expect_equal(fit$curve$raw_fdr, raw, tolerance = 1e-14)
      expect_equal(fit$fdr[ord], expected, tolerance = 1e-14)
      expect_equal(stream, fit, tolerance = 1e-14)
      expect_equal(fit$reject, fit$fdr <= .5)
    }
  }
})

test_that("LISA rank counts include ties at endpoints and outside the data", {
  x <- c(-1, 0, 0, 1, 1, 2)
  thresholds <- c(3, 2, 1, 0, -1, -2, .5)
  expect_equal(neurothresh:::.lisa_count_ge(x, thresholds),
               vapply(thresholds, function(s) sum(x >= s), integer(1)))
  expect_equal(neurothresh:::.lisa_count_ge(numeric(0), thresholds),
               rep(0, length(thresholds)))
})

test_that("public LISA streaming agrees with stored filtered-map batch inference", {
  set.seed(6)
  z <- array(rnorm(27), c(3, 3, 3))
  mask <- array(TRUE, dim(z))
  mask[1] <- FALSE
  ns <- null_from_maps(z, replicate(3, array(rnorm(27), dim(z)), simplify = FALSE),
                       mask = mask)
  fit <- lisa(null_spec = ns, tail = "two", store_filtered_nulls = TRUE,
               filter = lisa_bilateral(radius = 1, scale = "none"))
  batch <- lisa_fdr_map(fit$extras$filtered_stat, fit$extras$filtered_nulls,
                        mask = mask, tail = "two")
  expect_equal(as.numeric(fit$pvals)[mask], batch$fdr, tolerance = 1e-14)
  expect_equal(as.logical(fit$reject)[mask], batch$reject)
})

test_that("RFT cluster probabilities are invariant to physical unit scaling (#7)", {
  d <- c(10, 10, 10)
  a <- array(0, d)
  a[2:9, 2:9, 2:9] <- 4
  run <- function(s) {
    vol <- neuroim2::NeuroVol(a, neuroim2::NeuroSpace(d, spacing = c(2, 3, 4) * s))
    fw <- rft_cluster_fwer(vol, smoothness = c(6, 9, 12) * s)
    fd <- cluster_fdr(vol, p_method = "rft", smoothness = c(6, 9, 12) * s)
    expect_equal(fw$clusters$p_unc, fd$clusters$p_unc)
    c(fw$clusters$size_resels, fw$clusters$p_unc, fw$clusters$p_fwer,
      fd$clusters$p_unc)
  }
  ref <- run(1)
  expect_equal(run(10), ref, tolerance = 1e-12)
  expect_equal(run(.1), ref, tolerance = 1e-12)
  expect_equal(run(c(2, 3, 4)), ref, tolerance = 1e-12)
})

test_that("RFT extent tail approaches the independent EC moment approximation", {
  # SPM's spm_P_RF.m uses beta = (Gamma(5/2) / E[K])^(2/3),
  # E[K] = rho0/rho3. The package retains the high-u asymptote of this rule.
  # Compute rho3 via Gaussian density and rho0 via survival probability,
  # independently of the production extent helper.
  u <- c(3, 10, 20)
  rho3 <- (4 * log(2))^(3 / 2) / (2 * pi)^(3 / 2) *
    (u^2 - 1) * dnorm(u)
  mean_extent <- pnorm(u, lower.tail = FALSE) / rho3
  beta_ec <- (gamma(5 / 2) / mean_extent)^(2 / 3)
  beta_asym <- -log(neurothresh:::.rft_cluster_extent_tail(1, u))
  discrepancy <- abs(beta_ec / beta_asym - 1)
  expect_true(all(diff(discrepancy) < 0))
  expect_lt(discrepancy[3], 2e-5)
  # Integrating the extent survival curve must recover the high-u mean.
  mean_numerical <- integrate(function(k) {
    neurothresh:::.rft_cluster_extent_tail(k, 3)
  }, 0, Inf)$value
  mean_asym <- (2 * pi / (4 * log(2)))^(3 / 2) / 3^3
  expect_equal(mean_numerical, mean_asym, tolerance = 1e-6)
})

test_that("octree metadata identifies the null actually consumed (#8)", {
  z <- array(0, c(2, 2, 2))
  calls <- 0L
  draw <- function(b) { calls <<- calls + 1L; rep(0, 8) }
  ns <- null_from_function(draw, z, 2, kind = "synthetic_custom")
  for (scan in list(octree_scan_fwer, octree_scan_stepdown)) {
    calls <- 0L
    fit <- scan(z, null_spec = ns)
    expect_gt(calls, 0L)
    expect_identical(fit$params$null, ns$kind)
    fit <- scan(z, null_fun = draw, n_perm = 2)
    expect_identical(fit$params$null, "custom")
    expect_warning(fit <- scan(z, n_perm = 2), "signflip_voxel")
    expect_identical(fit$params$null, "signflip_voxel")
  }
})
