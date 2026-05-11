test_that("scale_space includes global, gaussian, and identity endpoints", {
  ss <- scale_space(c(8, 4), include_identity = TRUE, include_global = TRUE)

  expect_s3_class(ss, "scale_space_spec")
  expect_equal(ss$scales$type, c("global", "gaussian", "gaussian", "identity"))
  expect_equal(ss$scales$fwhm_mm, c(Inf, 4, 8, 0))
  expect_equal(ss$scales$sigma_mm[2], 4 / sqrt(8 * log(2)), tolerance = 1e-12)
})

test_that("channel contracts handle sign and energy metamorphisms", {
  z <- c(-2, -1, 0, 3)
  dims <- c(2, 2, 1)
  mask <- rep(TRUE, 4)
  mask_idx <- seq_along(mask)
  template <- array(0, dims)
  id <- data.frame(type = "identity", fwhm_mm = 0, sigma_mm = 0)

  e1 <- neurothresh:::.ss_eval_channel_scale(z, local_energy(), id, dims, mask, mask_idx, template)
  e2 <- neurothresh:::.ss_eval_channel_scale(-z, local_energy(), id, dims, mask, mask_idx, template)
  s1 <- neurothresh:::.ss_eval_channel_scale(z, signed_coherence("pos"), id, dims, mask, mask_idx, template)
  s2 <- neurothresh:::.ss_eval_channel_scale(-z, signed_coherence("pos"), id, dims, mask, mask_idx, template)
  ppos <- neurothresh:::.ss_eval_channel_scale(z, positive_energy("pos"), id, dims, mask, mask_idx, template)
  pneg <- neurothresh:::.ss_eval_channel_scale(z, positive_energy("neg"), id, dims, mask, mask_idx, template)
  lp <- neurothresh:::.ss_eval_channel_scale(z, lowpass_energy("local"), id, dims, mask, mask_idx, template)
  glp <- neurothresh:::.ss_eval_channel_scale(z, global_lowpass_energy(), id, dims, mask, mask_idx, template)

  expect_equal(e1, e2)
  expect_equal(s1, -s2)
  expect_equal(ppos, pmax(z, 0)^2)
  expect_equal(pneg, pmax(-z, 0)^2)
  expect_equal(lp, z^2)
  expect_equal(glp, rep(sum(z^2), length(z)))
})

test_that("scale_space_energy returns plus-one corrected maps and summaries", {
  set.seed(1)
  imgs <- replicate(6, array(rnorm(27), c(3, 3, 3)), simplify = FALSE)
  nf <- make_null_fun_onesample(imgs, n_perm = 8, seed = 2)

  fit <- scale_space_energy(
    nf$stat_vol,
    mask = nf$mask,
    scales = scale_space(2, include_global = TRUE, include_identity = TRUE),
    statistics = list(local_energy(), signed_coherence("two")),
    null_fun = nf$null_fun,
    n_perm = 8,
    alpha = 0.2
  )

  expect_s3_class(fit, "scale_space_energy_result")
  expect_equal(length(fit$max_null), 8)
  expect_true(all(vapply(fit$p_corr, function(p) all(p >= 1 / 9 & p <= 1), logical(1))))
  expect_true(all(c("best", "finest_detected", "stable_range") %in% names(fit$selected)))
  expect_equal(nrow(summary(fit)), 6)
  expect_true("effective_basis" %in% names(summary(fit)))
  expect_true(all(c("channel_role", "frequency_readout", "curve_use", "scale_label") %in% names(summary(fit))))
  expect_equal(
    summary(fit)$effective_basis[summary(fit)$scale_type == "global"][1],
    1
  )
})

test_that("studentization moments treat observed and null maps symmetrically", {
  set.seed(13)
  dims <- c(3L, 3L, 3L)
  n_full <- prod(dims)
  B <- 12L
  maps <- matrix(rnorm((B + 1L) * n_full), nrow = B + 1L)
  stat_vol <- array(maps[1L, ], dims)
  null_fun <- function(b) maps[b + 1L, ]

  fit <- scale_space_energy(
    stat_vol,
    scales = scale_space(numeric(0), include_global = FALSE, include_identity = TRUE),
    statistics = list(local_energy()),
    null_fun = null_fun,
    n_perm = B,
    studentize = TRUE
  )

  raw <- maps^2
  mu <- colMeans(raw)
  sd <- apply(raw, 2L, stats::sd)
  obs_t <- (raw[1L, ] - mu) / sd
  null_t <- sweep(raw[-1L, , drop = FALSE], 2L, mu, "-")
  null_t <- sweep(null_t, 2L, sd, "/")
  max_null <- apply(null_t, 1L, max)
  expected_p <- neurothresh:::fwer_p_from_maxnull(obs_t, max_null)

  expect_equal(fit$observed$statistic[[1L]], obs_t, tolerance = 1e-12)
  expect_equal(fit$p_corr[[1L]], expected_p, tolerance = 1e-12)
})

test_that("raw-image wrappers match manual low-level calls", {
  set.seed(3)
  imgs <- replicate(5, array(rnorm(8), c(2, 2, 2)), simplify = FALSE)
  scales <- scale_space(2, include_global = FALSE, include_identity = TRUE)
  stats <- list(local_energy())

  nf <- make_null_fun_onesample(imgs, n_perm = 5, seed = 9)
  manual <- scale_space_energy(nf$stat_vol, mask = nf$mask, scales = scales,
                               statistics = stats, null_fun = nf$null_fun,
                               n_perm = 5)
  wrapped <- scale_space_onesample(imgs, scales = scales, statistics = stats,
                                   n_perm = 5, seed = 9)

  expect_equal(wrapped$summary$min_p, manual$summary$min_p, tolerance = 1e-12)
})

test_that("GLM Freedman-Lane null generator produces deterministic contrast maps", {
  set.seed(4)
  n <- 8
  group <- factor(rep(c("control", "case"), each = 4), levels = c("control", "case"))
  age <- seq_len(n)
  imgs <- replicate(n, array(rnorm(8), c(2, 2, 2)), simplify = FALSE)

  nf <- make_null_fun_glm(
    imgs,
    design = ~ group + age,
    data = data.frame(group = group, age = age),
    contrast = "groupcase",
    permutation = permutation_glm(n_perm = 6, seed = 1)
  )

  expect_equal(length(nf$null_fun(2)), sum(nf$mask))
  expect_equal(nf$null_fun(2), nf$null_fun(2), tolerance = 0)
  expect_s3_class(permutation_glm(n_perm = 2), "neurothresh_permutation_glm")
})

test_that("graph_heat matches dense heat kernel on tiny graphs", {
  dims <- c(3L, 2L, 1L)
  mask <- array(TRUE, dims)
  mask_idx <- which(mask)
  template <- array(0, dims)
  z <- c(-1.5, 0.25, 2, -0.5, 1.25, 0.75)
  fwhm <- 1.8
  scale_row <- data.frame(
    type = "gaussian",
    fwhm_mm = fwhm,
    sigma_mm = fwhm / sqrt(8 * log(2))
  )
  cache <- neurothresh:::.ss_build_graph_cache(mask, dims, template)

  got <- neurothresh:::.ss_apply_scale(
    z, scale_row, dims, mask, mask_idx, template,
    backend = "graph_heat", graph_cache = cache
  )

  heat_time <- (scale_row$sigma_mm^2) / 2
  ref <- as.numeric(Matrix::expm(-heat_time * cache$laplacian) %*% z)
  expect_equal(got, ref, tolerance = 1e-8)
})

test_that("graph_heat preserves constants and disconnected mask components", {
  dims <- c(3L, 1L, 1L)
  mask <- array(c(TRUE, FALSE, TRUE), dims)
  mask_idx <- which(mask)
  template <- array(0, dims)
  scale_row <- data.frame(
    type = "gaussian",
    fwhm_mm = 4,
    sigma_mm = 4 / sqrt(8 * log(2))
  )
  cache <- neurothresh:::.ss_build_graph_cache(mask, dims, template)

  expect_equal(cache$degree, c(0, 0))
  expect_equal(
    neurothresh:::.ss_apply_scale(
      c(10, -3), scale_row, dims, mask, mask_idx, template,
      backend = "graph_heat", graph_cache = cache
    ),
    c(10, -3)
  )

  full_mask <- array(TRUE, c(4L, 1L, 1L))
  full_idx <- which(full_mask)
  full_template <- array(0, dim(full_mask))
  full_cache <- neurothresh:::.ss_build_graph_cache(full_mask, dim(full_mask), full_template)
  const <- rep(7, length(full_idx))
  expect_equal(
    neurothresh:::.ss_apply_scale(
      const, scale_row, dim(full_mask), full_mask, full_idx, full_template,
      backend = "graph_heat", graph_cache = full_cache
    ),
    const,
    tolerance = 1e-8
  )
})

test_that("scale_space_energy accepts calibrated graph_heat backend", {
  set.seed(5)
  imgs <- replicate(5, array(rnorm(27), c(3, 3, 3)), simplify = FALSE)
  nf <- make_null_fun_onesample(imgs, n_perm = 4, seed = 5)

  fit <- scale_space_energy(
    nf$stat_vol,
    mask = nf$mask,
    scales = scale_space(2, include_global = FALSE, include_identity = TRUE),
    statistics = list(local_energy()),
    null_fun = nf$null_fun,
    n_perm = 4,
    backend = "graph_heat"
  )

  expect_s3_class(fit, "scale_space_energy_result")
  expect_equal(fit$params$backend, "graph_heat")
  expect_true(all(is.finite(fit$max_null)))
  expect_true(all(vapply(fit$p_corr, function(p) all(p >= 1 / 5 & p <= 1), logical(1))))
})

test_that("grid_gaussian backend matches a hand-built mask-normalized Gaussian", {
  # Differential test on a tiny grid: the exposed grid_gaussian path must
  # reproduce a separately constructed mask-normalized Gaussian convolution
  # at the requested FWHM. Reference is computed from the same C++ smoother
  # used internally; this checks that the dispatch from
  # scale_space_energy() -> .ss_apply_scale() preserves the spec.
  set.seed(11)
  dims <- c(5L, 5L, 5L)
  z <- array(rnorm(prod(dims)), dims)
  mask <- array(TRUE, dims)
  mask[1, , ] <- FALSE
  mask_idx <- which(mask)
  z_mask <- as.numeric(z)[mask_idx]

  fwhm <- 2
  sigma_mm <- fwhm / sqrt(8 * log(2))
  scale_row <- data.frame(type = "gaussian", fwhm_mm = fwhm, sigma_mm = sigma_mm)

  got <- neurothresh:::.ss_apply_scale(
    z_mask, scale_row, dims, mask, mask_idx, template = z,
    backend = "grid_gaussian", graph_cache = NULL
  )

  full <- numeric(prod(dims))
  full[mask_idx] <- z_mask
  sm <- gaussian_smooth3d_cpp(full, as.integer(dims), sigma_xyz = rep(sigma_mm, 3))
  w_full <- numeric(prod(dims))
  w_full[mask_idx] <- 1
  ww <- gaussian_smooth3d_cpp(w_full, as.integer(dims), sigma_xyz = rep(sigma_mm, 3))
  ref <- sm[mask_idx] / pmax(ww[mask_idx], .Machine$double.eps)
  ref[!is.finite(ref)] <- 0

  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("signed_coherence('two') is |S_rho z|, not S_rho |z|", {
  # Two adjacent opposite-sign blobs should cancel under smoothing
  # (correct two-sided coherence is |S_rho z|), while the buggy
  # S_rho |z| construction reinforces magnitudes. The discriminating
  # check is the relative magnitude: with a coarse-enough FWHM, |S_rho z|
  # is much smaller than S_rho |z| for opposite-sign blobs and matches
  # it for same-sign blobs.
  dims <- c(11L, 1L, 1L)
  z <- array(0, dims)
  z[5, 1, 1] <- 5
  z[6, 1, 1] <- -5
  mask <- array(TRUE, dims)
  mask_idx <- seq_len(prod(dims))
  fwhm <- 4
  sigma_mm <- fwhm / sqrt(8 * log(2))
  scale_row <- data.frame(type = "gaussian", fwhm_mm = fwhm, sigma_mm = sigma_mm)

  two <- neurothresh:::.ss_eval_channel_scale(
    as.numeric(z), signed_coherence("two"), scale_row,
    dims, mask, mask_idx, template = z,
    backend = "grid_gaussian", graph_cache = NULL
  )
  alt <- neurothresh:::.ss_apply_scale(
    abs(as.numeric(z)), scale_row, dims, mask, mask_idx, template = z,
    backend = "grid_gaussian", graph_cache = NULL
  )

  # Opposite-sign cancellation: signed coherence is far below the
  # magnitude-first alternative.
  expect_lt(max(two) * 3, max(alt))
  expect_lt(max(two), max(abs(z)) / 2)

  # Same-sign reinforcement: signed coherence and magnitude-first agree
  # exactly because |z| == z when all signs agree.
  z2 <- array(0, dims)
  z2[5, 1, 1] <- 5
  z2[6, 1, 1] <- 5
  two2 <- neurothresh:::.ss_eval_channel_scale(
    as.numeric(z2), signed_coherence("two"), scale_row,
    dims, mask, mask_idx, template = z2,
    backend = "grid_gaussian", graph_cache = NULL
  )
  alt2 <- neurothresh:::.ss_apply_scale(
    abs(as.numeric(z2)), scale_row, dims, mask, mask_idx, template = z2,
    backend = "grid_gaussian", graph_cache = NULL
  )
  expect_equal(two2, alt2, tolerance = 1e-12)
  expect_gt(max(two2), max(two))
})

test_that("identity scale endpoints reproduce voxelwise channels exactly", {
  z <- c(-2.1, 0.4, 1.7, -0.8, 3.2, 0)
  dims <- c(2L, 3L, 1L)
  mask <- rep(TRUE, prod(dims))
  mask_idx <- seq_along(mask)
  template <- array(0, dims)
  id <- data.frame(type = "identity", fwhm_mm = 0, sigma_mm = 0)

  expect_equal(
    neurothresh:::.ss_eval_channel_scale(z, local_energy(), id, dims, mask, mask_idx, template),
    z^2
  )
  expect_equal(
    neurothresh:::.ss_eval_channel_scale(z, signed_coherence("pos"), id, dims, mask, mask_idx, template),
    z
  )
  expect_equal(
    neurothresh:::.ss_eval_channel_scale(z, signed_coherence("neg"), id, dims, mask, mask_idx, template),
    -z
  )
  expect_equal(
    neurothresh:::.ss_eval_channel_scale(z, signed_coherence("two"), id, dims, mask, mask_idx, template),
    abs(z)
  )
  expect_equal(
    neurothresh:::.ss_eval_channel_scale(z, lowpass_energy("local"), id, dims, mask, mask_idx, template),
    z^2
  )
  expect_equal(
    neurothresh:::.ss_eval_channel_scale(z, global_lowpass_energy(), id, dims, mask, mask_idx, template),
    rep(sum(z^2), length(z))
  )
})

test_that("effective_basis is monotone non-increasing in FWHM", {
  set.seed(23)
  imgs <- replicate(5, array(rnorm(64), c(4, 4, 4)), simplify = FALSE)
  nf <- make_null_fun_onesample(imgs, n_perm = 4, seed = 23)

  fit <- scale_space_energy(
    nf$stat_vol,
    mask = nf$mask,
    scales = scale_space(c(2, 3, 4, 6), include_global = TRUE, include_identity = TRUE),
    statistics = list(local_energy()),
    null_fun = nf$null_fun,
    n_perm = 4
  )

  s <- summary(fit)
  s <- s[order(s$scale_index), ]
  # Within Gaussian scales, FWHM is sorted ascending, so effective basis
  # must be non-increasing.
  g <- s[s$scale_type == "gaussian", ]
  expect_true(all(diff(g$effective_basis) <= 0))
  # Identity must be the largest basis count; global must be 1.
  expect_equal(s$effective_basis[s$scale_type == "global"], 1)
  expect_equal(s$effective_basis[s$scale_type == "identity"], sum(nf$mask))
})

test_that("one-sample sign-flip H0 FWER is not anti-conservative", {
  # Empirical FWER simulation: under H0 (mean-zero subjects), the omnibus
  # corrected p_FWER should obey p ~ uniform on the discrete plus-one grid,
  # so the empirical fraction of p <= alpha should approximate alpha.
  # Use a generous upper bound so the test is reliable across seeds.
  set.seed(101)
  K <- 25
  alpha <- 0.20
  reject <- 0L
  for (k in seq_len(K)) {
    imgs <- replicate(8, array(rnorm(27), c(3, 3, 3)), simplify = FALSE)
    nf <- make_null_fun_onesample(imgs, n_perm = 39, seed = 1000 + k)
    fit <- scale_space_energy(
      nf$stat_vol, mask = nf$mask,
      scales = scale_space(2, include_global = FALSE, include_identity = TRUE),
      statistics = list(local_energy()),
      null_fun = nf$null_fun, n_perm = 39, alpha = alpha
    )
    if (min(fit$summary$min_p) <= alpha) reject <- reject + 1L
  }
  # Expectation under H0 is K * alpha = 5. 5-sigma upper bound under
  # binomial(K, alpha) is K * alpha + 5 * sqrt(K * alpha * (1 - alpha)) ~ 14.
  expect_lte(reject, 14L)
})

test_that("scale_space_energy rejects unknown backends with a clear message", {
  set.seed(7)
  imgs <- replicate(4, array(rnorm(8), c(2, 2, 2)), simplify = FALSE)
  nf <- make_null_fun_onesample(imgs, n_perm = 3, seed = 7)
  expect_error(
    scale_space_energy(
      nf$stat_vol, mask = nf$mask,
      scales = scale_space(2, include_global = FALSE, include_identity = TRUE),
      statistics = list(local_energy()),
      null_fun = nf$null_fun, n_perm = 3,
      backend = "spectral"
    )
  )
})

test_that("make_null_fun_glm refuses multi-column contrasts", {
  set.seed(31)
  n <- 8
  group <- factor(rep(c("a", "b", "c", "d"), each = 2), levels = c("a", "b", "c", "d"))
  imgs <- replicate(n, array(rnorm(8), c(2, 2, 2)), simplify = FALSE)
  expect_error(
    make_null_fun_glm(
      imgs,
      design = ~ group,
      data = data.frame(group = group),
      contrast = c(0, 1, -1, 0),
      permutation = permutation_glm(n_perm = 3, seed = 1)
    ),
    "single-coefficient"
  )
})

test_that("matrix input without template warns about voxel spacing", {
  mask <- array(TRUE, c(2, 2, 2))
  X <- matrix(rnorm(2 * 8), nrow = 2, ncol = 8)
  expect_warning(
    as_subject_matrix(X, mask = mask),
    "voxel spacing"
  )
})

test_that("bandpass energy is explicit and finite on gaussian scales", {
  z <- array(rnorm(27), c(3, 3, 3))
  null_fun <- make_null_fun_voxel_signflip(as.numeric(z))

  fit <- scale_space_energy(
    z,
    scales = scale_space(2, include_global = TRUE, include_identity = TRUE),
    statistics = list(bandpass_energy()),
    null_fun = null_fun,
    n_perm = 4,
    studentize = FALSE
  )

  expect_true("bandpass_energy" %in% fit$channels)
  expect_true(all(is.finite(unlist(fit$observed$statistic))))
  expect_equal(fit$observed$raw[[1]], rep(0, length(fit$observed$raw[[1]])))
  expect_equal(fit$observed$raw[[3]], rep(0, length(fit$observed$raw[[3]])))
})

test_that("spectral energy channels expose frequency-readout diagnostics", {
  set.seed(41)
  dims <- c(5L, 5L, 1L)
  z <- array(rnorm(prod(dims)), dims)
  null <- null_voxel_signflip(z, n_perm = 5, seed = 41)

  fit <- scale_space_energy(
    z,
    scales = scale_space(c(2, 4), include_global = TRUE, include_identity = TRUE),
    statistics = list(local_energy(), global_lowpass_energy(), global_band_energy()),
    null_fun = null$draw,
    n_perm = 5,
    studentize = FALSE
  )

  expect_true(all(c("global_lowpass_energy", "global_band_energy") %in% fit$channels))
  expect_true(all(fit$summary$frequency_readout[fit$summary$channel %in% c("global_lowpass_energy", "global_band_energy")]))
  expect_false(any(fit$summary$frequency_readout[fit$summary$channel == "local_energy"]))

  diag <- scale_space_diagnostics(fit)
  expect_s3_class(diag, "data.frame")
  expect_true("interpretation" %in% names(diag))
  expect_true(diag$frequency_readout[diag$channel == "global_band_energy"])

  curve <- plot_resolution_curve(fit)
  expect_true("resolution_order" %in% names(curve))
  expect_equal(curve$scale_label[curve$scale_type == "global"][1], "global")
  expect_equal(curve$scale_label[curve$scale_type == "identity"][1], "voxel")

  # Global band energy has no meaningful global or identity endpoint.
  gb <- fit$observed$raw[fit$family$channel == "global_band_energy"]
  expect_equal(gb[[1]], rep(0, prod(dims)))
  expect_equal(gb[[length(gb)]], rep(0, prod(dims)))
})

test_that("energy_mass_map avoids spatially constant global endpoint", {
  dims <- c(5L, 5L, 1L)
  z <- array(0, dims)
  z[2:4, 2:4, 1] <- 3
  null <- null_voxel_signflip(z, n_perm = 5, seed = 42)

  fit <- scale_space_energy(
    z,
    scales = scale_space(2, include_global = TRUE, include_identity = TRUE),
    statistics = list(global_lowpass_energy(), local_energy()),
    null_fun = null$draw,
    n_perm = 5,
    studentize = FALSE
  )

  expect_error(
    energy_mass_map(fit, scale_index = fit$summary$scale_index[fit$summary$scale_type == "global"][1]),
    "global endpoint"
  )

  mass <- as.numeric(energy_mass_map(fit))
  mass <- mass[is.finite(mass)]
  expect_equal(sum(mass), 1, tolerance = 1e-12)
  expect_gt(diff(range(mass)), 0)
})
