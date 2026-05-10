#' LISA local adaptive inference
#'
#' Applies post-statistic bilateral filtering to an observed statistic map and
#' matched null maps, then estimates voxelwise FDR in filtered-statistic space.
#' This is a local-adaptive companion to scale-space energy inference.
#'
#' @param stat_vol Optional 3D Z-equivalent statistic map. Defaults to
#'   `null_spec$stat_vol`. For convenience, `stat_vol` may also be a
#'   [`null_spec`] object, in which case `null_spec` is inferred.
#' @param mask Optional mask override. Defaults to `null_spec$mask`.
#' @param q FDR level.
#' @param tail One of `"pos"`, `"neg"`, or `"two"`.
#' @param null_spec A [`null_spec`] object. Its `stat_vol`, `mask`, `draw`, and
#'   `n_perm` slots define the default observed map and null maps.
#' @param n_perm Optional number of null draws to use. Defaults to all draws in
#'   `null_spec`; may not exceed `null_spec$n_perm`.
#' @param seed Optional random seed set before consuming null draws. Most
#'   `null_spec` factories are deterministic by draw index, but this keeps LISA
#'   aligned with other permutation-facing neurothresh methods.
#' @param filter LISA filter specification from [lisa_bilateral()].
#' @param p0 Conservative null proportion for the empirical FDR estimator.
#' @param monotone Whether to enforce monotone q-values over evidence strength.
#' @param scale_n_null Number of early null maps used for fixed range-scale
#'   estimation when `filter$scale = "sd_first_30_nulls"`.
#' @param backend Filtering backend.
#' @param store_filtered_nulls Whether to store the full permutation-by-voxel
#'   filtered null matrix in `extras`. Defaults to `FALSE`; inference uses
#'   streaming tail counts either way.
#'
#' @return A `neurothresh_result` with LISA q-values in `pvals` and the
#'   filtered statistic map in `extras$filtered_stat`.
#' @export
lisa <- function(stat_vol = NULL,
                 mask = NULL,
                 q = 0.05,
                 tail = c("pos", "neg", "two"),
                 null_spec = NULL,
                 n_perm = NULL,
                 seed = NULL,
                 filter = lisa_bilateral(),
                 p0 = 1,
                 monotone = TRUE,
                 scale_n_null = 30L,
                 backend = c("neuroim2", "reference", "adjoin_sparse"),
                 store_filtered_nulls = FALSE) {
  call_capture <- match.call()
  tail <- match.arg(tail)
  backend <- match.arg(backend)
  if (is_null_spec(stat_vol) && is.null(null_spec)) {
    null_spec <- stat_vol
    stat_vol <- NULL
  }
  if (!is_null_spec(null_spec)) stop("null_spec must be a neurothresh_null_spec (see ?null_spec)")
  validate_null_spec(null_spec)
  if (!is_lisa_filter(filter)) stop("filter must be created by lisa_bilateral()")
  if (!is.null(seed)) set.seed(seed)

  if (is.null(n_perm)) {
    n_perm <- null_spec$n_perm
  } else {
    if (!is.numeric(n_perm) || length(n_perm) != 1L || n_perm < 1L) {
      stop("n_perm must be a positive integer")
    }
    n_perm <- as.integer(n_perm)
    if (n_perm > null_spec$n_perm) {
      stop("n_perm may not exceed null_spec$n_perm")
    }
  }

  if (is.null(stat_vol)) stat_vol <- null_spec$stat_vol
  if (is.null(mask)) mask <- null_spec$mask

  dims <- dim(stat_vol)
  if (is.null(dims) || length(dims) != 3L) stop("stat_vol must be a 3D volume-like object")
  n_full <- prod(dims)
  if (is.null(mask)) {
    mask <- array(is.finite(as.numeric(stat_vol)), dim = dims)
  } else {
    if (length(mask) != n_full) stop("mask must match stat_vol dimensions")
    mask <- array(as.logical(mask), dim = dims)
  }
  mask_idx <- which(mask)
  if (!length(mask_idx)) stop("mask is empty")
  if (any(!is.finite(as.numeric(stat_vol)[mask_idx]))) {
    stop("stat_vol contains non-finite values in mask")
  }

  range_scale <- switch(
    filter$scale,
    sd_first_30_nulls = .lisa_estimate_range_scale(null_spec, mask_idx, n_full, n = scale_n_null),
    per_map = NULL,
    none = 1
  )

  filtered_obs <- .lisa_apply_filter(stat_vol, mask, filter, backend, range_scale = range_scale)
  obs_vec <- as.numeric(filtered_obs)[mask_idx]

  counter <- .lisa_fdr_counter_init(obs_vec, tail = tail)
  null_mat <- if (isTRUE(store_filtered_nulls)) {
    matrix(NA_real_, nrow = n_perm, ncol = length(mask_idx))
  } else {
    NULL
  }
  for (b in seq_len(n_perm)) {
    draw <- .coerce_null_draw(null_spec$draw(b), n_full = n_full, mask_idx = mask_idx)
    draw_vol <- .array_to_vol(draw, dims, stat_vol)
    filt_b <- .lisa_apply_filter(draw_vol, mask, filter, backend, range_scale = range_scale)
    filt_vec <- as.numeric(filt_b)[mask_idx]
    counter <- .lisa_fdr_counter_update(counter, filt_vec, tail = tail)
    if (isTRUE(store_filtered_nulls)) null_mat[b, ] <- filt_vec
  }

  fdr <- .lisa_fdr_counter_finalize(counter, q = q, p0 = p0, monotone = monotone)

  p_full <- rep(1, n_full)
  p_full[mask_idx] <- fdr$fdr
  reject_full <- rep(FALSE, n_full)
  reject_full[mask_idx] <- fdr$reject

  new_neurothresh_result(
    method = "lisa",
    stat_vol = stat_vol,
    mask = mask,
    threshold = fdr$threshold,
    reject = .array_to_vol(reject_full, dims, stat_vol),
    pvals = .array_to_vol(p_full, dims, stat_vol),
    clusters = NULL,
    null_max = NULL,
    params = list(
      q = q,
      tail = tail,
      n_perm = n_perm,
      seed = seed,
      p0 = p0,
      backend = backend,
      filter = filter,
      null_kind = null_spec$kind,
      range_scale = range_scale,
      store_filtered_nulls = isTRUE(store_filtered_nulls)
    ),
    call = call_capture,
    extras = Filter(Negate(is.null), list(
      filtered_stat = filtered_obs,
      fdr_curve = fdr$curve,
      null_tail_counts = counter$null_counts,
      filtered_nulls = null_mat
    ))
  )
}

#' Extract the LISA filtered statistic map
#'
#' @param fit A `neurothresh_result` returned by [lisa()].
#' @return Filtered statistic map.
#' @export
lisa_filtered_stat <- function(fit) {
  if (!inherits(fit, "neurothresh_result") || !identical(fit$method, "lisa")) {
    stop("fit must be a LISA neurothresh_result")
  }
  fit$extras$filtered_stat
}

#' LISA one-sample wrapper
#'
#' @inheritParams null_one_sample
#' @param ... Additional arguments passed to [lisa()].
#' @export
lisa_onesample <- function(images, mask = NULL, n_perm = 5000, seed = NULL,
                           z_equiv = TRUE, template = NULL, ...) {
  ns <- null_one_sample(images = images, mask = mask, n_perm = n_perm,
                        seed = seed, z_equiv = z_equiv, template = template)
  lisa(stat_vol = ns$stat_vol, mask = ns$mask, null_spec = ns,
       n_perm = n_perm, seed = seed, ...)
}

#' LISA two-sample wrapper
#'
#' @inheritParams null_two_sample
#' @param ... Additional arguments passed to [lisa()].
#' @export
lisa_twosample <- function(images, group, mask = NULL, n_perm = 5000,
                           seed = NULL, variance = c("welch", "pooled"),
                           z_equiv = TRUE, template = NULL,
                           levels = NULL, contrast = NULL, ...) {
  variance <- match.arg(variance)
  ns <- null_two_sample(images = images, group = group, mask = mask,
                        n_perm = n_perm, seed = seed, variance = variance,
                        z_equiv = z_equiv, template = template,
                        levels = levels, contrast = contrast)
  lisa(stat_vol = ns$stat_vol, mask = ns$mask, null_spec = ns,
       n_perm = n_perm, seed = seed, ...)
}

#' LISA GLM wrapper
#'
#' @inheritParams null_glm
#' @param ... Additional arguments passed to [lisa()].
#' @export
lisa_glm <- function(images, design, contrast, mask = NULL, data = NULL,
                     n_perm = 5000, seed = NULL,
                     statistic = c("t", "z"), template = NULL, ...) {
  statistic <- match.arg(statistic)
  ns <- null_glm(images = images, design = design, contrast = contrast,
                 mask = mask, data = data, n_perm = n_perm, seed = seed,
                 statistic = statistic, template = template)
  lisa(stat_vol = ns$stat_vol, mask = ns$mask, null_spec = ns,
       n_perm = n_perm, seed = seed, ...)
}
