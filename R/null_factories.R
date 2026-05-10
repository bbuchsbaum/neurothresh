#' Voxel-level sign-flip null
#'
#' Permutes the sign of each voxel independently. This null is anti-conservative
#' in the presence of spatial correlation; prefer subject-level sign-flips
#' ([null_one_sample()] / [null_subject_signflip()]) when subject data are
#' available.
#'
#' @param stat_vol Observed (canonicalized) Z-equivalent statistic map.
#' @param mask Optional logical 3D mask. Defaults to finite voxels in
#'   `stat_vol`.
#' @param n_perm Number of sign-flip draws.
#' @param seed Optional random seed; the resulting `draw(b)` is deterministic
#'   by index.
#' @return A [`null_spec`].
#' @export
null_voxel_signflip <- function(stat_vol, mask = NULL, n_perm = 5000, seed = NULL) {
  if (n_perm < 1L) stop("n_perm must be >= 1")
  if (is.null(mask)) mask <- is.finite(as.numeric(stat_vol))
  mask_arr <- if (is.array(mask) || inherits(mask, "NeuroVol")) {
    array(as.logical(mask), dim = dim(mask))
  } else {
    array(as.logical(mask), dim = dim(stat_vol))
  }
  mask_idx <- which(mask_arr)
  z_obs <- as.numeric(stat_vol)[mask_idx]

  if (!is.null(seed)) set.seed(seed)
  signs <- matrix(sample(c(-1, 1), length(mask_idx) * n_perm, replace = TRUE),
                  nrow = n_perm, ncol = length(mask_idx))

  draw <- function(b) {
    if (b < 1L || b > n_perm) stop("Permutation index out of range")
    as.numeric(signs[b, ]) * z_obs
  }

  new_null_spec(
    stat_vol = stat_vol, mask = mask_arr, draw = draw,
    n_perm = n_perm, df = NA_real_, kind = "voxel_signflip",
    meta = list(n_mask = length(mask_idx))
  )
}

#' One-sample subject sign-flip null
#'
#' Builds an observed one-sample t (or Z-equivalent) map and a deterministic
#' subject-level sign-flip null indexed by `b`. Each subject map is multiplied
#' by a single \eqn{\pm 1} per draw, preserving within-subject spatial
#' covariance.
#'
#' @param images Subject-level images. Accepts any input handled by
#'   [as_subject_matrix()]: list of 3D arrays/`NeuroVol`s, 4D array,
#'   `NeuroVec`, `SparseNeuroVec`, or a subject-by-voxel matrix
#'   (with `template`/`mask`).
#' @param mask Optional logical 3D mask.
#' @param n_perm Number of sign-flip draws.
#' @param seed Optional random seed.
#' @param z_equiv If `TRUE` (default), convert the t statistic to a Z
#'   equivalent.
#' @param template Optional template (only used for matrix input).
#' @return A [`null_spec`].
#' @export
null_one_sample <- function(images, mask = NULL, n_perm = 5000, seed = NULL,
                            z_equiv = TRUE, template = NULL) {
  if (n_perm < 1L) stop("n_perm must be >= 1")
  prep <- as_subject_matrix(images, mask = mask, template = template)
  X <- prep$X
  n <- nrow(X)
  if (n < 2L) stop("at least two subjects required")

  if (!is.null(seed)) set.seed(seed)
  signs <- matrix(sample(c(-1, 1), n * n_perm, replace = TRUE),
                  nrow = n_perm, ncol = n)

  stat_from_signs <- function(s) {
    t <- .ss_one_sample_t(X * as.numeric(s))
    if (isTRUE(z_equiv)) .ss_t_to_z(t, df = n - 1) else t
  }

  z_obs <- stat_from_signs(rep(1, n))
  stat_vol <- .ss_stat_template(prep, z_obs)

  draw <- function(b) {
    if (b < 1L || b > n_perm) stop("Permutation index out of range")
    stat_from_signs(signs[b, ])
  }

  new_null_spec(
    stat_vol = stat_vol, mask = prep$mask, draw = draw,
    n_perm = n_perm, df = n - 1, kind = "one_sample",
    meta = list(prep = prep, signs = signs,
                statistic = if (isTRUE(z_equiv)) "z" else "t")
  )
}

#' Subject sign-flip null (alias of one-sample)
#'
#' The "subject signflip" name describes the mechanism. The factory under the
#' design name is [null_one_sample()]; this is an alias kept for callers who
#' prefer the mechanism-first name.
#'
#' @inheritParams null_one_sample
#' @return A [`null_spec`].
#' @export
null_subject_signflip <- function(images, mask = NULL, n_perm = 5000,
                                  seed = NULL, z_equiv = TRUE,
                                  template = NULL) {
  null_one_sample(images = images, mask = mask, n_perm = n_perm, seed = seed,
                  z_equiv = z_equiv, template = template)
}

#' Two-sample label-permutation null
#'
#' Builds an observed two-group t (or Z-equivalent) map and deterministic
#' label-permutation null draws. Default contrast direction is **level 2 minus
#' level 1**; pin a different direction with `levels` or `contrast`.
#'
#' @param images Subject-level images (see [as_subject_matrix()]).
#' @param group Vector of group labels (factor, character, or integer); coerced
#'   to a 2-level factor.
#' @param mask Optional logical 3D mask.
#' @param n_perm Number of label permutations.
#' @param seed Optional random seed.
#' @param variance `"welch"` (default) or `"pooled"`.
#' @param z_equiv If `TRUE` (default), convert t to Z-equivalent.
#' @param template Optional template (only used for matrix input).
#' @param levels Optional length-2 character giving the group level order
#'   used for the contrast (e.g., `c("control", "case")` → contrast is
#'   `case - control`). Overrides `as.factor(group)` ordering.
#' @param contrast Optional length-2 character of the form
#'   `c("positive", "negative")` specifying that the contrast is
#'   `mean(positive) - mean(negative)`. Equivalent to `levels = rev(contrast)`.
#' @return A [`null_spec`].
#' @export
null_two_sample <- function(images, group, mask = NULL,
                            n_perm = 5000, seed = NULL,
                            variance = c("welch", "pooled"),
                            z_equiv = TRUE, template = NULL,
                            levels = NULL, contrast = NULL) {
  variance <- match.arg(variance)
  if (n_perm < 1L) stop("n_perm must be >= 1")
  prep <- as_subject_matrix(images, mask = mask, template = template)
  X <- prep$X

  if (!is.null(contrast)) {
    if (!is.character(contrast) || length(contrast) != 2L) {
      stop("contrast must be a length-2 character (positive, negative)")
    }
    levels <- rev(contrast)
  }
  if (!is.null(levels)) {
    if (!is.character(levels) || length(levels) != 2L) {
      stop("levels must be a length-2 character")
    }
    group <- factor(as.character(group), levels = levels)
  } else {
    group <- as.factor(group)
  }

  if (length(group) != nrow(X)) stop("length(group) must equal number of images")
  if (anyNA(group)) {
    stop("group contains missing labels or labels not present in levels")
  }
  if (nlevels(group) != 2L) stop("group must have exactly two levels")

  if (!is.null(seed)) set.seed(seed)
  perm_index <- replicate(n_perm, sample.int(nrow(X)), simplify = FALSE)

  stat_from_group <- function(g) {
    t <- .ss_two_sample_t(X, g, variance = variance)
    if (!isTRUE(z_equiv)) return(t)
    df <- .ss_two_sample_df(X, g, variance = variance)
    .ss_t_to_z(t, df = df)
  }

  z_obs <- stat_from_group(group)
  stat_vol <- .ss_stat_template(prep, z_obs)

  draw <- function(b) {
    if (b < 1L || b > n_perm) stop("Permutation index out of range")
    stat_from_group(group[perm_index[[b]]])
  }

  df <- if (variance == "pooled") {
    nrow(X) - 2L
  } else {
    NA_real_
  }

  new_null_spec(
    stat_vol = stat_vol, mask = prep$mask, draw = draw,
    n_perm = n_perm, df = df, kind = "two_sample",
    meta = list(prep = prep, group = group, perm_index = perm_index,
                variance = variance,
                statistic = if (isTRUE(z_equiv)) "z" else "t",
                contrast = paste(rev(levels(group)), collapse = " minus "))
  )
}

#' GLM Freedman-Lane null
#'
#' Wraps [make_null_fun_glm()] (single-coefficient contrasts) and returns a
#' [`null_spec`].
#'
#' @inheritParams make_null_fun_glm
#' @param n_perm Number of permutations (default 5000).
#' @param seed Optional random seed.
#' @return A [`null_spec`].
#' @export
null_glm <- function(images, design, contrast, mask = NULL,
                     data = NULL, n_perm = 5000, seed = NULL,
                     statistic = c("t", "z"), template = NULL) {
  statistic <- match.arg(statistic)
  perm <- permutation_glm(method = "freedman_lane", n_perm = n_perm, seed = seed)
  fit <- make_null_fun_glm(images = images, design = design, contrast = contrast,
                           mask = mask, data = data, permutation = perm,
                           statistic = statistic, template = template)
  draw <- function(b) {
    if (b < 1L || b > n_perm) stop("Permutation index out of range")
    fit$null_fun(b)
  }
  new_null_spec(
    stat_vol = fit$stat_vol, mask = fit$mask, draw = draw,
    n_perm = n_perm, df = NA_real_, kind = "glm",
    meta = list(prep = fit$prep, design = fit$design,
                reduced_design = fit$reduced_design,
                contrast = fit$contrast, perm_index = fit$perm_index,
                statistic = statistic, method = "freedman_lane")
  )
}

#' Monte-Carlo FWHM null (Gaussian-smoothed white noise)
#'
#' Generates null draws by smoothing white noise with a Gaussian kernel that
#' yields the requested FWHM. Returns a [`null_spec`] whose `draw(b)` is
#' length `sum(mask)`.
#'
#' @param template 3D template carrying spatial metadata (`NeuroVol`, array,
#'   or anything with `dim()` of length 3).
#' @param mask Logical 3D mask.
#' @param smoothness Either a [`smoothness_fwhm()`] object, or a numeric FWHM
#'   in voxel units (length 1 or 3) for backwards compatibility.
#' @param n_perm Maximum number of indexed draws supported by the resulting spec.
#' @param seed Optional random seed used to precompute per-index RNG seeds; the
#'   resulting `draw(b)` is deterministic and restores the caller's RNG state.
#' @return A [`null_spec`].
#' @export
null_mc_fwhm <- function(template, mask, smoothness, n_perm = 5000, seed = NULL) {
  if (n_perm < 1L) stop("n_perm must be >= 1")
  dims3 <- dim(template)
  if (is.null(dims3) || length(dims3) != 3L) stop("template must be 3D")
  side <- max(dims3)
  mask_arr <- array(as.logical(mask), dim = dims3)
  mask_full_idx <- which(mask_arr)

  fwhm_vox <- .resolve_smoothness_to_fwhm_vox(smoothness, template = template)

  seeds <- .null_index_seeds(n_perm, seed = seed)
  bare <- make_null_fun_mc_fwhm(x0 = NULL, y0 = NULL, z0 = NULL,
                                mask_full_idx0 = mask_full_idx,
                                side = side, fwhm_vox = fwhm_vox)
  draw <- function(b) {
    .check_perm_index(b, n_perm)
    .draw_with_seed(seeds[b], bare, b)
  }

  new_null_spec(
    stat_vol = template, mask = mask_arr, draw = draw,
    n_perm = n_perm, df = NA_real_, kind = "mc_fwhm",
    meta = list(smoothness = smoothness, side = side, seeds = seeds)
  )
}

#' Monte-Carlo ACF null (AFNI-style autocorrelation)
#'
#' Generates null draws by simulating fields with a mixed Gaussian +
#' monoexponential autocorrelation model. Returns a [`null_spec`].
#'
#' @inheritParams null_mc_fwhm
#' @param smoothness A [`smoothness_acf()`] object, a numeric vector
#'   `c(a, b, c)` in millimetres, or a legacy list with `a`, `fwhm_vox`, and
#'   `lambda_vox`.
#' @return A [`null_spec`].
#' @export
null_mc_acf <- function(template, mask, smoothness, n_perm = 5000, seed = NULL) {
  if (n_perm < 1L) stop("n_perm must be >= 1")
  dims3 <- dim(template)
  if (is.null(dims3) || length(dims3) != 3L) stop("template must be 3D")
  side <- max(dims3)
  mask_arr <- array(as.logical(mask), dim = dims3)
  mask_full_idx <- which(mask_arr)

  acf_params <- .resolve_smoothness_to_acf(smoothness, template = template)

  seeds <- .null_index_seeds(n_perm, seed = seed)
  bare <- make_null_fun_mc_acf(x0 = NULL, y0 = NULL, z0 = NULL,
                               mask_full_idx0 = mask_full_idx,
                               side = side, acf_params = acf_params)
  draw <- function(b) {
    .check_perm_index(b, n_perm)
    .draw_with_seed(seeds[b], bare, b)
  }

  new_null_spec(
    stat_vol = template, mask = mask_arr, draw = draw,
    n_perm = n_perm, df = NA_real_, kind = "mc_acf",
    meta = list(smoothness = smoothness, side = side, seeds = seeds)
  )
}

#' Parcel-level permutation null
#'
#' Used internally by [hier_scan()]. For each draw, voxel signs are flipped
#' independently at the voxel level (the same family used by
#' [null_voxel_signflip()] but exposed for hierarchical workflows).
#'
#' @param stat_vol Observed (canonicalized) Z-equivalent map.
#' @param parcels A `ClusteredNeuroVol`, integer label array, or vector
#'   aligned with `stat_vol`.
#' @param mask Optional logical 3D mask.
#' @param n_perm Number of draws.
#' @param seed Optional random seed.
#' @return A [`null_spec`] with `kind = "parcel_perm"`.
#' @export
null_parcel_perm <- function(stat_vol, parcels = NULL, mask = NULL,
                             n_perm = 5000, seed = NULL) {
  ns <- null_voxel_signflip(stat_vol = stat_vol, mask = mask,
                            n_perm = n_perm, seed = seed)
  ns$kind <- "parcel_perm"
  ns$meta$parcels <- parcels
  ns
}

#' External observed/null statistic maps
#'
#' Wraps an observed statistic map and externally generated null maps in a
#' [`null_spec`]. This is the generic bridge for MVPA, searchlight, PLS, ASCA,
#' or any workflow that already has matched observed and permuted statistic
#' maps.
#'
#' @param stat_vol Observed canonicalized Z-equivalent statistic map.
#' @param null_maps Null maps. Accepts a list of 3D maps, a 4D array whose
#'   fourth dimension indexes permutations, or a numeric matrix with one row per
#'   permutation and either `sum(mask)` or `prod(dim(stat_vol))` columns.
#' @param mask Optional logical 3D mask. Defaults to finite voxels in
#'   `stat_vol`.
#' @param kind Character kind tag for the resulting [`null_spec`].
#' @param meta Optional named metadata list.
#'
#' @return A [`null_spec`] with deterministic `draw(b)` over supplied maps.
#' @export
null_from_maps <- function(stat_vol,
                           null_maps,
                           mask = NULL,
                           kind = "external_maps",
                           meta = list()) {
  dims <- dim(stat_vol)
  if (is.null(dims) || length(dims) != 3L) stop("stat_vol must be a 3D volume-like object")
  n_full <- prod(dims)

  if (is.null(mask)) {
    mask_arr <- array(is.finite(as.numeric(stat_vol)), dim = dims)
  } else {
    if (length(mask) != n_full) stop("mask must match stat_vol dimensions")
    mask_arr <- array(as.logical(mask), dim = dims)
  }
  mask_idx <- which(mask_arr)
  if (!length(mask_idx)) stop("mask is empty")
  if (any(!is.finite(as.numeric(stat_vol)[mask_idx]))) {
    stop("stat_vol contains non-finite values in mask")
  }
  if (!is.list(meta)) stop("meta must be a list")

  source <- NULL
  n_perm <- NULL
  draw <- NULL

  if (is.list(null_maps)) {
    if (!length(null_maps)) stop("null_maps list must be non-empty")
    source <- "list"
    n_perm <- length(null_maps)
    draws <- lapply(seq_len(n_perm), function(b) {
      .coerce_null_draw(null_maps[[b]], n_full = n_full, mask_idx = mask_idx)[mask_idx]
    })
    draw <- function(b) {
      if (b < 1L || b > n_perm) stop("Permutation index out of range")
      draws[[b]]
    }
  } else if (is.array(null_maps) && length(dim(null_maps)) == 4L) {
    d4 <- dim(null_maps)
    if (!identical(as.integer(d4[1:3]), as.integer(dims))) {
      stop("4D null_maps spatial dimensions must match stat_vol")
    }
    source <- "array4d"
    n_perm <- d4[4L]
    draw <- function(b) {
      if (b < 1L || b > n_perm) stop("Permutation index out of range")
      as.numeric(null_maps[, , , b])[mask_idx]
    }
  } else if (is.matrix(null_maps)) {
    source <- "matrix"
    n_perm <- nrow(null_maps)
    if (n_perm < 1L) stop("null_maps matrix must have at least one row")
    if (!ncol(null_maps) %in% c(length(mask_idx), n_full)) {
      stop("null_maps matrix must have n_mask or full-volume columns")
    }
    draw <- function(b) {
      if (b < 1L || b > n_perm) stop("Permutation index out of range")
      row <- as.numeric(null_maps[b, ])
      if (length(row) == n_full) row[mask_idx] else row
    }
  } else {
    stop("null_maps must be a list, 4D array, or numeric matrix")
  }

  test <- draw(1L)
  if (!is.numeric(test) || length(test) != length(mask_idx)) {
    stop("null_maps draw did not produce a numeric vector of length n_mask")
  }
  if (any(!is.finite(test))) stop("null_maps contain non-finite values in mask")

  new_null_spec(
    stat_vol = stat_vol,
    mask = mask_arr,
    draw = draw,
    n_perm = n_perm,
    df = NA_real_,
    kind = kind,
    meta = c(meta, list(source = source, n_mask = length(mask_idx)))
  )
}

#' @keywords internal
.resolve_smoothness_to_fwhm_vox <- function(s, template) {
  if (inherits(s, "neurothresh_smoothness_fwhm")) {
    spacing <- .vox_size_mm(template)
    return(.fwhm_mm_to_vox(s$fwhm_mm, spacing))
  }
  if (inherits(s, "neurothresh_smoothness")) {
    stop("null_mc_fwhm() requires a smoothness_fwhm() (got ",
         class(s)[1L], ")")
  }
  if (is.numeric(s)) {
    if (length(s) == 1L) s <- rep(s, 3L)
    if (length(s) != 3L) stop("numeric smoothness must have length 1 or 3")
    return(s)
  }
  stop("smoothness must be a smoothness_fwhm() or a length-1/3 numeric")
}

#' @keywords internal
.resolve_smoothness_to_acf <- function(s, template) {
  if (inherits(s, "neurothresh_smoothness_acf")) {
    spacing <- .vox_size_mm(template)
    return(list(
      a = s$a,
      fwhm_vox = .fwhm_mm_to_vox(s$b, spacing),
      lambda_vox = as.numeric(s$c) / spacing
    ))
  }
  if (inherits(s, "neurothresh_smoothness")) {
    stop("null_mc_acf() requires a smoothness_acf() (got ",
         class(s)[1L], ")")
  }
  if (is.list(s) && all(c("a", "fwhm_vox", "lambda_vox") %in% names(s))) {
    return(s)
  }
  if (is.numeric(s) && length(s) == 3L) {
    spacing <- .vox_size_mm(template)
    return(list(
      a = s[1L],
      fwhm_vox = .fwhm_mm_to_vox(s[2L], spacing),
      lambda_vox = as.numeric(s[3L]) / spacing
    ))
  }
  stop("smoothness must be a smoothness_acf(), a list(a, fwhm_vox, lambda_vox), or a length-3 numeric c(a, b, c)")
}

#' @keywords internal
.check_perm_index <- function(b, n_perm) {
  if (length(b) != 1L || is.na(b) || b < 1L || b > n_perm) {
    stop("Permutation index out of range")
  }
  invisible(as.integer(b))
}

#' @keywords internal
.null_index_seeds <- function(n_perm, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sample.int(.Machine$integer.max, n_perm)
}

#' @keywords internal
.draw_with_seed <- function(seed, draw_fun, b) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  draw_fun(b)
}
