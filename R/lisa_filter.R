#' LISA bilateral filter specification
#'
#' Creates a filter specification for LISA local adaptive inference. LISA uses
#' the `neuroim2::bilateral_filter()` implementation, whose neighborhoods are
#' mask-normalized so only in-mask, in-bounds neighbors contribute to each
#' filtered statistic.
#'
#' @param radius Integer neighborhood half-width in voxels.
#' @param sigma_space Spatial Gaussian sigma in voxel units.
#' @param sigma_range Range/intensity sigma. With a fixed range scale, the
#'   effective range sigma is `sigma_range * range_scale`.
#' @param n_iter Number of repeated filtering iterations.
#' @param scale Range-scale policy. `"sd_first_30_nulls"` estimates one fixed
#'   range scale from early null maps. `"per_map"` uses each map's own masked
#'   standard deviation. `"none"` uses a unit range scale.
#'
#' @return Object of class `"neurothresh_lisa_filter"`.
#' @export
lisa_bilateral <- function(radius = 2L,
                           sigma_space = 2,
                           sigma_range = 2,
                           n_iter = 2L,
                           scale = c("sd_first_30_nulls", "per_map", "none")) {
  scale <- match.arg(scale)
  if (!is.numeric(radius) || length(radius) != 1L || radius < 1L) {
    stop("radius must be a positive integer")
  }
  if (!is.numeric(sigma_space) || length(sigma_space) != 1L || sigma_space <= 0) {
    stop("sigma_space must be positive")
  }
  if (!is.numeric(sigma_range) || length(sigma_range) != 1L || sigma_range <= 0) {
    stop("sigma_range must be positive")
  }
  if (!is.numeric(n_iter) || length(n_iter) != 1L || n_iter < 1L) {
    stop("n_iter must be a positive integer")
  }
  structure(
    list(
      type = "bilateral",
      radius = as.integer(radius),
      sigma_space = as.numeric(sigma_space),
      sigma_range = as.numeric(sigma_range),
      n_iter = as.integer(n_iter),
      scale = scale
    ),
    class = "neurothresh_lisa_filter"
  )
}

#' Test whether an object is a LISA filter specification
#'
#' @param x Any R object.
#' @return Logical scalar.
#' @export
is_lisa_filter <- function(x) inherits(x, "neurothresh_lisa_filter")

#' Apply a LISA filter
#'
#' Applies a LISA filter specification to a 3D statistic map. This function is
#' mostly intended for diagnostics and tests; `lisa()` applies it internally to
#' observed and null maps.
#'
#' @param stat_vol 3D statistic map.
#' @param mask Optional logical 3D mask.
#' @param filter A LISA filter specification from [lisa_bilateral()].
#' @param range_scale Optional fixed range scale. If `NULL`, the filter uses
#'   the map's masked standard deviation.
#'
#' @return Filtered volume-like object aligned to `stat_vol`.
#' @export
lisa_apply_filter <- function(stat_vol,
                              mask = NULL,
                              filter = lisa_bilateral(),
                              range_scale = NULL) {
  .lisa_apply_filter(stat_vol, mask = mask, filter = filter,
                     range_scale = range_scale)
}

#' @keywords internal
.lisa_apply_filter <- function(stat_vol, mask, filter, range_scale = NULL) {
  if (!is_lisa_filter(filter)) stop("filter must be created by lisa_bilateral()")
  .lisa_filter_neuroim2(stat_vol, mask, filter, range_scale = range_scale)
}

#' @keywords internal
.lisa_filter_neuroim2 <- function(stat_vol, mask, filter, range_scale = NULL) {
  formals <- names(formals(neuroim2::bilateral_filter))
  if (!("range_scale" %in% formals)) {
    stop("LISA requires neuroim2::bilateral_filter() with range_scale support and mask-normalized neighborhoods; install neuroim2 >= 0.16.0")
  }
  vol <- .lisa_as_neurovol(stat_vol)
  mask_vol <- .lisa_as_maskvol(mask, vol)
  out <- vol
  for (i in seq_len(filter$n_iter)) {
    out <- neuroim2::bilateral_filter(
      out,
      mask = mask_vol,
      spatial_sigma = .lisa_sigma_space_mm(filter$sigma_space, out),
      intensity_sigma = filter$sigma_range,
      window = filter$radius,
      range_scale = range_scale
    )
  }
  .array_to_vol(as.numeric(out), dim(stat_vol), stat_vol)
}

#' @keywords internal
.lisa_as_neurovol <- function(stat_vol) {
  if (inherits(stat_vol, "NeuroVol")) return(stat_vol)
  dims <- dim(stat_vol)
  sp <- neuroim2::NeuroSpace(dims, spacing = rep(1, 3))
  neuroim2::NeuroVol(array(as.numeric(stat_vol), dim = dims), sp)
}

#' @keywords internal
.lisa_as_maskvol <- function(mask, vol) {
  dims <- dim(vol)
  if (is.null(mask)) {
    mask <- array(is.finite(as.numeric(vol)), dim = dims)
  } else {
    mask <- array(as.logical(mask), dim = dims)
  }
  neuroim2::LogicalNeuroVol(mask, neuroim2::space(vol))
}

#' @keywords internal
.lisa_sigma_space_mm <- function(sigma_space, vol) {
  spacing <- tryCatch(neuroim2::spacing(vol), error = function(e) rep(1, 3))
  sigma_space * mean(spacing[seq_len(3L)])
}

#' @keywords internal
.lisa_estimate_range_scale <- function(null_spec, mask_idx, n_full, n = 30L) {
  n <- min(as.integer(n), null_spec$n_perm)
  if (n < 1L) return(1)
  sds <- numeric(n)
  for (b in seq_len(n)) {
    draw <- .coerce_null_draw(null_spec$draw(b), n_full = n_full, mask_idx = mask_idx)
    sds[b] <- stats::sd(draw[mask_idx])
  }
  out <- mean(sds[is.finite(sds) & sds > 0])
  if (!is.finite(out) || out <= 0) 1 else out
}
