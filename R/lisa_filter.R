#' LISA bilateral filter specification
#'
#' Creates a filter specification for LISA local adaptive inference. The
#' default backend delegates to `neuroim2::bilateral_filter()`; a small
#' reference backend is available for tests and exact mask-aware semantics on
#' small arrays.
#'
#' @param radius Integer neighborhood half-width in voxels.
#' @param sigma_space Spatial Gaussian sigma in voxel units.
#' @param sigma_range Range/intensity sigma. With the `reference` backend and a
#'   supplied range scale, the effective range sigma is
#'   `sigma_range * range_scale`.
#' @param n_iter Number of repeated filtering iterations.
#' @param scale Range-scale policy. `"sd_first_30_nulls"` estimates one fixed
#'   range scale from early null maps. `"per_map"` uses each map's own masked
#'   standard deviation. `"none"` uses a unit range scale.
#' @param boundary Boundary policy. `"neuroim2"` faithfully delegates to
#'   `neuroim2::bilateral_filter()`. `"mask_normalized"` excludes out-of-mask
#'   neighbors in the reference backend. `"lisa2018"` is reserved for the
#'   paper-style median fallback and currently errors.
#'
#' @return Object of class `"neurothresh_lisa_filter"`.
#' @export
lisa_bilateral <- function(radius = 2L,
                           sigma_space = 2,
                           sigma_range = 2,
                           n_iter = 2L,
                           scale = c("sd_first_30_nulls", "per_map", "none"),
                           boundary = c("neuroim2", "mask_normalized", "lisa2018")) {
  scale <- match.arg(scale)
  boundary <- match.arg(boundary)
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
      scale = scale,
      boundary = boundary
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
#' @param backend Filtering backend: `"neuroim2"` or `"reference"`.
#' @param range_scale Optional fixed range scale for the reference backend.
#'
#' @return Filtered volume-like object aligned to `stat_vol`.
#' @export
lisa_apply_filter <- function(stat_vol,
                              mask = NULL,
                              filter = lisa_bilateral(),
                              backend = c("neuroim2", "reference", "adjoin_sparse"),
                              range_scale = NULL) {
  backend <- match.arg(backend)
  .lisa_apply_filter(stat_vol, mask = mask, filter = filter,
                     backend = backend, range_scale = range_scale)
}

#' @keywords internal
.lisa_apply_filter <- function(stat_vol, mask, filter, backend, range_scale = NULL) {
  if (!is_lisa_filter(filter)) stop("filter must be created by lisa_bilateral()")
  if (filter$boundary == "lisa2018") {
    stop("boundary = 'lisa2018' is not implemented yet")
  }
  if (backend == "adjoin_sparse") {
    stop("backend = 'adjoin_sparse' is planned but not implemented yet")
  }
  if (backend == "neuroim2" && filter$boundary != "neuroim2") {
    backend <- "reference"
  }
  if (backend == "neuroim2") {
    return(.lisa_filter_neuroim2(stat_vol, mask, filter))
  }
  .lisa_filter_reference(stat_vol, mask, filter, range_scale = range_scale)
}

#' @keywords internal
.lisa_filter_neuroim2 <- function(stat_vol, mask, filter) {
  vol <- .lisa_as_neurovol(stat_vol)
  mask_vol <- .lisa_as_maskvol(mask, vol)
  out <- vol
  for (i in seq_len(filter$n_iter)) {
    out <- neuroim2::bilateral_filter(
      out,
      mask = mask_vol,
      spatial_sigma = .lisa_sigma_space_mm(filter$sigma_space, out),
      intensity_sigma = filter$sigma_range,
      window = filter$radius
    )
  }
  .array_to_vol(as.numeric(out), dim(stat_vol), stat_vol)
}

#' @keywords internal
.lisa_filter_reference <- function(stat_vol, mask, filter, range_scale = NULL) {
  dims <- dim(stat_vol)
  if (is.null(dims) || length(dims) != 3L) stop("stat_vol must be 3D")
  z <- as.numeric(stat_vol)
  if (is.null(mask)) mask <- array(is.finite(z), dim = dims)
  mask <- array(as.logical(mask), dim = dims)
  mask_idx <- which(mask)
  if (!length(mask_idx)) stop("mask is empty")

  spacing <- .ss_voxel_spacing(stat_vol)
  sigma_space <- filter$sigma_space * mean(spacing)
  cur <- z
  cur[!is.finite(cur)] <- 0
  for (iter in seq_len(filter$n_iter)) {
    scale <- range_scale
    if (is.null(scale)) {
      scale <- stats::sd(cur[mask_idx])
    }
    if (!is.finite(scale) || scale <= 0) scale <- 1
    cur <- .lisa_filter_reference_once(
      cur, dims = dims, mask = mask, radius = filter$radius,
      sigma_space = sigma_space,
      sigma_range = filter$sigma_range * scale,
      spacing = spacing
    )
  }
  .array_to_vol(cur, dims, stat_vol)
}

#' @keywords internal
.lisa_filter_reference_once <- function(z, dims, mask, radius, sigma_space,
                                        sigma_range, spacing) {
  out <- rep(0, length(z))
  coords <- which(mask, arr.ind = TRUE)
  mask_idx <- which(mask)
  offsets <- expand.grid(
    dx = seq.int(-radius, radius),
    dy = seq.int(-radius, radius),
    dz = seq.int(-radius, radius)
  )
  d2 <- (offsets$dx * spacing[1L])^2 +
    (offsets$dy * spacing[2L])^2 +
    (offsets$dz * spacing[3L])^2
  spatial_w <- exp(-d2 / (2 * sigma_space^2))
  range_denom <- 2 * sigma_range^2

  for (r in seq_along(mask_idx)) {
    cx <- coords[r, 1L]
    cy <- coords[r, 2L]
    cz <- coords[r, 3L]
    center <- z[mask_idx[r]]
    vals <- numeric(0)
    wts <- numeric(0)
    for (o in seq_len(nrow(offsets))) {
      x <- cx + offsets$dx[o]
      y <- cy + offsets$dy[o]
      zz <- cz + offsets$dz[o]
      if (x < 1L || x > dims[1L] || y < 1L || y > dims[2L] ||
          zz < 1L || zz > dims[3L]) {
        next
      }
      if (!mask[x, y, zz]) next
      idx <- x + (y - 1L) * dims[1L] + (zz - 1L) * dims[1L] * dims[2L]
      val <- z[idx]
      w <- spatial_w[o] * exp(-((val - center)^2) / range_denom)
      vals <- c(vals, val)
      wts <- c(wts, w)
    }
    out[mask_idx[r]] <- if (sum(wts) > 0) sum(wts * vals) / sum(wts) else center
  }
  out
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
