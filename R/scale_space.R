#' Define a Gaussian Scale-Space Grid
#'
#' @param fwhm_mm Numeric vector of Gaussian FWHM scales in millimeters.
#' @param include_identity Include the unsmoothed voxelwise endpoint.
#' @param include_global Include the global endpoint.
#' @param report Metadata requested by summaries.
#'
#' @return A \code{scale_space_spec} object.
#'
#' @export
scale_space <- function(fwhm_mm,
                        include_identity = TRUE,
                        include_global = TRUE,
                        report = c("fwhm", "effective_basis")) {
  if (missing(fwhm_mm) || is.null(fwhm_mm)) fwhm_mm <- numeric(0)
  if (!is.numeric(fwhm_mm) || any(!is.finite(fwhm_mm)) || any(fwhm_mm <= 0)) {
    stop("fwhm_mm must contain positive finite values")
  }
  fwhm_mm <- sort(unique(as.numeric(fwhm_mm)))

  scales <- data.frame(
    type = character(),
    fwhm_mm = numeric(),
    sigma_mm = numeric(),
    stringsAsFactors = FALSE
  )
  if (isTRUE(include_global)) {
    scales <- rbind(scales, data.frame(type = "global", fwhm_mm = Inf, sigma_mm = Inf))
  }
  if (length(fwhm_mm)) {
    scales <- rbind(
      scales,
      data.frame(
        type = "gaussian",
        fwhm_mm = fwhm_mm,
        sigma_mm = fwhm_mm / sqrt(8 * log(2))
      )
    )
  }
  if (isTRUE(include_identity)) {
    scales <- rbind(scales, data.frame(type = "identity", fwhm_mm = 0, sigma_mm = 0))
  }

  structure(
    list(scales = scales, report = report),
    class = "scale_space_spec"
  )
}

#' Local Energy Statistic Channel
#'
#' @export
local_energy <- function() {
  structure(list(name = "local_energy"), class = c("ss_local_energy", "ss_channel"))
}

#' Signed Coherence Statistic Channel
#'
#' @param tail Direction for signed evidence.
#'
#' @export
signed_coherence <- function(tail = c("pos", "neg", "two")) {
  tail <- match.arg(tail)
  structure(list(name = "signed_coherence", tail = tail),
            class = c("ss_signed_coherence", "ss_channel"))
}

#' Positive-Direction Local Energy Statistic Channel
#'
#' @param direction Direction of evidence to square. Use "pos" for positive
#'   evidence or "neg" for negative evidence.
#'
#' @export
positive_energy <- function(direction = c("pos", "neg")) {
  direction <- match.arg(direction)
  structure(list(name = "positive_energy", direction = direction),
            class = c("ss_positive_energy", "ss_channel"))
}

#' Band-Pass Local Energy Statistic Channel
#'
#' Explicit opt-in spatial band-pass channel using a difference of smoothers
#' around the requested scale.
#'
#' @param ratio Multiplicative scale ratio around the target FWHM.
#'
#' @export
bandpass_energy <- function(ratio = sqrt(2)) {
  if (!is.numeric(ratio) || length(ratio) != 1L || !is.finite(ratio) || ratio <= 1) {
    stop("ratio must be a finite scalar > 1")
  }
  structure(list(name = "bandpass_energy", ratio = ratio),
            class = c("ss_bandpass_energy", "ss_channel"))
}

#' @keywords internal
.ss_channel_name <- function(channel) {
  if (!inherits(channel, "ss_channel")) stop("statistics must be channel objects")
  if (inherits(channel, "ss_signed_coherence")) {
    return(paste0("signed_coherence_", channel$tail))
  }
  if (inherits(channel, "ss_positive_energy")) {
    return(paste0("positive_energy_", channel$direction))
  }
  if (inherits(channel, "ss_bandpass_energy")) {
    return("bandpass_energy")
  }
  channel$name
}

#' @keywords internal
.ss_apply_scale <- function(z_mask, scale, dims, mask, mask_idx, template,
                            backend = "grid_gaussian", graph_cache = NULL) {
  type <- scale$type
  if (identical(type, "identity")) return(z_mask)
  if (identical(type, "global")) return(rep(mean(z_mask), length(z_mask)))

  if (identical(backend, "graph_heat")) {
    if (is.null(graph_cache)) stop("graph_cache is required for graph_heat backend")
    return(.ss_graph_heat_smooth(z_mask, scale, graph_cache))
  }

  full <- numeric(prod(dims))
  full[mask_idx] <- z_mask
  sigma_vox <- .ss_sigma_vox(scale$sigma_mm, template)
  sm <- gaussian_smooth3d_cpp(full, as.integer(dims), sigma_xyz = sigma_vox)

  w_full <- numeric(prod(dims))
  w_full[mask_idx] <- 1
  ww <- gaussian_smooth3d_cpp(w_full, as.integer(dims), sigma_xyz = sigma_vox)
  out <- sm[mask_idx] / pmax(ww[mask_idx], .Machine$double.eps)
  out[!is.finite(out)] <- 0
  out
}

#' @keywords internal
.ss_sigma_vox <- function(sigma_mm, template) {
  if (!is.finite(sigma_mm) || sigma_mm <= 0) return(rep(0, 3))
  spacing <- .ss_voxel_spacing(template)
  rep(sigma_mm, 3) / spacing
}

#' @keywords internal
.ss_voxel_spacing <- function(template) {
  sp <- NULL
  sp <- tryCatch({
    if (requireNamespace("neuroim2", quietly = TRUE) &&
        inherits(template, "NeuroVol")) {
      neuroim2::spacing(neuroim2::space(template))
    } else {
      NULL
    }
  }, error = function(e) NULL)
  if (is.null(sp)) sp <- attr(template, "spacing")
  if (is.null(sp)) sp <- rep(1, 3)
  if (length(sp) == 1L) sp <- rep(sp, 3)
  as.numeric(sp)[1:3]
}

#' @keywords internal
.ss_eval_channel_scale <- function(z_mask, channel, scale, dims, mask, mask_idx, template,
                                   backend = "grid_gaussian", graph_cache = NULL) {
  if (inherits(channel, "ss_local_energy")) {
    if (identical(scale$type, "global")) {
      return(rep(sum(z_mask^2), length(z_mask)))
    }
    return(.ss_apply_scale(z_mask^2, scale, dims, mask, mask_idx, template,
                           backend = backend, graph_cache = graph_cache))
  }

  if (inherits(channel, "ss_positive_energy")) {
    qz <- if (channel$direction == "neg") -z_mask else z_mask
    e <- pmax(qz, 0)^2
    if (identical(scale$type, "global")) {
      return(rep(sum(e), length(z_mask)))
    }
    return(.ss_apply_scale(e, scale, dims, mask, mask_idx, template,
                           backend = backend, graph_cache = graph_cache))
  }

  if (inherits(channel, "ss_bandpass_energy")) {
    if (!identical(scale$type, "gaussian")) {
      return(rep(0, length(z_mask)))
    }
    fine <- scale
    coarse <- scale
    fine$fwhm_mm <- scale$fwhm_mm / channel$ratio
    coarse$fwhm_mm <- scale$fwhm_mm * channel$ratio
    fine$sigma_mm <- fine$fwhm_mm / sqrt(8 * log(2))
    coarse$sigma_mm <- coarse$fwhm_mm / sqrt(8 * log(2))
    bf <- .ss_apply_scale(z_mask, fine, dims, mask, mask_idx, template,
                          backend = backend, graph_cache = graph_cache)
    bc <- .ss_apply_scale(z_mask, coarse, dims, mask, mask_idx, template,
                          backend = backend, graph_cache = graph_cache)
    return((bf - bc)^2)
  }

  if (inherits(channel, "ss_signed_coherence")) {
    if (channel$tail == "two") {
      # Two-sided signed coherence: |S_rho z|. Smooth the signed field, then
      # take the magnitude. This is what makes coherent same-sign blobs add
      # under smoothing while opposite-sign blobs cancel. Computing
      # S_rho |z| (magnitude before smoothing) is a different statistic that
      # behaves like a linear magnitude-energy channel.
      if (identical(scale$type, "global")) {
        return(rep(abs(sum(z_mask)) / sqrt(length(z_mask)), length(z_mask)))
      }
      smoothed <- .ss_apply_scale(z_mask, scale, dims, mask, mask_idx, template,
                                  backend = backend, graph_cache = graph_cache)
      return(abs(smoothed))
    }
    zz <- if (channel$tail == "neg") -z_mask else z_mask
    if (identical(scale$type, "global")) {
      return(rep(sum(zz) / sqrt(length(zz)), length(zz)))
    }
    return(.ss_apply_scale(zz, scale, dims, mask, mask_idx, template,
                           backend = backend, graph_cache = graph_cache))
  }

  stop("unsupported scale-space statistic channel")
}

#' @keywords internal
.ss_eval_all <- function(z_mask, scales, channels, dims, mask, mask_idx, template,
                         backend = "grid_gaussian", graph_cache = NULL) {
  vals <- vector("list", nrow(scales) * length(channels))
  meta <- vector("list", length(vals))
  k <- 1L
  for (si in seq_len(nrow(scales))) {
    for (ci in seq_along(channels)) {
      vals[[k]] <- .ss_eval_channel_scale(z_mask, channels[[ci]], scales[si, ],
                                          dims, mask, mask_idx, template,
                                          backend = backend, graph_cache = graph_cache)
      meta[[k]] <- data.frame(
        scale_index = si,
        channel_index = ci,
        scale_type = scales$type[si],
        fwhm_mm = scales$fwhm_mm[si],
        channel = .ss_channel_name(channels[[ci]]),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  list(values = vals, meta = do.call(rbind, meta))
}

#' Build a weighted 6-connected graph cache for a mask.
#'
#' Edge weights are inverse squared voxel spacings, so the graph Laplacian is
#' on the same physical millimetre scale as the Gaussian FWHM specification.
#'
#' @keywords internal
#' @noRd
.ss_build_graph_cache <- function(mask, dims, template = NULL) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("graph_heat backend requires the Matrix package")
  }
  dims <- as.integer(dims)
  if (length(dims) != 3L) stop("dims must have length 3")
  if (length(mask) != prod(dims)) stop("mask must match dims")

  mask_idx <- which(mask)
  n <- length(mask_idx)
  spacing <- .ss_voxel_spacing(template)
  pos_full <- integer(prod(dims))
  pos_full[mask_idx] <- seq_len(n)
  coords <- arrayInd(mask_idx, .dim = dims)

  ii <- integer(0)
  jj <- integer(0)
  xx <- numeric(0)
  strides <- c(1L, dims[1], dims[1] * dims[2])
  for (axis in seq_len(3L)) {
    source <- mask_idx[coords[, axis] < dims[axis]]
    if (!length(source)) next
    target <- source + strides[axis]
    keep <- pos_full[target] > 0L
    if (!any(keep)) next
    src_pos <- pos_full[source[keep]]
    tgt_pos <- pos_full[target[keep]]
    weight <- 1 / (spacing[axis]^2)
    ii <- c(ii, src_pos, tgt_pos)
    jj <- c(jj, tgt_pos, src_pos)
    xx <- c(xx, rep(weight, 2L * length(src_pos)))
  }

  A <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n))
  deg <- Matrix::rowSums(A)
  L <- Matrix::Diagonal(x = as.numeric(deg)) - A
  lambda_max_bound <- if (n) 2 * max(as.numeric(deg)) else 0
  list(
    adjacency = A,
    laplacian = L,
    degree = as.numeric(deg),
    lambda_max_bound = lambda_max_bound,
    spacing = spacing,
    dims = dims,
    n = n
  )
}

#' Continuous-time graph heat smoother.
#'
#' Applies \code{expm(-t L) z}, where \code{L} is the weighted graph Laplacian
#' and \code{t = sigma_mm^2 / 2}. On a regular grid away from boundaries, this
#' matches the Gaussian heat-equation calibration in physical millimetres.
#'
#' @keywords internal
#' @noRd
.ss_graph_heat_smooth <- function(z_mask, scale, graph_cache) {
  if (scale$type != "gaussian") return(z_mask)
  if (length(z_mask) != graph_cache$n) stop("z_mask length must match graph cache")
  heat_time <- .ss_graph_heat_time(scale)
  .ss_heat_apply_chebyshev(
    graph_cache$laplacian,
    as.numeric(z_mask),
    time = heat_time,
    lambda_max = graph_cache$lambda_max_bound
  )
}

#' @keywords internal
#' @noRd
.ss_graph_heat_time <- function(scale) {
  if (!is.finite(scale$sigma_mm) || scale$sigma_mm <= 0) return(0)
  scale$sigma_mm^2 / 2
}

#' @keywords internal
#' @noRd
.ss_heat_apply_chebyshev <- function(L, z, time, lambda_max,
                                    tol = 1e-10, min_degree = 24L) {
  if (!is.finite(time) || time <= 0 || !is.finite(lambda_max) || lambda_max <= 0) {
    return(as.numeric(z))
  }
  a <- time * lambda_max / 2
  if (a <= 0) return(as.numeric(z))

  matvec <- function(v) {
    as.numeric((2 / lambda_max) * (L %*% v) - v)
  }

  max_degree <- max(min_degree, ceiling(a + 12 * sqrt(a + 1) + 30))
  t0 <- as.numeric(z)
  out <- besselI(a, 0, expon.scaled = TRUE) * t0

  t1 <- matvec(t0)
  coef <- -2 * besselI(a, 1, expon.scaled = TRUE)
  out <- out + coef * t1
  small_tail <- as.integer(abs(coef) < tol)

  for (k in seq.int(2L, max_degree)) {
    tk <- 2 * matvec(t1) - t0
    coef <- 2 * if (k %% 2L) -1 else 1
    coef <- coef * besselI(a, k, expon.scaled = TRUE)
    out <- out + coef * tk
    t0 <- t1
    t1 <- tk

    if (k > a && abs(coef) < tol) {
      small_tail <- small_tail + 1L
      if (small_tail >= 6L) break
    } else {
      small_tail <- 0L
    }
  }

  out[!is.finite(out)] <- 0
  as.numeric(out)
}
