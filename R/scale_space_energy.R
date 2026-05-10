#' Scale-Space Localized Energy Inference
#'
#' Scans a statistic map across a prespecified Gaussian scale-space and controls
#' FWER over scale, location, and statistic channel using synchronized null maps.
#'
#' @param stat_vol 3D statistic map, or a [`null_spec`] when `null_spec` is not
#'   supplied separately.
#' @param mask Optional logical analysis mask.
#' @param scales A \code{scale_space_spec} object.
#' @param statistics List of channel objects such as \code{local_energy()} and
#'   \code{signed_coherence()}.
#' @param null_fun Function \code{function(b)} returning a null statistic draw.
#'   It should be deterministic by index because calibration reuses draw
#'   \code{b} after estimating common observed-plus-null moments.
#' @param null_spec Optional [`null_spec`] supplying `stat_vol`, `mask`,
#'   deterministic `draw(b)`, and `n_perm`.
#' @param n_perm Number of permutations/null draws.
#' @param studentize Whether to studentize transformed fields using null mean
#'   and standard deviation estimated symmetrically from the observed map plus
#'   the null draws. Defaults to TRUE.
#' @param alpha FWER level used for summaries.
#' @param seed Optional random seed. Only meaningful when \code{null_fun} draws
#'   on the global RNG; the bundled \code{make_null_fun_*} helpers are
#'   deterministic by index and ignore this seed.
#' @param backend Smoothing backend. \code{"grid_gaussian"} uses a separable,
#'   mask-normalized 3D Gaussian convolution. \code{"graph_heat"} uses a
#'   weighted 6-neighbor mask graph and applies the continuous-time heat
#'   operator \code{expm(-t L)} with \code{t = sigma_mm^2 / 2}.
#'
#' @return An object of class \code{scale_space_energy_result}.
#'
#' @export
scale_space_energy <- function(stat_vol = NULL,
                               mask = NULL,
                               scales = scale_space(c(6, 10, 14, 18, 24, 32, 48)),
                               statistics = list(local_energy(), signed_coherence(tail = "two")),
                               null_fun = NULL,
                               null_spec = NULL,
                               n_perm = NULL,
                               studentize = TRUE,
                               alpha = 0.05,
                               seed = NULL,
                               backend = c("grid_gaussian", "graph_heat")) {
  backend <- match.arg(backend)

  if (is_null_spec(stat_vol) && is.null(null_spec)) {
    null_spec <- stat_vol
    stat_vol <- NULL
  }
  if (!is.null(null_spec)) {
    if (!is_null_spec(null_spec)) stop("null_spec must be a neurothresh_null_spec (see ?null_spec)")
    validate_null_spec(null_spec)
    if (is.null(stat_vol)) stat_vol <- null_spec$stat_vol
    if (is.null(mask)) mask <- null_spec$mask
    if (is.null(null_fun)) null_fun <- null_spec$draw
    if (is.null(n_perm)) n_perm <- null_spec$n_perm
  }
  if (is.null(null_fun)) {
    stop("scale_space_energy() requires null_spec or null_fun for calibrated inference")
  }
  if (is.null(n_perm)) n_perm <- 1000L
  if (!inherits(scales, "scale_space_spec")) stop("scales must be created by scale_space()")
  if (!is.list(statistics) || !length(statistics)) stop("statistics must be a non-empty list")
  if (!is.numeric(n_perm) || length(n_perm) != 1L || n_perm < 1L) stop("n_perm must be >= 1")
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1)")
  }
  if (!is.null(seed)) set.seed(seed)

  dims <- dim(stat_vol)
  if (length(dims) != 3L) stop("stat_vol must be a 3D volume-like object")
  z_full <- as.numeric(stat_vol)
  if (is.null(mask)) {
    mask <- array(is.finite(z_full), dim = dims)
  } else {
    mask <- array(as.logical(mask), dim = dims)
  }
  if (length(mask) != length(z_full)) stop("mask must match stat_vol dimensions")
  mask_idx <- which(mask)
  if (!length(mask_idx)) stop("mask is empty")
  z_obs <- z_full[mask_idx]
  if (any(!is.finite(z_obs))) stop("stat_vol contains non-finite values in mask")

  scale_df <- scales$scales
  channels <- statistics
  graph_cache <- if (backend == "graph_heat") .ss_build_graph_cache(mask, dims, stat_vol) else NULL
  obs_raw <- .ss_eval_all(z_obs, scale_df, channels, dims, mask, mask_idx, stat_vol,
                          backend = backend, graph_cache = graph_cache)
  n_families <- length(obs_raw$values)
  n_mask <- length(mask_idx)
  effective_basis <- .ss_effective_basis(scale_df, n_mask, stat_vol)

  means <- vector("list", n_families)
  m2 <- vector("list", n_families)
  for (j in seq_len(n_families)) {
    means[[j]] <- obs_raw$values[[j]]
    m2[[j]] <- numeric(n_mask)
  }

  for (b in seq_len(n_perm)) {
    z_perm <- .coerce_null_draw(null_fun(b), n_full = prod(dims), mask_idx = mask_idx)[mask_idx]
    cur <- .ss_eval_all(z_perm, scale_df, channels, dims, mask, mask_idx, stat_vol,
                        backend = backend, graph_cache = graph_cache)$values
    for (j in seq_len(n_families)) {
      delta <- cur[[j]] - means[[j]]
      means[[j]] <- means[[j]] + delta / (b + 1L)
      m2[[j]] <- m2[[j]] + delta * (cur[[j]] - means[[j]])
    }
  }

  sds <- vector("list", n_families)
  for (j in seq_len(n_families)) {
    varj <- m2[[j]] / n_perm
    sds[[j]] <- sqrt(pmax(varj, .Machine$double.eps))
  }

  obs_stat <- vector("list", n_families)
  for (j in seq_len(n_families)) {
    if (isTRUE(studentize)) {
      obs_stat[[j]] <- (obs_raw$values[[j]] - means[[j]]) / sds[[j]]
    } else {
      obs_stat[[j]] <- obs_raw$values[[j]]
    }
  }

  max_null <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    z_perm <- .coerce_null_draw(null_fun(b), n_full = prod(dims), mask_idx = mask_idx)[mask_idx]
    cur <- .ss_eval_all(z_perm, scale_df, channels, dims, mask, mask_idx, stat_vol,
                        backend = backend, graph_cache = graph_cache)$values
    mx <- -Inf
    for (j in seq_len(n_families)) {
      tj <- if (isTRUE(studentize)) (cur[[j]] - means[[j]]) / sds[[j]] else cur[[j]]
      mx <- max(mx, tj, na.rm = TRUE)
    }
    max_null[b] <- mx
  }

  p_corr <- lapply(obs_stat, fwer_p_from_maxnull, max_null = max_null)
  threshold <- fwer_threshold_from_maxnull(max_null, alpha = alpha)
  summary <- .ss_summary_table(obs_stat, p_corr, obs_raw$meta, alpha, effective_basis)

  out <- list(
    method = "scale_space_energy",
    observed = list(raw = obs_raw$values, statistic = obs_stat),
    null_moments = list(mean = means, sd = sds),
    p_corr = p_corr,
    max_null = max_null,
    threshold = threshold,
    scales = scales,
    channels = vapply(channels, .ss_channel_name, character(1)),
    family = obs_raw$meta,
    summary = summary,
    selected = .ss_selected_scales(summary, alpha),
    mask = mask,
    mask_idx = mask_idx,
    dims = dims,
    template = stat_vol,
    params = list(alpha = alpha, n_perm = n_perm, studentize = studentize, backend = backend)
  )
  class(out) <- "scale_space_energy_result"
  out
}

#' @keywords internal
.ss_summary_table <- function(obs_stat, p_corr, meta, alpha, effective_basis = NULL) {
  rows <- vector("list", length(obs_stat))
  for (j in seq_along(obs_stat)) {
    p <- p_corr[[j]]
    i_best <- which.min(p)
    rows[[j]] <- cbind(
      meta[j, , drop = FALSE],
      min_p = p[i_best],
      max_stat = max(obs_stat[[j]], na.rm = TRUE),
      any_supported = any(p <= alpha),
      n_supported = sum(p <= alpha)
    )
  }
  out <- do.call(rbind, rows)
  if (!is.null(effective_basis)) {
    out$effective_basis <- effective_basis[out$scale_index]
  }
  rownames(out) <- NULL
  out
}

#' @keywords internal
.ss_effective_basis <- function(scales, n_mask, template) {
  out <- numeric(nrow(scales))
  spacing <- .ss_voxel_spacing(template)
  voxel_volume <- prod(spacing)
  for (i in seq_len(nrow(scales))) {
    if (scales$type[i] == "global") {
      out[i] <- 1
    } else if (scales$type[i] == "identity") {
      out[i] <- n_mask
    } else {
      # First-pass resolution-element approximation. It is monotone in FWHM
      # and bounded, while exact trace-based estimates are reserved for graph
      # or spectral backends.
      resel_volume <- max(scales$fwhm_mm[i]^3, voxel_volume)
      out[i] <- max(1, min(n_mask, n_mask * voxel_volume / resel_volume))
    }
  }
  out
}

#' @keywords internal
.ss_selected_scales <- function(summary, alpha) {
  best_i <- which.min(summary$min_p)
  supported <- summary[summary$any_supported, , drop = FALSE]
  finite_supported <- supported[is.finite(supported$fwhm_mm), , drop = FALSE]
  finest <- if (nrow(finite_supported)) {
    finite_supported[which.min(finite_supported$fwhm_mm), , drop = FALSE]
  } else {
    supported[0, , drop = FALSE]
  }
  list(
    best = summary[best_i, , drop = FALSE],
    finest_detected = finest,
    stable_range = .ss_stable_range(supported)
  )
}

#' @keywords internal
.ss_stable_range <- function(supported) {
  g <- supported[is.finite(supported$fwhm_mm) & supported$scale_type == "gaussian", , drop = FALSE]
  if (nrow(g) < 2L) return(NULL)
  data.frame(
    min_fwhm_mm = min(g$fwhm_mm),
    max_fwhm_mm = max(g$fwhm_mm),
    n_scales = length(unique(g$fwhm_mm))
  )
}

#' @export
print.scale_space_energy_result <- function(x, ...) {
  cat("Scale-space localized energy inference\n")
  cat("  permutations:", x$params$n_perm, "\n")
  cat("  channels:", paste(x$channels, collapse = ", "), "\n")
  cat("  omnibus p_FWER:", format(min(x$summary$min_p), digits = 4), "\n")
  invisible(x)
}

#' @export
summary.scale_space_energy_result <- function(object, ...) {
  object$summary
}

#' Resolution Evidence Curve
#'
#' @param fit A \code{scale_space_energy_result}.
#' @param ... Reserved for future plotting options.
#'
#' @return A data frame with resolution-evidence values.
#'
#' @export
plot_resolution_curve <- function(fit, ...) {
  if (!inherits(fit, "scale_space_energy_result")) stop("fit must be a scale_space_energy_result")
  out <- fit$summary
  out$neg_log10_p <- -log10(out$min_p)
  out
}

#' Extract a Scale-Space Corrected Map
#'
#' @param fit A \code{scale_space_energy_result}.
#' @param scale_index Scale index to extract.
#' @param channel Channel name to extract.
#' @param what One of "p_corr", "statistic", or "raw".
#'
#' @return A volume-like object aligned to the input statistic map.
#'
#' @export
plot_scale_space_map <- function(fit,
                                 scale_index = NULL,
                                 channel = NULL,
                                 what = c("p_corr", "statistic", "raw")) {
  if (!inherits(fit, "scale_space_energy_result")) stop("fit must be a scale_space_energy_result")
  what <- match.arg(what)
  fam <- fit$family
  if (is.null(scale_index)) scale_index <- fit$selected$best$scale_index[1L]
  if (is.null(channel)) channel <- fit$selected$best$channel[1L]
  idx <- which(fam$scale_index == scale_index & fam$channel == channel)[1L]
  if (!length(idx) || is.na(idx)) stop("requested scale/channel not found")
  vals <- switch(what,
    p_corr = fit$p_corr[[idx]],
    statistic = fit$observed$statistic[[idx]],
    raw = fit$observed$raw[[idx]]
  )
  full <- rep(NA_real_, prod(fit$dims))
  full[fit$mask_idx] <- vals
  .array_to_vol(full, fit$dims, fit$template)
}

#' Above-Null Energy Mass Map
#'
#' @param fit A \code{scale_space_energy_result}.
#' @param scale_index Scale index to extract.
#' @param mass Optional cumulative mass cutoff.
#'
#' @return A volume-like normalized above-null energy mass map.
#'
#' @export
energy_mass_map <- function(fit, scale_index = NULL, mass = NULL) {
  if (!inherits(fit, "scale_space_energy_result")) stop("fit must be a scale_space_energy_result")
  fam <- fit$family
  if (is.null(scale_index)) scale_index <- fit$selected$best$scale_index[1L]
  idx <- which(fam$scale_index == scale_index & fam$channel == "local_energy")[1L]
  if (!length(idx) || is.na(idx)) stop("local_energy not available at requested scale")
  vals <- pmax(fit$observed$statistic[[idx]], 0)
  if (sum(vals) > 0) vals <- vals / sum(vals)
  if (!is.null(mass)) {
    ord <- order(vals, decreasing = TRUE)
    keep <- rep(FALSE, length(vals))
    keep[ord[cumsum(vals[ord]) <= mass]] <- TRUE
    if (!any(keep) && length(ord)) keep[ord[1L]] <- TRUE
    vals[!keep] <- 0
  }
  full <- rep(NA_real_, prod(fit$dims))
  full[fit$mask_idx] <- vals
  .array_to_vol(full, fit$dims, fit$template)
}
