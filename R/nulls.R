#' @keywords internal
.fwhm_to_sigma_vox <- function(fwhm_vox) {
  # FWHM = sqrt(8 log 2) * sigma
  fwhm_vox / sqrt(8 * log(2))
}

#' Estimate Gaussian ACF FWHM from a single map (map-only)
#'
#' Map-only smoothness estimator based on robust sign correlations at short lags.
#' This assumes approximate stationarity and that signal occupies a relatively
#' small fraction of the mask.
#'
#' @param z_vol 3D volume-like object.
#' @param mask Optional logical mask aligned to \code{z_vol}. If NULL, uses
#'   \code{is.finite(as.numeric(z_vol))}.
#' @param clip Clip value applied to \code{z_vol} prior to taking \code{sign},
#'   to reduce sensitivity to large outliers (default 3).
#' @param isotropic Logical; if TRUE, returns a single averaged FWHM replicated
#'   over x/y/z.
#' @param z_trim If not NULL, restricts estimation to voxels with
#'   \code{|Z| <= z_trim} (default 1.5) as a heuristic to reduce contamination
#'   from true signal when only a single statistic map is available.
#'
#' @return Numeric vector length 3: estimated \code{fwhm_vox} in voxel units.
#'
#' @details
#' This estimates an effective stationary smoothness parameter from the
#' statistic field itself (not from residuals). It is intended for
#' model-based Monte Carlo calibration when subject-level/residual resampling
#' is not available.
#'
#' @export
estimate_fwhm_vox <- function(z_vol,
                              mask = NULL,
                              clip = 3,
                              isotropic = FALSE,
                              z_trim = 1.5) {
  dims <- dim(z_vol)
  if (length(dims) != 3) stop("z_vol must be a 3D volume-like object")

  z_full <- as.numeric(z_vol)
  if (length(z_full) != prod(dims)) stop("z_vol length does not match dim(z_vol)")

  if (is.null(mask)) {
    mask_full <- is.finite(z_full)
  } else {
    mask_full <- as.logical(mask)
  }
  if (length(mask_full) != length(z_full)) stop("mask must match z_vol dimensions")

  if (!is.null(z_trim)) {
    if (!is.finite(z_trim) || z_trim <= 0) stop("z_trim must be positive or NULL")
    # Heuristic: treat likely-null voxels as those with small |Z|.
    mask_full <- mask_full & is.finite(z_full) & (abs(z_full) <= z_trim)
  }

  dims_i <- as.integer(dims)

  lags <- list(
    c(1L, 0L, 0L),
    c(0L, 1L, 0L),
    c(0L, 0L, 1L)
  )

  rho_hat <- numeric(3)
  for (i in 1:3) {
    out <- acf_sign_lag_mean_cpp(
      z_full = z_full,
      mask_full = mask_full,
      dims = dims_i,
      lag_xyz = as.integer(lags[[i]]),
      clip = clip
    )
    m <- out$mean
    if (!is.finite(m)) {
      rho_hat[i] <- NA_real_
      next
    }
    # For bivariate normal: E[sign(X)sign(Y)] = (2/pi) asin(rho)
    rho <- sin((pi / 2) * m)
    rho_hat[i] <- rho
  }

  # Convert lag-1 correlation to Gaussian ACF sigma via rho = exp(-1^2 / (2 sigma^2))
  sigma <- rep(NA_real_, 3)
  for (i in 1:3) {
    rho <- rho_hat[i]
    if (!is.finite(rho) || rho <= 0 || rho >= 1) next
    sigma[i] <- sqrt(-1 / (2 * log(rho)))
  }

  fwhm <- sqrt(8 * log(2)) * sigma

  if (isTRUE(isotropic)) {
    f0 <- mean(fwhm, na.rm = TRUE)
    if (!is.finite(f0)) stop("Could not estimate FWHM from map")
    return(rep(f0, 3))
  }

  if (any(!is.finite(fwhm))) stop("Could not estimate FWHM from map (non-finite)")
  fwhm
}

#' @keywords internal
.fwhm_mm_to_vox <- function(fwhm_mm, spacing_mm) {
  if (length(fwhm_mm) == 1) fwhm_mm <- rep(fwhm_mm, 3)
  if (length(fwhm_mm) != 3) stop("fwhm_mm must be length 1 or 3")
  if (length(spacing_mm) == 1) spacing_mm <- rep(spacing_mm, 3)
  if (length(spacing_mm) != 3) stop("spacing_mm must be length 1 or 3")
  as.numeric(fwhm_mm) / as.numeric(spacing_mm)
}

#' @keywords internal
.prep_fwhm_vox <- function(fwhm_vox) {
  if (is.null(fwhm_vox)) stop("fwhm_vox must be provided for null='mc_fwhm'")
  if (length(fwhm_vox) == 1) fwhm_vox <- rep(fwhm_vox, 3)
  if (length(fwhm_vox) != 3) stop("fwhm_vox must be length 1 or 3")
  if (any(!is.finite(fwhm_vox)) || any(fwhm_vox < 0)) stop("fwhm_vox must be finite and >= 0")
  as.numeric(fwhm_vox)
}

#' @keywords internal
.full_index0 <- function(x0, y0, z0, side) {
  x0 + side * (y0 + side * z0) + 1L
}

#' @keywords internal
.standardize_mask <- function(z_mask) {
  mu <- mean(z_mask)
  s <- stats::sd(z_mask)
  if (!is.finite(s) || s <= 0) {
    return(z_mask - mu)
  }
  (z_mask - mu) / s
}

#' @keywords internal
make_null_fun_mc_fwhm <- function(x0, y0, z0, mask_full_idx0, side, fwhm_vox) {
  fwhm_vox <- .prep_fwhm_vox(fwhm_vox)
  # Interpret fwhm_vox as the target ACF (field) smoothness; smoothing white noise
  # by a Gaussian kernel yields an ACF with sigma_field = sqrt(2) * sigma_kernel.
  sigma_field <- .fwhm_to_sigma_vox(fwhm_vox)
  sigma <- sigma_field / sqrt(2)
  dims_full <- as.integer(c(side, side, side))

  function(b) {
    eps <- stats::rnorm(side * side * side)
    # neuroim2 also provides gaussian_blur(), but its current implementation
    # computes a full 3D neighborhood kernel per voxel (O(N * window^3)),
    # which is typically too slow for large Monte Carlo loops here.
    sm <- gaussian_smooth3d_cpp(eps, dims_full, sigma_xyz = sigma)
    z_mask <- sm[mask_full_idx0]
    .standardize_mask(z_mask)
  }
}

#' @keywords internal
make_null_fun_mc_acf <- function(x0, y0, z0, mask_full_idx0, side, acf_params) {
  if (is.null(acf_params)) stop("acf_params must be provided for null='mc_acf'")
  a <- acf_params$a
  fwhm_vox <- acf_params$fwhm_vox
  lambda_vox <- acf_params$lambda_vox
  if (!is.numeric(a) || length(a) != 1 || !is.finite(a) || a <= 0 || a >= 1) {
    stop("acf_params$a must be a scalar in (0, 1)")
  }
  fwhm_vox <- .prep_fwhm_vox(fwhm_vox)
  if (length(lambda_vox) == 1) lambda_vox <- rep(lambda_vox, 3)
  if (length(lambda_vox) != 3) stop("acf_params$lambda_vox must be length 1 or 3")
  if (any(!is.finite(lambda_vox)) || any(lambda_vox <= 0)) stop("acf_params$lambda_vox must be finite and > 0")
  lambda_vox <- as.numeric(lambda_vox)

  sigma_field <- .fwhm_to_sigma_vox(fwhm_vox)
  sigma_kernel <- sigma_field / sqrt(2)
  dims_full <- as.integer(c(side, side, side))

  function(b) {
    eps1 <- stats::rnorm(side * side * side)
    eps2 <- stats::rnorm(side * side * side)

    g <- gaussian_smooth3d_cpp(eps1, dims_full, sigma_xyz = sigma_kernel)
    e <- exp_smooth3d_cpp(eps2, dims_full, lambda_xyz = lambda_vox)

    sm <- sqrt(a) * g + sqrt(1 - a) * e
    z_mask <- sm[mask_full_idx0]
    .standardize_mask(z_mask)
  }
}

#' @keywords internal
make_null_fun_voxel_signflip <- function(z_vec) {
  n <- length(z_vec)
  function(b) {
    signs <- sample(c(-1, 1), n, replace = TRUE)
    signs * z_vec
  }
}
