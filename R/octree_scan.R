#' Octree Max-Scan with permutation FWER (U0-only)
#'
#' Computes a multiscale dyadic-cube scan statistic using the variance-stabilized
#' prior-weighted mean
#' \deqn{U_0(R) = \sum_{v \in R}\pi(v)Z(v)/\sqrt{\sum_{v \in R}\pi(v)^2}}
#' over all dyadic cubes in an implicit octree, and calibrates a global
#' threshold by a maxT permutation (sign-flip) procedure.
#'
#' @param z_vol A 3D volume-like object with finite values in-brain.
#' @param prior_vol Optional prior weight volume aligned to \code{z_vol}.
#'   If NULL, uniform weights are used.
#' @param mask Optional logical mask (same shape as \code{z_vol}). If NULL,
#'   uses all finite voxels in \code{z_vol}.
#' @param alpha Significance level for strong FWER control (default 0.05).
#' @param n_perm Number of permutations/sign-flips (default 1000).
#' @param null Null calibration method.
#'   \describe{
#'     \item{\code{"mc_fwhm"}}{Correlated Monte Carlo null via Gaussian smoothing
#'       of white noise with smoothness set by \code{fwhm_vox} (or \code{fwhm_mm}
#'       for \code{NeuroVol}s). If both are NULL, an effective \code{fwhm_vox} is
#'       estimated from the statistic map using \code{\link{estimate_fwhm_vox}}.}
#'     \item{\code{"mc_acf"}}{Correlated Monte Carlo null using a simple mixed-ACF
#'       model (Gaussian + long-tailed exponential) with user-supplied parameters
#'       \code{acf_params}.}
#'     \item{\code{"signflip_voxel"}}{Independent voxelwise sign-flips of the
#'       group map (fast but generally invalid for regional/multiscale statistics;
#'       mainly useful for debugging).}
#'   }
#' @param null_fun Optional function \code{function(b)} returning a length
#'   \code{n_mask} Z-vector for permutation \code{b}. Use this to supply
#'   subject-level randomization that preserves spatial covariance.
#' @param fwhm_vox For \code{null="mc_fwhm"}, Gaussian FWHM in voxel units
#'   (length 1 or 3).
#' @param fwhm_mm For \code{null="mc_fwhm"}, Gaussian FWHM in mm (length 1 or 3).
#'   Requires \code{z_vol} to be a \code{neuroim2::NeuroVol} so voxel spacing is
#'   available.
#' @param acf_params For \code{null="mc_acf"}, a list with components:
#'   \describe{
#'     \item{a}{Mixing weight in (0, 1) for the Gaussian component.}
#'     \item{fwhm_vox}{Gaussian ACF FWHM in voxel units (length 1 or 3).}
#'     \item{lambda_vox}{Exponential kernel scale in voxel units (length 1 or 3).}
#'   }
#' @param seed Optional RNG seed for reproducibility.
#' @param two_sided Logical; if TRUE, uses \code{abs(z_vol)} prior to sign-flips.
#' @param prior_eta Mixing weight in \[0, 1\] to shrink the prior toward uniform
#'   mass (default 0.9).
#'
#' @return A list with components:
#'   \describe{
#'     \item{u}{Global maxT threshold at level \code{alpha}.}
#'     \item{M_obs}{Observed max scan statistic.}
#'     \item{M_null}{Numeric vector of permutation maxima.}
#'     \item{nodes}{Data frame of all non-empty dyadic nodes with columns
#'       \code{level,i,j,k,score,scale,x0,x1,y0,y1,z0,z1}.}
#'     \item{sig_nodes}{Subset of \code{nodes} with \code{score > u}.}
#'     \item{params}{List of run parameters.}
#'   }
#'
#' @details
#' This implements the simplest global calibration: a single-step maxT
#' permutation threshold over the full dyadic-cube family.
#'
#' If only a single group-level map is available, \code{null="mc_fwhm"} is a
#' safer default than voxelwise sign-flips because it preserves an assumed
#' spatial autocorrelation model. Exact resampling-based inference typically
#' requires subject-level or residual maps; supply \code{null_fun} to use such
#' randomizations.
#'
#' For hierarchical step-down over sibling families during descent, see
#' \code{\link{octree_scan_stepdown}}.
#'
#' @examples
#' # Map-only (model-based) calibration using a Gaussian ACF null:
#' z <- array(rnorm(20 * 20 * 20), dim = c(20, 20, 20))
#' res <- octree_scan_fwer(z, n_perm = 200, null = "mc_fwhm")
#'
#' # User-supplied mixed ACF parameters (Gaussian + long tail):
#' acf <- list(a = 0.7, fwhm_vox = c(4, 4, 4), lambda_vox = c(2, 2, 2))
#' res2 <- octree_scan_fwer(z, n_perm = 200, null = "mc_acf", acf_params = acf)
#'
#' # Subject-level sign-flips (recommended when subject maps are available):
#' subj <- replicate(8, array(rnorm(20 * 20 * 20), dim = c(20, 20, 20)), simplify = FALSE)
#' sf <- make_null_fun_subject_signflip(subj)
#' res3 <- octree_scan_fwer(sf$z_vol, mask = sf$mask, n_perm = 200, null_fun = sf$null_fun)
#'
#' @export
octree_scan_fwer <- function(z_vol,
                             prior_vol = NULL,
                             mask = NULL,
                             alpha = 0.05,
                             n_perm = 1000,
                             null = c("signflip_voxel", "mc_fwhm", "mc_acf"),
                             null_fun = NULL,
                             fwhm_vox = NULL,
                             fwhm_mm = NULL,
                             acf_params = NULL,
                             seed = NULL,
                             two_sided = FALSE,
                             prior_eta = 0.9) {
  null <- match.arg(null)
  if (!is.null(seed)) set.seed(seed)

  dims <- dim(z_vol)
  if (length(dims) != 3) stop("z_vol must be a 3D volume-like object")

  if (is.null(mask)) {
    mask <- is.finite(as.numeric(z_vol))
  } else {
    mask <- as.logical(mask)
  }

  mask_idx <- which(mask)
  n_mask <- length(mask_idx)
  if (n_mask == 0) stop("No valid voxels in mask")

  z_vec <- as.numeric(z_vol)[mask_idx]
  if (two_sided) z_vec <- abs(z_vec)

  if (is.null(prior_vol)) {
    pi_vec <- rep(1.0, n_mask)
  } else {
    pi_vec <- as.numeric(prior_vol)[mask_idx]
  }
  pi_vec <- .prep_prior(pi_vec, eta = prior_eta)

  grid <- .index_to_grid(mask_idx, dims)
  x0 <- as.integer(grid[, 1] - 1L)
  y0 <- as.integer(grid[, 2] - 1L)
  z0 <- as.integer(grid[, 3] - 1L)

  side <- .next_pow2(max(dims))
  den_xptr <- build_den_pyramid_xptr_cpp(x0, y0, z0, pi_vec, side)
  mask_full_idx0 <- .full_index0(x0, y0, z0, side)

  obs <- pyramid_scores_u0_cpp(z_vec, pi_vec, x0, y0, z0, den_xptr, do_abs = FALSE)
  M_obs <- obs$max_score

  if (is.null(null_fun)) {
    if (null == "signflip_voxel") {
      warning(
        "null='signflip_voxel' flips signs independently per voxel and does not preserve spatial covariance; ",
        "this can be anti-conservative for regional/multiscale statistics. Prefer null_fun (subject-level resampling) ",
        "or null='mc_fwhm' with an appropriate fwhm_vox.",
        call. = FALSE
      )
      null_fun <- make_null_fun_voxel_signflip(z_vec)
    } else if (null == "mc_fwhm") {
      if (!is.null(fwhm_mm)) {
        if (inherits(z_vol, "NeuroVol")) {
          spacing_mm <- neuroim2::spacing(z_vol)
          fwhm_vox <- .fwhm_mm_to_vox(fwhm_mm, spacing_mm)
        } else {
          stop("fwhm_mm requires z_vol to be a neuroim2::NeuroVol (for voxel spacing)")
        }
      }
      if (is.null(fwhm_vox)) {
        fwhm_vox <- estimate_fwhm_vox(z_vol, mask = mask, isotropic = FALSE)
      }
      null_fun <- make_null_fun_mc_fwhm(x0, y0, z0, mask_full_idx0, side, fwhm_vox = fwhm_vox)
    } else {
      null_fun <- make_null_fun_mc_acf(x0, y0, z0, mask_full_idx0, side, acf_params = acf_params)
    }
  }

  M_null <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    Zb <- null_fun(b)
    M_null[b] <- pyramid_max_u0_cpp(Zb, pi_vec, x0, y0, z0, den_xptr, do_abs = FALSE)
  }

  u <- unname(stats::quantile(M_null, probs = 1 - alpha, names = FALSE))

  nodes <- obs$nodes
  nodes$scale <- 2^nodes$level
  bbox <- .nodes_to_bbox(nodes, dims)
  nodes <- cbind(nodes, bbox)

  sig_nodes <- nodes[nodes$score > u, , drop = FALSE]

  list(
    u = u,
    M_obs = M_obs,
    M_null = M_null,
    nodes = nodes,
    sig_nodes = sig_nodes,
    params = list(
      alpha = alpha,
      n_perm = n_perm,
      two_sided = two_sided,
      side = side,
      null = null,
      fwhm_vox = fwhm_vox,
      fwhm_mm = fwhm_mm
    )
  )
}

#' @keywords internal
.next_pow2 <- function(n) {
  if (n <= 1) return(1L)
  2L^(ceiling(log2(as.numeric(n))))
}

#' @keywords internal
.nodes_to_bbox <- function(nodes, dims) {
  s <- as.integer(nodes$scale)
  x0 <- (nodes$i - 1L) * s + 1L
  y0 <- (nodes$j - 1L) * s + 1L
  z0 <- (nodes$k - 1L) * s + 1L

  x1 <- pmin.int(x0 + s - 1L, dims[1])
  y1 <- pmin.int(y0 + s - 1L, dims[2])
  z1 <- pmin.int(z0 + s - 1L, dims[3])

  x0 <- pmax.int(x0, 1L)
  y0 <- pmax.int(y0, 1L)
  z0 <- pmax.int(z0, 1L)

  data.frame(x0 = x0, x1 = x1, y0 = y0, y1 = y1, z0 = z0, z1 = z1)
}
