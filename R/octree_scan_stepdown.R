#' Octree scan with hierarchical step-down inference (U0-only)
#'
#' Performs a multiscale dyadic-cube search using the variance-stabilized
#' prior-weighted mean statistic \eqn{U_0(R)} and applies a Westfall-Young
#' step-down procedure *within each sibling family* during octree descent
#' (alpha-spending across the tree).
#'
#' Compared to \code{\link{octree_scan_fwer}}, this avoids a single global maxT
#' threshold over all dyadic cubes. Instead, it leverages the nesting of the
#' tree and tests only small sibling families (<= 8) at each split.
#'
#' Note: this is not the same as a global Westfall-Young step-down over the full
#' set of dyadic cubes (which would require the joint null distribution over
#' all hypotheses). Use \code{\link{octree_scan_fwer}} for a single-step maxT
#' procedure over the full dyadic-cube family.
#'
#' @param z_vol A 3D volume-like object containing Z-equivalent values.
#' @param prior_vol Optional prior weight volume aligned to \code{z_vol}.
#'   If NULL, uniform weights are used.
#' @param mask Optional logical mask (same shape as \code{z_vol}). If NULL,
#'   uses all finite voxels in \code{z_vol}.
#' @param alpha FWER target (default 0.05).
#' @param n_perm Number of sign-flip permutations per tested family (default 1000).
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
#' @param gamma Fraction of the local alpha budget used to test the current
#'   sibling family (default 0.5). Remaining budget is split among rejected
#'   children by prior mass.
#' @param min_voxels Minimum child size to consider for descent (default 8).
#' @param min_alpha Minimum alpha budget to continue descent (default 1e-6).
#' @param report How to select a non-redundant set of regions from the discovered
#'   significant nodes: \code{"coarsest"}, \code{"finest"}, or \code{"greedy"}.
#' @param seed Optional RNG seed.
#' @param two_sided Logical; if TRUE, uses \code{abs(z_vol)} for the observed map
#'   and permutations (two-sided via absolute Z).
#' @param prior_eta Mixing weight in \[0, 1\] to shrink the prior toward uniform
#'   mass (default 0.9).
#'
#' @return A list with components:
#'   \describe{
#'     \item{hits}{Data frame of all discovered significant nodes (may include nested nodes).}
#'     \item{selected}{Data frame of selected nodes per \code{report}.}
#'     \item{sig_mask}{Logical mask volume for \code{selected} nodes.}
#'     \item{sig_indices}{1-based linear indices into \code{z_vol} for \code{sig_mask}.}
#'     \item{params}{Run parameters.}
#'   }
#'
#' @export
octree_scan_stepdown <- function(z_vol,
                                 prior_vol = NULL,
                                 mask = NULL,
                                 alpha = 0.05,
                                 n_perm = 1000,
                                 null = c("signflip_voxel", "mc_fwhm", "mc_acf"),
                                 null_fun = NULL,
                                 fwhm_vox = NULL,
                                 fwhm_mm = NULL,
                                 acf_params = NULL,
                                 gamma = 0.5,
                                 min_voxels = 8,
                                 min_alpha = 1e-6,
                                 report = c("coarsest", "finest", "greedy"),
                                 seed = NULL,
                                 two_sided = FALSE,
                                 prior_eta = 0.9) {
  report <- match.arg(report)
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
  x <- as.integer(grid[, 1])
  y <- as.integer(grid[, 2])
  z <- as.integer(grid[, 3])

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
      x0 <- as.integer(x - 1L)
      y0 <- as.integer(y - 1L)
      z0 <- as.integer(z - 1L)
      side <- .next_pow2(max(dim(z_vol)))
      mask_full_idx0 <- .full_index0(x0, y0, z0, side)
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
      x0 <- as.integer(x - 1L)
      y0 <- as.integer(y - 1L)
      z0 <- as.integer(z - 1L)
      side <- .next_pow2(max(dim(z_vol)))
      mask_full_idx0 <- .full_index0(x0, y0, z0, side)
      null_fun <- make_null_fun_mc_acf(x0, y0, z0, mask_full_idx0, side, acf_params = acf_params)
    }
  }

  perm_fun <- function(b) null_fun(b)

  id_env <- new.env(parent = emptyenv())
  id_env$next_id <- 1L
  .next_id <- function() {
    out <- id_env$next_id
    id_env$next_id <- out + 1L
    out
  }

  idx_root <- seq_len(n_mask)
  bbox_root <- bbox_from_indices_cpp(idx_root, x, y, z)

  hits_list <- .octree_descend_u0_stepdown(
    idx_parent = idx_root,
    bbox = bbox_root,
    alpha_budget = alpha,
    z_vec = z_vec,
    pi_vec = pi_vec,
    x = x, y = y, z = z,
    perm_fun = perm_fun,
    n_perm = n_perm,
    gamma = gamma,
    min_alpha = min_alpha,
    min_voxels = min_voxels,
    parent_id = NA_integer_,
    depth = 0L,
    next_id = .next_id
  )

  hits <- .hits_to_table(hits_list)
  selected <- .select_nodes(hits, report = report)

  sig_idx_mask <- integer(0)
  if (nrow(selected) > 0) {
    sig_idx_mask <- unique(unlist(selected$indices, use.names = FALSE))
  }
  sig_idx_vol <- mask_idx[sig_idx_mask]

  sig_any <- rep(FALSE, prod(dims))
  sig_any[sig_idx_vol] <- TRUE
  sig_mask <- .array_to_vol(sig_any, dims, z_vol)

  list(
    hits = hits,
    selected = selected,
    sig_mask = sig_mask,
    sig_indices = sig_idx_vol,
    params = list(
      alpha = alpha,
      n_perm = n_perm,
      null = null,
      fwhm_vox = fwhm_vox,
      fwhm_mm = fwhm_mm,
      gamma = gamma,
      min_voxels = min_voxels,
      min_alpha = min_alpha,
      report = report,
      two_sided = two_sided
    )
  )
}

#' @keywords internal
.octree_descend_u0_stepdown <- function(idx_parent, bbox, alpha_budget,
                                        z_vec, pi_vec, x, y, z,
                                        perm_fun, n_perm,
                                        gamma, min_alpha, min_voxels,
                                        parent_id, depth, next_id) {
  hits <- list()

  if (alpha_budget < min_alpha || length(idx_parent) < min_voxels) {
    return(hits)
  }

  split <- octree_split_info_cpp(idx_parent, bbox, x, y, z, pi_vec)
  if (split$m == 0) return(hits)

  alpha_test <- gamma * alpha_budget
  alpha_desc <- alpha_budget - alpha_test

  t_obs <- score_children_onepass_cpp(
    z_vec, idx_parent, split$w_parent, split$child_id,
    split$log_den1, split$inv_sqrt_den2,
    kappa_grid_pos = numeric(0),
    do_abs = FALSE
  )

  t_null <- matrix(NA_real_, nrow = n_perm, ncol = split$m)
  for (b in seq_len(n_perm)) {
    z_perm <- perm_fun(b)
    t_null[b, ] <- score_children_onepass_cpp(
      z_perm, idx_parent, split$w_parent, split$child_id,
      split$log_den1, split$inv_sqrt_den2,
      kappa_grid_pos = numeric(0),
      do_abs = FALSE
    )
  }

  sd <- wy_stepdown(t_obs, t_null, alpha = alpha_test)
  rej <- which(sd$rejected)
  if (length(rej) == 0) return(hits)

  res <- split_indices_bbox_mass_cpp(idx_parent, split$child_id, split$m, x, y, z, pi_vec)

  child_node_ids <- integer(length(rej))

  for (k in seq_along(rej)) {
    i <- rej[k]
    node_id <- next_id()
    child_node_ids[k] <- node_id

    idx_child <- res$idx[[i]]
    bbox_child <- as.integer(res$bbox[i, ])

    hits[[length(hits) + 1]] <- list(
      node_id = node_id,
      parent_id = parent_id,
      depth = depth + 1L,
      indices = idx_child,
      bbox = bbox_child,
      score = t_obs[i],
      p_adj = sd$p_adj[i],
      alpha_test = alpha_test,
      n_voxels = length(idx_child),
      mass = res$den1[i]
    )
  }

  if (alpha_desc < min_alpha) return(hits)

  child_mass <- res$den1[rej]
  w <- child_mass / sum(child_mass)

  for (k in seq_along(rej)) {
    i <- rej[k]
    idx_child <- res$idx[[i]]
    if (length(idx_child) < min_voxels) next

    bbox_child <- as.integer(res$bbox[i, ])

    child_hits <- .octree_descend_u0_stepdown(
      idx_parent = idx_child,
      bbox = bbox_child,
      alpha_budget = alpha_desc * w[k],
      z_vec = z_vec,
      pi_vec = pi_vec,
      x = x, y = y, z = z,
      perm_fun = perm_fun,
      n_perm = n_perm,
      gamma = gamma,
      min_alpha = min_alpha,
      min_voxels = min_voxels,
      parent_id = child_node_ids[k],
      depth = depth + 1L,
      next_id = next_id
    )
    hits <- c(hits, child_hits)
  }

  hits
}

#' @keywords internal
.hits_to_table <- function(hits_list) {
  if (length(hits_list) == 0) {
    return(data.frame(
      node_id = integer(0),
      parent_id = integer(0),
      depth = integer(0),
      score = numeric(0),
      p_adj = numeric(0),
      alpha_test = numeric(0),
      n_voxels = integer(0),
      mass = numeric(0),
      x0 = integer(0), x1 = integer(0),
      y0 = integer(0), y1 = integer(0),
      z0 = integer(0), z1 = integer(0),
      indices = I(list())
    ))
  }

  node_id <- vapply(hits_list, `[[`, integer(1), "node_id")
  parent_id <- vapply(hits_list, function(x) {
    pid <- x$parent_id
    if (is.na(pid)) NA_integer_ else as.integer(pid)
  }, integer(1))
  depth <- vapply(hits_list, `[[`, integer(1), "depth")
  score <- vapply(hits_list, `[[`, numeric(1), "score")
  p_adj <- vapply(hits_list, `[[`, numeric(1), "p_adj")
  alpha_test <- vapply(hits_list, `[[`, numeric(1), "alpha_test")
  n_voxels <- vapply(hits_list, `[[`, integer(1), "n_voxels")
  mass <- vapply(hits_list, `[[`, numeric(1), "mass")
  bbox_mat <- do.call(rbind, lapply(hits_list, `[[`, "bbox"))
  indices <- lapply(hits_list, `[[`, "indices")

  data.frame(
    node_id = node_id,
    parent_id = parent_id,
    depth = depth,
    score = score,
    p_adj = p_adj,
    alpha_test = alpha_test,
    n_voxels = n_voxels,
    mass = mass,
    x0 = bbox_mat[, 1], x1 = bbox_mat[, 2],
    y0 = bbox_mat[, 3], y1 = bbox_mat[, 4],
    z0 = bbox_mat[, 5], z1 = bbox_mat[, 6],
    indices = I(indices)
  )
}

#' @keywords internal
.select_nodes <- function(hits, report = c("coarsest", "finest", "greedy")) {
  report <- match.arg(report)
  if (nrow(hits) == 0) return(hits)

  ids <- hits$node_id
  parent_ids <- hits$parent_id

  if (report == "coarsest") {
    keep <- is.na(parent_ids) | !(parent_ids %in% ids)
    return(hits[keep, , drop = FALSE])
  }

  if (report == "finest") {
    has_sig_child <- ids %in% parent_ids
    return(hits[!has_sig_child, , drop = FALSE])
  }

  ord <- order(hits$score, decreasing = TRUE)
  kept <- integer(0)

  bbox_overlap <- function(a, b) {
    !(a["x1"] < b["x0"] || b["x1"] < a["x0"] ||
        a["y1"] < b["y0"] || b["y1"] < a["y0"] ||
        a["z1"] < b["z0"] || b["z1"] < a["z0"])
  }

  for (ii in ord) {
    if (length(kept) == 0) {
      kept <- c(kept, ii)
      next
    }

    a <- hits[ii, c("x0", "x1", "y0", "y1", "z0", "z1")]
    ok <- TRUE
    for (jj in kept) {
      b <- hits[jj, c("x0", "x1", "y0", "y1", "z0", "z1")]
      if (bbox_overlap(a, b)) {
        ok <- FALSE
        break
      }
    }
    if (ok) kept <- c(kept, ii)
  }

  hits[kept, , drop = FALSE]
}
