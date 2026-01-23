#' Hierarchical Testing with Alpha-Spending
#'
#' Implements top-down hierarchical testing with alpha-spending for
#' multi-scale inference.
#'
#' @name hierarchical
NULL

#' Test a tree of regions with alpha-spending
#'
#' Recursive function that tests children only if the parent is rejected,
#' with alpha budget allocated proportionally to prior mass.
#'
#' @param idx_parent Integer vector of parent mask-space indices
#' @param bbox Integer vector of length 6 for parent bounding box
#' @param alpha_budget Alpha budget for this subtree
#' @param z_vec Numeric vector of Z-scores (mask-space)
#' @param pi_vec Numeric vector of prior weights (mask-space)
#' @param x,y,z Integer vectors of coordinates (mask-space)
#' @param perm_fun Function that takes permutation index and returns Z vector
#' @param n_perm Number of permutations
#' @param kappa_grid Numeric vector of positive kappa values
#' @param gamma Fraction of alpha for testing at each level (default 0.5)
#' @param min_alpha Minimum alpha to continue testing (default 1e-6)
#' @param min_voxels Minimum region size to test (default 8)
#'
#' @return List of discovered regions with scores and p-values
#'
#' @keywords internal
hier_descend <- function(idx_parent, bbox, alpha_budget,
                         z_vec, pi_vec, x, y, z,
                         perm_fun, n_perm, kappa_grid,
                         gamma = 0.5, min_alpha = 1e-6,
                         min_voxels = 8,
                         parallel = FALSE) {

  hits <- list()

  if (alpha_budget < min_alpha || length(idx_parent) < min_voxels) {
    return(hits)
  }

  kappa_pos <- kappa_grid[kappa_grid > 0]

  # Build split info
  split <- octree_split_info_cpp(idx_parent, bbox, x, y, z, pi_vec)

  if (split$m == 0) {
    return(hits)
  }

  # Allocate alpha
  alpha_test <- gamma * alpha_budget
  alpha_desc <- alpha_budget - alpha_test

  # Observed child scores
  t_obs <- score_children_onepass_cpp(
    z_vec, idx_parent, split$w_parent, split$child_id,
    split$log_den1, split$inv_sqrt_den2, kappa_pos, do_abs = FALSE
  )

  # Null child scores
  if (parallel) {
    t_null_list <- future.apply::future_lapply(seq_len(n_perm), function(b) {
      z_perm <- perm_fun(b)
      score_children_onepass_cpp(
        z_perm, idx_parent, split$w_parent, split$child_id,
        split$log_den1, split$inv_sqrt_den2, kappa_pos, do_abs = FALSE
      )
    }, future.seed = TRUE)
    t_null <- do.call(rbind, t_null_list)
  } else {
    t_null <- matrix(NA_real_, nrow = n_perm, ncol = split$m)
    for (b in seq_len(n_perm)) {
      z_perm <- perm_fun(b)
      t_null[b, ] <- score_children_onepass_cpp(
        z_perm, idx_parent, split$w_parent, split$child_id,
        split$log_den1, split$inv_sqrt_den2, kappa_pos, do_abs = FALSE
      )
    }
  }

  # Westfall-Young step-down
  sd <- wy_stepdown(t_obs, t_null, alpha = alpha_test)
  rej <- which(sd$rejected)

  if (length(rej) == 0) {
    return(hits)
  }

  # Get child indices and bboxes
  res <- split_indices_bbox_mass_cpp(idx_parent, split$child_id, split$m,
                                     x, y, z, pi_vec)

  # Record rejected children
  for (i in rej) {
    hits[[length(hits) + 1]] <- list(
      indices = res$idx[[i]],
      score = t_obs[i],
      p_adj = sd$p_adj[i],
      alpha_test = alpha_test,
      level = "child"
    )
  }

  # Allocate descendant alpha by prior mass
  child_mass <- res$den1[rej]
  w <- child_mass / sum(child_mass)

  # Recurse into rejected children
  for (k in seq_along(rej)) {
    i <- rej[k]
    idx_child <- res$idx[[i]]
    bbox_child <- res$bbox[i, ]

    if (length(idx_child) >= min_voxels) {
      child_hits <- hier_descend(
        idx_child, bbox_child,
        alpha_budget = alpha_desc * w[k],
        z_vec, pi_vec, x, y, z,
        perm_fun, n_perm, kappa_grid,
        gamma, min_alpha, min_voxels,
        parallel = parallel
      )
      hits <- c(hits, child_hits)
    }
  }

  hits
}

#' Alpha-spending hierarchical descent without step-down
#'
#' Tests each node against its own permutation distribution and only
#' descends when the node is significant at its allocated alpha.
#'
#' @param idx_parent Integer vector of parent mask-space indices
#' @param bbox Integer vector of length 6 for parent bounding box
#' @param alpha_budget Alpha budget for this subtree
#' @param z_vec Numeric vector of Z-scores (mask-space)
#' @param pi_vec Numeric vector of prior weights (mask-space)
#' @param x,y,z Integer vectors of coordinates (mask-space)
#' @param perm_fun Function that takes permutation index and returns Z vector
#' @param n_perm Number of permutations
#' @param kappa_grid Numeric vector of positive kappa values
#' @param gamma Fraction of alpha for testing at each level (default 0.5)
#' @param min_alpha Minimum alpha to continue testing (default 1e-6)
#' @param min_voxels Minimum region size to test (default 8)
#'
#' @return List of discovered regions with scores and p-values
#'
#' @keywords internal
hier_descend_alpha <- function(idx_parent, bbox, alpha_budget,
                               z_vec, pi_vec, x, y, z,
                               perm_fun, n_perm, kappa_grid,
                               gamma = 0.5, min_alpha = 1e-6,
                               min_voxels = 8,
                               parallel = FALSE) {

  hits <- list()

  if (alpha_budget < min_alpha || length(idx_parent) < min_voxels) {
    return(hits)
  }

  obs <- score_set_omnibus(idx_parent, z_vec, pi_vec, kappa_grid)
  s_obs <- obs$T_omnibus

  ge <- 0
  if (parallel) {
    perm_scores <- future.apply::future_lapply(seq_len(n_perm), function(b) {
      z_perm <- perm_fun(b)
      score_set_omnibus(idx_parent, z_perm, pi_vec, kappa_grid)$T_omnibus
    }, future.seed = TRUE)
    ge <- sum(unlist(perm_scores) >= s_obs)
  } else {
    for (b in seq_len(n_perm)) {
      z_perm <- perm_fun(b)
      s_perm <- score_set_omnibus(idx_parent, z_perm, pi_vec, kappa_grid)$T_omnibus
      if (s_perm >= s_obs) ge <- ge + 1
    }
  }

  p_val <- (ge + 1) / (n_perm + 1)
  alpha_test <- gamma * alpha_budget

  if (p_val > alpha_test) {
    return(hits)
  }

  hits[[length(hits) + 1]] <- list(
    indices = idx_parent,
    score = s_obs,
    p_adj = p_val,
    alpha_test = alpha_test,
    level = "node",
    kappa_best = obs$kappa_best
  )

  alpha_desc <- alpha_budget - alpha_test

  split <- octree_split_info_cpp(idx_parent, bbox, x, y, z, pi_vec)
  if (split$m == 0 || alpha_desc < min_alpha) {
    return(hits)
  }

  res <- split_indices_bbox_mass_cpp(idx_parent, split$child_id, split$m,
                                     x, y, z, pi_vec)
  child_mass <- res$den1
  child_mass <- child_mass / sum(child_mass)

  for (k in seq_len(split$m)) {
    idx_child <- res$idx[[k]]
    if (length(idx_child) < min_voxels) next
    bbox_child <- res$bbox[k, ]
    child_hits <- hier_descend_alpha(
      idx_child, bbox_child,
      alpha_budget = alpha_desc * child_mass[k],
      z_vec, pi_vec, x, y, z,
      perm_fun, n_perm, kappa_grid,
      gamma, min_alpha, min_voxels,
      parallel = parallel
    )
    hits <- c(hits, child_hits)
  }

  hits
}

#' Test at parcel level with step-down
#'
#' Tests all parcels simultaneously using Westfall-Young step-down,
#' then returns significant parcels for further descent.
#'
#' @param z_vec Numeric vector of Z-scores (mask-space)
#' @param pi_vec Numeric vector of prior weights (mask-space)
#' @param lab_vec Integer vector of parcel labels (mask-space)
#' @param n_perm Number of permutations
#' @param kappa_grid Numeric vector of positive kappa values
#' @param alpha Significance level
#' @param seed Random seed (optional)
#'
#' @return List with components:
#'   \describe{
#'     \item{stepdown}{Data frame from wy_stepdown}
#'     \item{sig_parcels}{Integer vector of significant parcel labels}
#'     \item{parcel_indices}{List of index vectors for each parcel}
#'   }
#'
#' @keywords internal
test_parcels <- function(z_vec, pi_vec, lab_vec, n_perm, kappa_grid,
                         alpha = 0.05, seed = NULL,
                         parallel = FALSE) {

  # Precompute parcel constants
  pi2_vec <- pi_vec * pi_vec
  den1 <- as.numeric(rowsum(pi_vec, group = lab_vec, reorder = FALSE))
  den2 <- as.numeric(rowsum(pi2_vec, group = lab_vec, reorder = FALSE))
  log_den1 <- log(den1)
  inv_sqrt_den2 <- 1 / sqrt(den2)

  parcels <- sort(unique(lab_vec))
  n_parcels <- length(parcels)

  # Observed parcel scores
  t_obs <- .score_parcels_vec(z_vec, pi_vec, lab_vec,
                              log_den1, inv_sqrt_den2, kappa_grid)

  # Null parcel scores
  if (!is.null(seed)) set.seed(seed)
  perm_fun <- make_perm_fun(z_vec, n_perm, seed = NULL)

  if (parallel) {
    t_null_list <- future.apply::future_lapply(seq_len(n_perm), function(b) {
      z_perm <- perm_fun(b)
      .score_parcels_vec(z_perm, pi_vec, lab_vec,
                         log_den1, inv_sqrt_den2, kappa_grid)
    }, future.seed = TRUE)
    t_null <- do.call(rbind, t_null_list)
  } else {
    t_null <- matrix(NA_real_, nrow = n_perm, ncol = n_parcels)
    for (b in seq_len(n_perm)) {
      z_perm <- perm_fun(b)
      t_null[b, ] <- .score_parcels_vec(z_perm, pi_vec, lab_vec,
                                        log_den1, inv_sqrt_den2, kappa_grid)
    }
  }

  # Step-down
  sd <- wy_stepdown(t_obs, t_null, alpha = alpha)

  # Get significant parcels
  sig_idx <- which(sd$rejected)
  sig_parcels <- parcels[sig_idx]

  # Build parcel index lists
  parcel_indices <- lapply(parcels, function(p) which(lab_vec == p))
  names(parcel_indices) <- as.character(parcels)

  list(
    stepdown = sd,
    sig_parcels = sig_parcels,
    parcel_indices = parcel_indices,
    den1 = den1
  )
}
