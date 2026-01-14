#' Threshold-Free Cluster Enhancement (TFCE)
#'
#' Implementation of FSL-style TFCE for threshold-free cluster
#' enhancement with permutation-based FWER control.
#'
#' @name tfce-methods
NULL

#' TFCE Transform
#'
#' Computes the Threshold-Free Cluster Enhancement transform:
#' \deqn{TFCE(v) = \int_0^{h_v} e(h)^E \cdot h^H \, dh}
#'
#' where e(h) is the extent of the cluster containing voxel v
#' at threshold h.
#'
#' @param stat_vol NeuroVol containing statistic map
#' @param mask Optional LogicalNeuroVol for analysis domain
#' @param H Height exponent (default 2.0)
#' @param E Extent exponent (default 0.5)
#' @param dh Threshold increment (default 0.1)
#' @param tail "pos" for activations, "neg" for deactivations
#'
#' @return NeuroVol with TFCE-transformed values
#'
#' @details
#' TFCE provides threshold-free inference by integrating cluster
#' extent over a range of thresholds. Default parameters (H=2, E=0.5)
#' are those recommended by Smith & Nichols (2009).
#'
#' @references
#' Smith, S. M., & Nichols, T. E. (2009). Threshold-free cluster
#' enhancement: Addressing problems of smoothing, threshold dependence
#' and localisation in cluster inference. NeuroImage, 44(1), 83-98.
#'
#' @examples
#' \dontrun{
#' tfce_map <- tfce_transform(z_map, H = 2, E = 0.5)
#' }
#'
#' @export
tfce_transform <- function(stat_vol, mask = NULL,
                           H = 2.0, E = 0.5, dh = 0.1,
                           tail = c("pos", "neg")) {
  tail <- match.arg(tail)

  if (is.null(mask)) {
    mask <- is.finite(as.numeric(stat_vol))
  }

  z <- as.numeric(stat_vol)
  if (tail == "neg") z <- -z

  dims <- dim(stat_vol)
  n <- prod(dims)

  # Only positive support
  z[!mask] <- 0
  z[z < 0] <- 0

  h_max <- max(z, na.rm = TRUE)
  if (h_max <= 0) {
    tfce_vals <- rep(0, n)
    return(.array_to_vol(tfce_vals, dims, stat_vol))
  }

  hs <- seq(0, h_max, by = dh)
  tfce_vals <- rep(0, n)

  for (h in hs) {
    if (h == 0) next

    # Find clusters at this threshold
    supra <- z >= h
    clusters <- .find_clusters(supra, dims)

    if (length(clusters$sizes) == 0) next

    # For each cluster, add contribution
    for (cid in seq_along(clusters$sizes)) {
      extent <- clusters$sizes[cid]
      voxels <- which(clusters$labels == cid)
      contribution <- (extent^E) * (h^H) * dh
      tfce_vals[voxels] <- tfce_vals[voxels] + contribution
    }
  }

  .array_to_vol(tfce_vals, dims, stat_vol)
}

#' TFCE with Permutation FWER
#'
#' Computes TFCE with permutation-based FWER control using the
#' maxT procedure.
#'
#' @param stat_vol NeuroVol containing statistic map
#' @param mask Optional LogicalNeuroVol for analysis domain
#' @param n_perm Number of permutations (default 5000)
#' @param H Height exponent (default 2.0)
#' @param E Extent exponent (default 0.5)
#' @param dh Threshold increment (default 0.1)
#' @param alpha Significance level (default 0.05)
#' @param seed Random seed for reproducibility
#' @param tail "pos", "neg", or "two"
#'
#' @return List with components:
#'   \describe{
#'     \item{tfce}{TFCE-transformed observed map}
#'     \item{threshold}{FWER-corrected threshold}
#'     \item{sig_mask}{Logical mask of significant voxels}
#'     \item{p_corr}{Voxelwise FWER-corrected p-values}
#'   }
#'
#' @export
tfce_fwer <- function(stat_vol, mask = NULL, n_perm = 5000,
                      H = 2.0, E = 0.5, dh = 0.1,
                      alpha = 0.05, seed = NULL,
                      tail = c("pos", "neg", "two")) {
  tail <- match.arg(tail)

  if (!is.null(seed)) set.seed(seed)

  if (is.null(mask)) {
    mask <- is.finite(as.numeric(stat_vol))
  }

  z_obs <- as.numeric(stat_vol)
  dims <- dim(stat_vol)
  n <- prod(dims)
  mask_idx <- which(mask)

  if (tail == "two") {
    return(.tfce_fwer_twosided(stat_vol, mask, n_perm, H, E, dh, alpha, seed))
  }

  # Observed TFCE
  tfce_obs <- tfce_transform(stat_vol, mask, H, E, dh, tail = tail)
  tfce_obs_vec <- as.numeric(tfce_obs)

  # Null distribution via sign-flipping
  max_null <- numeric(n_perm)

  for (b in seq_len(n_perm)) {
    signs <- sample(c(-1, 1), n, replace = TRUE)
    z_perm <- z_obs * signs

    # Create permuted volume (simple array for now)
    tfce_perm <- .tfce_from_vec(z_perm, dims, mask, H, E, dh, tail)
    max_null[b] <- max(tfce_perm[mask_idx], na.rm = TRUE)
  }

  # Threshold at (1-alpha) quantile
  threshold <- quantile(max_null, probs = 1 - alpha, na.rm = TRUE)

  # Significance mask
  sig_mask <- mask & FALSE
  sig_mask[tfce_obs_vec >= threshold & mask] <- TRUE

  # Voxelwise corrected p-values
  p_corr <- rep(1, n)
  for (v in mask_idx) {
    p_corr[v] <- (1 + sum(max_null >= tfce_obs_vec[v])) / (n_perm + 1)
  }

  list(
    tfce = tfce_obs,
    threshold = threshold,
    sig_mask = sig_mask,
    p_corr = .array_to_vol(p_corr, dims, stat_vol),
    max_null = max_null,
    method = "tfce_fwer",
    params = list(alpha = alpha, H = H, E = E, dh = dh, n_perm = n_perm)
  )
}

#' Two-sided TFCE with combined max null
#' @keywords internal
.tfce_fwer_twosided <- function(stat_vol, mask, n_perm,
                                H, E, dh, alpha, seed) {
  if (!is.null(seed)) set.seed(seed)

  z_obs <- as.numeric(stat_vol)
  dims <- dim(stat_vol)
  n <- prod(dims)
  mask_idx <- which(mask)

  # Observed TFCE for both tails
  tfce_pos <- tfce_transform(stat_vol, mask, H, E, dh, tail = "pos")
  tfce_neg <- tfce_transform(stat_vol, mask, H, E, dh, tail = "neg")

  tfce_pos_vec <- as.numeric(tfce_pos)
  tfce_neg_vec <- as.numeric(tfce_neg)
  tfce_2 <- pmax(tfce_pos_vec, tfce_neg_vec)

  # Null distribution of max two-sided TFCE
  max_null <- numeric(n_perm)

  for (b in seq_len(n_perm)) {
    signs <- sample(c(-1, 1), n, replace = TRUE)
    z_perm <- z_obs * signs

    tfce_p <- .tfce_from_vec(z_perm, dims, mask, H, E, dh, "pos")
    tfce_n <- .tfce_from_vec(-z_perm, dims, mask, H, E, dh, "pos")
    max_null[b] <- max(pmax(tfce_p, tfce_n)[mask_idx], na.rm = TRUE)
  }

  threshold <- quantile(max_null, probs = 1 - alpha, na.rm = TRUE)

  sig_pos <- tfce_pos_vec >= threshold & tfce_pos_vec >= tfce_neg_vec & mask
  sig_neg <- tfce_neg_vec >= threshold & tfce_neg_vec > tfce_pos_vec & mask
  sig_any <- tfce_2 >= threshold & mask

  list(
    tfce_pos = tfce_pos,
    tfce_neg = tfce_neg,
    tfce_2 = .array_to_vol(tfce_2, dims, stat_vol),
    threshold = threshold,
    sig_mask = .array_to_vol(sig_any, dims, stat_vol),
    sig_pos = .array_to_vol(sig_pos, dims, stat_vol),
    sig_neg = .array_to_vol(sig_neg, dims, stat_vol),
    max_null = max_null,
    method = "tfce_fwer_twosided",
    params = list(alpha = alpha, H = H, E = E, dh = dh, n_perm = n_perm)
  )
}

#' TFCE from vector (internal helper)
#' @keywords internal
.tfce_from_vec <- function(z_vec, dims, mask, H, E, dh, tail) {
  z <- z_vec
  if (tail == "neg") z <- -z

  n <- prod(dims)
  z[!mask] <- 0
  z[z < 0] <- 0

  h_max <- max(z, na.rm = TRUE)
  if (h_max <= 0) return(rep(0, n))

  hs <- seq(0, h_max, by = dh)
  tfce_vals <- rep(0, n)

  for (h in hs) {
    if (h == 0) next
    supra <- z >= h
    clusters <- .find_clusters(supra, dims)
    if (length(clusters$sizes) == 0) next

    for (cid in seq_along(clusters$sizes)) {
      extent <- clusters$sizes[cid]
      voxels <- which(clusters$labels == cid)
      tfce_vals[voxels] <- tfce_vals[voxels] + (extent^E) * (h^H) * dh
    }
  }

  tfce_vals
}

#' Convert array to vol-like object (placeholder)
#' @keywords internal
.array_to_vol <- function(vals, dims, template) {
  if (length(vals) != prod(dims)) {
    stop("vals length does not match dims")
  }

  out <- template
  ok <- TRUE
  tryCatch({
    out[] <- vals
  }, error = function(e) {
    ok <<- FALSE
  })

  if (ok) {
    return(out)
  }

  array(vals, dim = dims)
}
