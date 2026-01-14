#' Cluster-Level FDR Control
#'
#' Implements topological FDR control at the cluster level using
#' Benjamini-Hochberg procedure applied to cluster p-values.
#'
#' @name cluster-fdr
NULL

#' Cluster-FDR with RFT or Permutation P-values
#'
#' Controls the False Discovery Rate at the cluster level by applying
#' the Benjamini-Hochberg procedure to cluster-extent p-values.
#'
#' @param stat_vol NeuroVol containing statistic map
#' @param mask Optional LogicalNeuroVol for analysis domain
#' @param cluster_thresh Cluster-forming threshold (Z-score)
#' @param q FDR level (default 0.05)
#' @param p_method Method for computing cluster p-values:
#'   "rft" (Random Field Theory) or "perm" (permutation)
#' @param fwhm_mm Smoothness in mm (required for p_method = "rft")
#' @param n_perm Number of permutations (for p_method = "perm")
#' @param tail "pos", "neg", or "two"
#' @param two_sided_policy For two-sided: "BH_all" (pool both tails) or
#'   "split_q" (q/2 per tail)
#' @param seed Random seed for permutation
#'
#' @return List with components:
#'   \describe{
#'     \item{sig_mask}{Logical mask of FDR-significant clusters}
#'     \item{table}{Data frame with cluster statistics}
#'     \item{n_clusters}{Number of significant clusters}
#'   }
#'
#' @details
#' This implements "topological FDR" as advocated by Chumbley & Friston,
#' where FDR is controlled over the finite set of clusters rather than
#' over all voxels.
#'
#' @references
#' Chumbley, J. R., & Friston, K. J. (2009). False discovery rate
#' revisited: FDR and topological inference using Gaussian random fields.
#' NeuroImage, 44(1), 62-70.
#'
#' @examples
#' \dontrun{
#' # RFT-based cluster FDR
#' result <- cluster_fdr(z_map, cluster_thresh = 3.0, q = 0.05,
#'                       p_method = "rft", fwhm_mm = 8)
#'
#' # Permutation-based cluster FDR
#' result <- cluster_fdr(z_map, cluster_thresh = 3.0, q = 0.05,
#'                       p_method = "perm", n_perm = 5000)
#' }
#'
#' @export
cluster_fdr <- function(stat_vol, mask = NULL,
                        cluster_thresh = 3.0, q = 0.05,
                        p_method = c("rft", "perm"),
                        fwhm_mm = NULL, n_perm = 5000,
                        tail = c("pos", "neg", "two"),
                        two_sided_policy = c("BH_all", "split_q"),
                        seed = NULL) {

  p_method <- match.arg(p_method)
  tail <- match.arg(tail)
  two_sided_policy <- match.arg(two_sided_policy)

  if (p_method == "rft" && is.null(fwhm_mm)) {
    stop("fwhm_mm required for p_method='rft'")
  }

  if (is.null(mask)) {
    mask <- is.finite(as.numeric(stat_vol))
  }

  if (tail == "two") {
    return(.cluster_fdr_twosided(stat_vol, mask, cluster_thresh, q,
                                 p_method, fwhm_mm, n_perm,
                                 two_sided_policy, seed))
  }

  z <- as.numeric(stat_vol)
  if (tail == "neg") z <- -z
  dims <- dim(stat_vol)

  # Find clusters
  supra <- (z >= cluster_thresh) & mask
  clusters <- .find_clusters(supra, dims)

  if (length(clusters$sizes) == 0) {
    return(list(
      sig_mask = mask & FALSE,
      table = data.frame(),
      n_clusters = 0,
      method = "cluster_fdr"
    ))
  }

  # Compute cluster p-values
  if (p_method == "rft") {
    vox_size <- .vox_size_mm(stat_vol)
    p_unc <- .cluster_pvals_rft(clusters$sizes, cluster_thresh, fwhm_mm, vox_size)
  } else {
    p_unc <- .cluster_pvals_perm(z, mask, cluster_thresh, clusters, n_perm, seed)
  }

  # Benjamini-Hochberg
  rejected <- .bh_reject(p_unc, q)

  # Build table
  cluster_table <- data.frame(
    cluster_id = seq_along(clusters$sizes),
    size_voxels = clusters$sizes,
    p_unc = p_unc,
    significant = rejected,
    tail = tail
  )

  # Create significance mask
  sig_mask <- mask & FALSE
  for (id in which(rejected)) {
    sig_mask[clusters$labels == id] <- TRUE
  }

  list(
    sig_mask = sig_mask,
    table = cluster_table,
    n_clusters = sum(rejected),
    cluster_labels = clusters$labels,
    method = "cluster_fdr",
    params = list(q = q, cluster_thresh = cluster_thresh,
                  p_method = p_method, tail = tail)
  )
}

#' Two-sided cluster FDR
#' @keywords internal
.cluster_fdr_twosided <- function(stat_vol, mask, cluster_thresh, q,
                                  p_method, fwhm_mm, n_perm,
                                  two_sided_policy, seed) {

  # Run for both tails
  pos <- cluster_fdr(stat_vol, mask, cluster_thresh, q,
                     p_method, fwhm_mm, n_perm, tail = "pos", seed = seed)
  neg <- cluster_fdr(stat_vol, mask, cluster_thresh, q,
                     p_method, fwhm_mm, n_perm, tail = "neg", seed = seed)

  if (two_sided_policy == "split_q") {
    # Already used q for each tail, now combine
    sig_mask <- pos$sig_mask | neg$sig_mask
  } else {
    # Pool p-values and apply BH once
    p_all <- c(pos$table$p_unc, neg$table$p_unc)
    if (length(p_all) > 0) {
      rej_all <- .bh_reject(p_all, q)
      n_pos <- nrow(pos$table)

      rej_pos <- rej_all[seq_len(n_pos)]
      rej_neg <- rej_all[(n_pos + 1):length(rej_all)]

      pos$table$significant <- rej_pos
      neg$table$significant <- rej_neg

      # Rebuild masks
      sig_pos <- mask & FALSE
      for (id in which(rej_pos)) {
        sig_pos[pos$cluster_labels == id] <- TRUE
      }

      sig_neg <- mask & FALSE
      for (id in which(rej_neg)) {
        sig_neg[neg$cluster_labels == id] <- TRUE
      }

      sig_mask <- sig_pos | sig_neg
    } else {
      sig_mask <- mask & FALSE
    }
  }

  # Combine tables
  table <- rbind(
    cbind(pos$table, tail = "pos"),
    cbind(neg$table, tail = "neg")
  )

  list(
    sig_mask = sig_mask,
    sig_pos = pos$sig_mask,
    sig_neg = neg$sig_mask,
    table = table,
    n_clusters = sum(table$significant),
    method = "cluster_fdr_twosided",
    params = list(q = q, cluster_thresh = cluster_thresh,
                  p_method = p_method, two_sided_policy = two_sided_policy)
  )
}

#' RFT cluster p-values
#' @keywords internal
.cluster_pvals_rft <- function(sizes, u0, fwhm_mm, vox_size) {
  if (length(fwhm_mm) == 1) fwhm_mm <- rep(fwhm_mm, 3)

  # Convert to resels
  s_resel <- sizes * prod(vox_size) / prod(fwhm_mm)

  # Constant c from Worsley
  c_const <- prod(fwhm_mm) * (2 * pi / u0)^(3/2) *
    (4 * log(2))^(-3/2) / gamma(5/2)

  # Uncorrected cluster p-value
  p_unc <- exp(-u0 * (s_resel / c_const)^(2/3))
  pmin(pmax(p_unc, 0), 1)
}

#' Permutation cluster p-values
#' @keywords internal
.cluster_pvals_perm <- function(z, mask, u0, obs_clusters, n_perm, seed) {
  if (!is.null(seed)) set.seed(seed)

  dims <- dim(z)
  n <- length(z)
  obs_sizes <- obs_clusters$sizes

  # Collect null cluster sizes
  null_sizes <- numeric(0)

  for (b in seq_len(n_perm)) {
    signs <- sample(c(-1, 1), n, replace = TRUE)
    z_perm <- z * signs

    supra <- (z_perm >= u0) & mask
    perm_clusters <- .find_clusters(supra, dims)

    if (length(perm_clusters$sizes) > 0) {
      null_sizes <- c(null_sizes, perm_clusters$sizes)
    }
  }

  # Compute p-values
  if (length(null_sizes) == 0) {
    p_unc <- rep(1 / (n_perm + 1), length(obs_sizes))
  } else {
    p_unc <- vapply(obs_sizes, function(s) {
      (1 + sum(null_sizes >= s)) / (1 + length(null_sizes))
    }, numeric(1))
  }

  p_unc
}

#' Benjamini-Hochberg rejection
#' @keywords internal
.bh_reject <- function(p, q) {
  m <- length(p)
  if (m == 0) return(logical(0))

  o <- order(p)
  ps <- p[o]
  thr <- (seq_len(m) / m) * q
  k <- max(which(ps <= thr), 0L)

  rej <- rep(FALSE, m)
  if (k > 0) rej[o[seq_len(k)]] <- TRUE
  rej
}
