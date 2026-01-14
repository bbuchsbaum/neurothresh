#' Random Field Theory (RFT) Thresholding Methods
#'
#' Classic GRF/RFT inference using expected Euler characteristic
#' for FWER control. These are provided as baseline comparisons
#' for LR-MFT.
#'
#' @name rft-methods
NULL

#' RFT Peak-Level FWER
#'
#' Computes the FWER-corrected threshold for peak (voxelwise)
#' inference using Random Field Theory.
#'
#' @param stat_vol NeuroVol containing statistic map
#' @param mask Optional LogicalNeuroVol for analysis domain
#' @param fwhm_mm Smoothness in mm. Scalar or length-3 vector for
#'   anisotropic smoothness.
#' @param alpha Significance level (default 0.05)
#' @param df Degrees of freedom. Use Inf for Gaussian field (Z-scores).
#' @param tail "pos" (activations), "neg" (deactivations), or "two"
#'
#' @return List with components:
#'   \describe{
#'     \item{threshold}{The FWER-corrected threshold}
#'     \item{sig_mask}{Logical mask of significant voxels}
#'     \item{method}{"rft_peak"}
#'   }
#'
#' @details
#' Uses the expected Euler characteristic approximation:
#' \deqn{P(\max Z > u) \approx \sum_{d=0}^{D} \text{Resels}_d \cdot EC_d(u)}
#'
#' For high thresholds in smooth Gaussian fields, this provides
#' accurate FWER control.
#'
#' @examples
#' \dontrun{
#' result <- rft_peak_fwer(z_map, fwhm_mm = 8, alpha = 0.05)
#' }
#'
#' @export
rft_peak_fwer <- function(stat_vol, mask = NULL, fwhm_mm, alpha = 0.05,
                          df = Inf, tail = c("pos", "neg", "two")) {
  tail <- match.arg(tail)

  if (is.null(mask)) {
    mask <- is.finite(as.numeric(stat_vol))
  }

  z <- as.numeric(stat_vol)
  if (tail == "neg") z <- -z
  if (tail == "two") z <- abs(z)

  # Estimate search volume in resels
  n_voxels <- sum(mask)

  if (length(fwhm_mm) == 1) {
    fwhm_mm <- rep(fwhm_mm, 3)
  }

  # Get voxel dimensions (assume 2mm isotropic if not available)
  vox_size <- .vox_size_mm(stat_vol)

  # Resels
  resel3 <- prod(vox_size) * n_voxels / prod(fwhm_mm)

  # Find threshold using RFT approximation
  # For Gaussian field: EC_3(u) = (2*pi)^(-2) * u * exp(-u^2/2)
  if (is.infinite(df)) {
    # Gaussian field
    p_rft <- function(u) {
      # 3D Euler characteristic density for Gaussian field
      ec3 <- (2 * pi)^(-2) * u * exp(-u^2 / 2)
      resel3 * ec3
    }

    # Bonferroni bound for comparison
    u_bonf <- qnorm(1 - alpha / n_voxels)

    # Find RFT threshold
    u_rft <- tryCatch({
      uniroot(function(u) p_rft(u) - alpha, c(2, 10))$root
    }, error = function(e) {
      u_bonf
    })

    threshold <- min(u_bonf, u_rft)
  } else {
    # t-field: use conservative Bonferroni
    threshold <- qt(1 - alpha / n_voxels, df)
  }

  sig_mask <- (z >= threshold) & mask

  list(
    threshold = threshold,
    sig_mask = sig_mask,
    n_significant = sum(sig_mask),
    method = "rft_peak",
    params = list(alpha = alpha, fwhm_mm = fwhm_mm, tail = tail)
  )
}

#' RFT Cluster-Level FWER
#'
#' Computes cluster-level FWER-corrected p-values using Random Field
#' Theory cluster extent distribution.
#'
#' @param stat_vol NeuroVol containing statistic map
#' @param mask Optional LogicalNeuroVol for analysis domain
#' @param fwhm_mm Smoothness in mm
#' @param cluster_thresh Cluster-forming threshold (Z-score)
#' @param alpha Significance level (default 0.05)
#' @param df Degrees of freedom. Use Inf for Gaussian field.
#' @param tail "pos", "neg", or "two"
#'
#' @return List with components:
#'   \describe{
#'     \item{sig_mask}{Logical mask of significant clusters}
#'     \item{table}{Data frame with cluster statistics}
#'     \item{n_clusters}{Number of significant clusters}
#'   }
#'
#' @details
#' First thresholds at cluster_thresh to form clusters, then tests
#' each cluster's extent against the RFT null distribution.
#'
#' @export
rft_cluster_fwer <- function(stat_vol, mask = NULL, fwhm_mm,
                             cluster_thresh = 3.0, alpha = 0.05,
                             df = Inf, tail = c("pos", "neg", "two")) {
  tail <- match.arg(tail)

  if (is.null(mask)) {
    mask <- is.finite(as.numeric(stat_vol))
  }

  z <- as.numeric(stat_vol)
  if (tail == "neg") z <- -z
  if (tail == "two") z <- abs(z)

  dims <- dim(stat_vol)

  # Find supra-threshold voxels
  supra <- (z >= cluster_thresh) & mask

  if (sum(supra) == 0) {
    return(list(
      sig_mask = mask & FALSE,
      table = data.frame(),
      n_clusters = 0,
      method = "rft_cluster"
    ))
  }

  # Find connected components
  clusters <- .find_clusters(supra, dims)

  if (length(clusters$sizes) == 0) {
    return(list(
      sig_mask = mask & FALSE,
      table = data.frame(),
      n_clusters = 0,
      method = "rft_cluster"
    ))
  }

  # Compute cluster p-values under RFT
  if (length(fwhm_mm) == 1) fwhm_mm <- rep(fwhm_mm, 3)
  vox_size <- .vox_size_mm(stat_vol)

  # Convert cluster sizes to resels
  s_resel <- clusters$sizes * prod(vox_size) / prod(fwhm_mm)

  # RFT cluster extent p-value (Worsley approximation)
  # P(S > s) ~ exp(-z0 * (s/c)^(2/D))
  # where c depends on field smoothness
  u0 <- cluster_thresh
  c_const <- prod(fwhm_mm) * (2 * pi / u0)^(3/2) *
    (4 * log(2))^(-3/2) / gamma(5/2)

  p_unc <- exp(-u0 * (s_resel / c_const)^(2/3))

  # Expected number of clusters for FWER correction
  n_voxels <- sum(mask)
  resel3 <- prod(vox_size) * n_voxels / prod(fwhm_mm)
  EK <- resel3 * (2 * pi)^(-2) * u0 * exp(-u0^2 / 2)

  p_fwer <- pmin(1, EK * p_unc)

  # Build output table
  cluster_table <- data.frame(
    cluster_id = seq_along(clusters$sizes),
    size_voxels = clusters$sizes,
    size_resels = s_resel,
    p_unc = p_unc,
    p_fwer = p_fwer,
    significant = p_fwer <= alpha
  )

  # Create significance mask
  sig_ids <- which(cluster_table$significant)
  sig_mask <- mask & FALSE
  for (id in sig_ids) {
    sig_mask[clusters$labels == id] <- TRUE
  }

  list(
    sig_mask = sig_mask,
    table = cluster_table,
    n_clusters = sum(cluster_table$significant),
    cluster_labels = clusters$labels,
    method = "rft_cluster",
    params = list(alpha = alpha, cluster_thresh = cluster_thresh,
                  fwhm_mm = fwhm_mm, tail = tail)
  )
}

#' @keywords internal
.vox_size_mm <- function(stat_vol) {
  vox_size <- NULL
  sp <- try(neuroim2::space(stat_vol), silent = TRUE)
  if (!inherits(sp, "try-error")) {
    vox_size <- tryCatch(sp@spacing, error = function(e) NULL)
    if (is.null(vox_size) && !is.null(sp$spacing)) {
      vox_size <- sp$spacing
    }
  }

  if (is.null(vox_size) || length(vox_size) < 3) {
    return(c(2, 2, 2))
  }

  as.numeric(abs(vox_size[1:3]))
}

#' Two-sided RFT wrapper
#'
#' @param stat_vol NeuroVol containing statistic map
#' @param method "peak" or "cluster"
#' @param ... Additional arguments passed to rft_peak_fwer or rft_cluster_fwer
#'
#' @export
rft_fwer_twosided <- function(stat_vol, method = c("peak", "cluster"), ...) {
  method <- match.arg(method)

  if (method == "peak") {
    pos <- rft_peak_fwer(stat_vol, tail = "pos", ...)
    neg <- rft_peak_fwer(stat_vol, tail = "neg", ...)
  } else {
    pos <- rft_cluster_fwer(stat_vol, tail = "pos", ...)
    neg <- rft_cluster_fwer(stat_vol, tail = "neg", ...)
  }

  sig_mask <- pos$sig_mask | neg$sig_mask

  list(
    sig_mask = sig_mask,
    pos = pos,
    neg = neg,
    method = paste0("rft_", method, "_twosided")
  )
}

#' Find connected components (simple 6-connectivity)
#' @keywords internal
.find_clusters <- function(binary_mask, dims) {
  # Simple flood-fill based clustering
  # For production, use neuroim2::conn_comp instead

  n <- prod(dims)
  labels <- integer(n)
  current_label <- 0

  linear_idx <- which(binary_mask)
  if (length(linear_idx) == 0) {
    return(list(labels = labels, sizes = integer(0)))
  }

  remaining <- linear_idx

  while (length(remaining) > 0) {
    current_label <- current_label + 1
    seed <- remaining[1]
    queue <- seed
    labels[seed] <- current_label

    while (length(queue) > 0) {
      current <- queue[1]
      queue <- queue[-1]

      # Get 6-connected neighbors
      neighbors <- .get_neighbors_6(current, dims)
      neighbors <- neighbors[binary_mask[neighbors] & labels[neighbors] == 0]

      labels[neighbors] <- current_label
      queue <- c(queue, neighbors)
    }

    remaining <- remaining[labels[remaining] == 0]
  }

  sizes <- tabulate(labels[labels > 0])

  list(labels = labels, sizes = sizes)
}

#' Get 6-connected neighbors
#' @keywords internal
.get_neighbors_6 <- function(idx, dims) {
  # Convert to i,j,k
  idx0 <- idx - 1
  i <- idx0 %% dims[1]
  j <- (idx0 %/% dims[1]) %% dims[2]
  k <- idx0 %/% (dims[1] * dims[2])

  neighbors <- integer(0)

  # x neighbors
  if (i > 0) neighbors <- c(neighbors, idx - 1)
  if (i < dims[1] - 1) neighbors <- c(neighbors, idx + 1)

  # y neighbors
  if (j > 0) neighbors <- c(neighbors, idx - dims[1])
  if (j < dims[2] - 1) neighbors <- c(neighbors, idx + dims[1])

  # z neighbors
  if (k > 0) neighbors <- c(neighbors, idx - dims[1] * dims[2])
  if (k < dims[3] - 1) neighbors <- c(neighbors, idx + dims[1] * dims[2])

  neighbors
}
