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
#' Computes the FWER-corrected threshold for peak (voxelwise) inference using
#' Random Field Theory.
#'
#' @param stat_vol 3D Z-equivalent statistic map (array or `NeuroVol`).
#'   Pass already-canonicalized input (use [canonicalize()] first for
#'   t-statistics or -log10(p) maps).
#' @param mask Optional logical 3D mask. Defaults to finite voxels in
#'   `stat_vol`.
#' @param alpha Significance level (default 0.05).
#' @param tail One of `"pos"`, `"neg"`, `"two"`. For `"two"`, each direction
#'   uses `alpha / 2` and the union is reported, targeting FWER `alpha`
#'   across both tails.
#' @param smoothness A [`smoothness_fwhm()`] object specifying field
#'   smoothness in millimetres. Numeric input is also accepted (length 1 or
#'   3, mm) and coerced internally.
#' @param df Degrees of freedom. Use `Inf` for Gaussian field (Z-scores).
#'
#' @return A [`neurothresh_result`] with `method = "rft_peak_fwer"`. Slots:
#'   `threshold` is the positive-side critical value; `reject` is the 3D
#'   logical rejection volume; `pvals` is `NA` (RFT does not produce per-voxel
#'   corrected p-values without further computation).
#'
#' @details
#' Uses the expected Euler characteristic approximation:
#' \deqn{P(\max Z > u) \approx \sum_{d=0}^{D} \mathrm{Resels}_d \cdot EC_d(u)}
#'
#' For high thresholds in smooth Gaussian fields, this provides
#' accurate FWER control.
#'
#' @seealso
#' [rft_cluster_fwer()] for cluster-level inference;
#' [tfce_fwer()] and [cluster_fdr()] for alternative methods.
#'
#' @export
rft_peak_fwer <- function(stat_vol, mask = NULL, alpha = 0.05,
                          tail = c("pos", "neg", "two"),
                          smoothness, df = Inf) {
  tail <- match.arg(tail)
  call_capture <- match.call()

  if (is.null(mask)) {
    mask <- is.finite(as.numeric(stat_vol))
  }

  dims <- dim(stat_vol)
  if (is.null(dims) || length(dims) != 3L) {
    stop("stat_vol must be a 3D volume-like object")
  }

  fwhm_mm <- .resolve_smoothness_fwhm_mm(smoothness)
  vox_size <- .vox_size_mm(stat_vol)
  n_voxels <- sum(mask)
  if (n_voxels == 0L) stop("No valid voxels in mask")

  alpha_tail <- if (tail == "two") alpha / 2 else alpha
  threshold <- .rft_peak_threshold(dims, fwhm_mm, vox_size, n_voxels, alpha_tail, df)

  z <- as.numeric(stat_vol)
  reject_vec <- switch(tail,
    pos = (z >= threshold) & mask,
    neg = (-z >= threshold) & mask,
    two = ((z >= threshold) | (-z >= threshold)) & mask
  )
  reject <- array(reject_vec, dim = dims)

  pvals <- array(NA_real_, dim = dims)

  new_neurothresh_result(
    method = "rft_peak_fwer",
    stat_vol = stat_vol,
    mask = mask,
    threshold = threshold,
    reject = reject,
    pvals = pvals,
    params = list(alpha = alpha, tail = tail,
                  smoothness = smoothness, fwhm_mm = fwhm_mm,
                  df = df, n_voxels = n_voxels),
    call = call_capture
  )
}

#' RFT Cluster-Level FWER
#'
#' Computes cluster-level FWER-corrected p-values using Random Field Theory
#' cluster extent distribution.
#'
#' @inheritParams rft_peak_fwer
#' @param cluster_thresh Cluster-forming threshold on the Z scale.
#'
#' @return A [`neurothresh_result`] with `method = "rft_cluster_fwer"`. The
#'   `clusters` slot is the per-cluster table (id, size, p_unc, p_fwer,
#'   significant). `threshold` is `cluster_thresh`.
#'
#' @details First thresholds at `cluster_thresh` to form clusters, then
#' tests each cluster's extent against the RFT null distribution. For
#' `tail = "two"`, the analysis runs separately on positive and negative
#' tails at `alpha / 2` and the union of significant clusters is reported.
#' The reported `p_fwer` values are doubled (capped at 1) to describe the
#' combined family; `p_unc` remains the one-tail cluster-extent probability.
#'
#' Extent probabilities use a three-dimensional, high-threshold Gaussian
#' approximation in dimensionless resels. This is an asymptotic baseline;
#' changing physical units does not change its probabilities.
#'
#' @seealso [rft_peak_fwer()], [cluster_fdr()].
#' @export
rft_cluster_fwer <- function(stat_vol, mask = NULL, alpha = 0.05,
                             tail = c("pos", "neg", "two"),
                             cluster_thresh = 3.0,
                             smoothness, df = Inf) {
  tail <- match.arg(tail)
  call_capture <- match.call()

  if (is.null(mask)) {
    mask <- is.finite(as.numeric(stat_vol))
  }
  dims <- dim(stat_vol)
  if (is.null(dims) || length(dims) != 3L) {
    stop("stat_vol must be a 3D volume-like object")
  }

  fwhm_mm <- .resolve_smoothness_fwhm_mm(smoothness)
  vox_size <- .vox_size_mm(stat_vol)

  if (tail %in% c("pos", "neg")) {
    res <- .rft_cluster_one_tail(stat_vol, mask, dims, fwhm_mm, vox_size,
                                  cluster_thresh, alpha, tail)
    return(new_neurothresh_result(
      method = "rft_cluster_fwer",
      stat_vol = stat_vol,
      mask = mask,
      threshold = cluster_thresh,
      reject = res$reject,
      pvals = array(NA_real_, dim = dims),
      clusters = res$table,
      params = list(alpha = alpha, tail = tail,
                    cluster_thresh = cluster_thresh,
                    smoothness = smoothness, fwhm_mm = fwhm_mm,
                    df = df),
      call = call_capture
    ))
  }

  pos_res <- .rft_cluster_one_tail(stat_vol, mask, dims, fwhm_mm, vox_size,
                                    cluster_thresh, alpha / 2, "pos")
  neg_res <- .rft_cluster_one_tail(stat_vol, mask, dims, fwhm_mm, vox_size,
                                    cluster_thresh, alpha / 2, "neg")
  reject <- array(as.logical(pos_res$reject) | as.logical(neg_res$reject),
                  dim = dims)

  pos_tab <- pos_res$table
  neg_tab <- neg_res$table
  if (nrow(pos_tab) > 0L) pos_tab$tail <- "pos"
  if (nrow(neg_tab) > 0L) {
    neg_tab$tail <- "neg"
    neg_tab$cluster_id <- neg_tab$cluster_id + nrow(pos_tab)
  }
  table <- rbind(pos_tab, neg_tab)
  table$p_fwer <- pmin(1, 2 * table$p_fwer)

  new_neurothresh_result(
    method = "rft_cluster_fwer",
    stat_vol = stat_vol,
    mask = mask,
    threshold = cluster_thresh,
    reject = reject,
    pvals = array(NA_real_, dim = dims),
    clusters = table,
    params = list(alpha = alpha, tail = "two",
                  cluster_thresh = cluster_thresh,
                  smoothness = smoothness, fwhm_mm = fwhm_mm,
                  df = df),
    call = call_capture
  )
}

#' @keywords internal
.rft_peak_threshold <- function(dims, fwhm_mm, vox_size, n_voxels, alpha, df) {
  if (is.infinite(df)) {
    p_rft <- function(u) {
      .rft_expected_ec_gaussian(u, dims = dims, fwhm_mm = fwhm_mm,
                                vox_size = vox_size)
    }
    u_bonf <- qnorm(1 - alpha / n_voxels)
    u_rft <- tryCatch(
      uniroot(function(u) p_rft(u) - alpha, c(2, 10))$root,
      error = function(e) u_bonf
    )
    max(u_bonf, u_rft)
  } else {
    qt(1 - alpha / n_voxels, df)
  }
}

#' @keywords internal
.rft_cluster_one_tail <- function(stat_vol, mask, dims, fwhm_mm, vox_size,
                                   cluster_thresh, alpha, tail) {
  z <- as.numeric(stat_vol)
  if (tail == "neg") z <- -z

  supra <- (z >= cluster_thresh) & mask
  empty_table <- data.frame(
    cluster_id = integer(0), size_voxels = integer(0),
    size_resels = numeric(0), p_unc = numeric(0),
    p_fwer = numeric(0), significant = logical(0)
  )

  if (sum(supra) == 0L) {
    return(list(reject = array(FALSE, dim = dims), table = empty_table))
  }

  clusters <- .find_clusters(supra, dims)
  if (length(clusters$sizes) == 0L) {
    return(list(reject = array(FALSE, dim = dims), table = empty_table))
  }

  s_resel <- clusters$sizes * prod(vox_size) / prod(fwhm_mm)
  u0 <- cluster_thresh
  p_unc <- .rft_cluster_extent_tail(s_resel, u0)
  EK <- .rft_expected_ec_gaussian(u0, dims = dims, fwhm_mm = fwhm_mm,
                                   vox_size = vox_size)
  p_fwer <- pmin(1, EK * p_unc)

  table <- data.frame(
    cluster_id = seq_along(clusters$sizes),
    size_voxels = clusters$sizes,
    size_resels = s_resel,
    p_unc = p_unc,
    p_fwer = p_fwer,
    significant = p_fwer <= alpha
  )

  sig_ids <- which(table$significant)
  reject_vec <- mask & FALSE
  for (id in sig_ids) {
    reject_vec[clusters$labels == id] <- TRUE
  }

  list(reject = array(reject_vec, dim = dims), table = table)
}

#' @noRd
.rft_cluster_extent_tail <- function(size_resels, u) {
  # Nosko's 3D high-threshold extent approximation (Friston et al., 1994,
  # Assessing the significance of focal activations using their spatial extent).
  # E[K] ~ (2*pi)^(3/2) / ((4*log(2))^(3/2) * u^3), in resels.
  # beta = (Gamma(5/2) / E[K])^(2/3); P(K >= k) ~ exp(-beta*k^(2/3)).
  # The extent has already been converted to resels: no physical FWHM factor
  # belongs in beta. This retains the existing high-u approximation.
  beta <- gamma(5 / 2)^(2 / 3) * (4 * log(2)) * u^2 / (2 * pi)
  exp(-beta * size_resels^(2 / 3))
}

#' @keywords internal
.resolve_smoothness_fwhm_mm <- function(s) {
  if (inherits(s, "neurothresh_smoothness_fwhm")) return(s$fwhm_mm)
  if (inherits(s, "neurothresh_smoothness")) {
    stop("rft_*_fwer() requires a smoothness_fwhm() (got ",
         class(s)[1L], ")")
  }
  if (is.numeric(s)) {
    if (length(s) == 1L) s <- rep(s, 3L)
    if (length(s) != 3L) stop("numeric smoothness must have length 1 or 3 (mm)")
    return(s)
  }
  stop("smoothness must be a smoothness_fwhm() or a length-1/3 numeric mm")
}

#' @noRd
.rft_ec_densities_gaussian <- function(u) {
  u <- as.numeric(u)
  e <- exp(-u^2 / 2)
  c1 <- sqrt(4 * log(2)) / (2 * pi)
  c2 <- (4 * log(2)) / ((2 * pi)^(3 / 2))
  c3 <- (4 * log(2))^(3 / 2) / ((2 * pi)^2)
  c(
    rho0 = 1 - pnorm(u),
    rho1 = c1 * e,
    rho2 = c2 * u * e,
    rho3 = c3 * (u^2 - 1) * e
  )
}

#' @noRd
.rft_lkc_box <- function(dims, fwhm_mm, vox_size) {
  if (length(dims) != 3) stop("dims must be length 3")
  if (length(fwhm_mm) == 1) fwhm_mm <- rep(fwhm_mm, 3)
  if (length(fwhm_mm) != 3) stop("fwhm_mm must be length 1 or 3")
  if (length(vox_size) == 1) vox_size <- rep(vox_size, 3)
  if (length(vox_size) != 3) stop("vox_size must be length 3")

  Lx <- (dims[1] * vox_size[1]) / fwhm_mm[1]
  Ly <- (dims[2] * vox_size[2]) / fwhm_mm[2]
  Lz <- (dims[3] * vox_size[3]) / fwhm_mm[3]

  c(
    lkc0 = 1,
    lkc1 = Lx + Ly + Lz,
    lkc2 = Lx * Ly + Lx * Lz + Ly * Lz,
    lkc3 = Lx * Ly * Lz
  )
}

#' @noRd
.rft_expected_ec_gaussian <- function(u, dims, fwhm_mm, vox_size) {
  lkc <- .rft_lkc_box(dims = dims, fwhm_mm = fwhm_mm, vox_size = vox_size)
  rho <- .rft_ec_densities_gaussian(u)
  sum(lkc * rho)
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

#' Find connected components (simple 6-connectivity)
#' @keywords internal
.find_clusters <- function(binary_mask, dims) {
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
  idx0 <- idx - 1
  i <- idx0 %% dims[1]
  j <- (idx0 %/% dims[1]) %% dims[2]
  k <- idx0 %/% (dims[1] * dims[2])

  neighbors <- integer(0)

  if (i > 0) neighbors <- c(neighbors, idx - 1)
  if (i < dims[1] - 1) neighbors <- c(neighbors, idx + 1)

  if (j > 0) neighbors <- c(neighbors, idx - dims[1])
  if (j < dims[2] - 1) neighbors <- c(neighbors, idx + dims[1])

  if (k > 0) neighbors <- c(neighbors, idx - dims[1] * dims[2])
  if (k < dims[3] - 1) neighbors <- c(neighbors, idx + dims[1] * dims[2])

  neighbors
}
