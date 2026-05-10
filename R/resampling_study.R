#' Simulation Study: Resampling Choice and Null Rejection Rates
#'
#' Estimates the probability of at least one false positive under a global null
#' for permutation-based methods, comparing:
#' - default voxelwise sign-flip calibration, and
#' - subject-level sign-flip calibration via \code{null_fun}.
#'
#' Under the complete null, this "any discovery" rate coincides with the FDR
#' metric for \code{cluster_fdr()} and is directly comparable to target level.
#'
#' @param n_sims Number of independent null simulations.
#' @param dims 3D volume dimensions.
#' @param n_subj Number of subject maps used to construct the group map.
#' @param n_perm Number of permutations per method call.
#' @param alpha Target level used for TFCE (\code{alpha}) and cluster-FDR
#'   (\code{q}).
#' @param cluster_thresh Cluster-forming threshold for \code{cluster_fdr()}.
#' @param fwhm_vox Smoothness (voxel units) used when simulating subject maps.
#'   Scalar or length-3 vector.
#' @param tail Tail for both methods: \code{"pos"} or \code{"two"}.
#' @param two_sided_policy Two-sided policy for \code{cluster_fdr()}.
#' @param seed RNG seed.
#' @param verbose If TRUE, prints simulation progress.
#'
#' @return A list with components:
#'   \describe{
#'     \item{draws}{Per-simulation results with columns
#'       \code{sim}, \code{method}, \code{null_model}, \code{reject},
#'       \code{time_sec}.}
#'     \item{summary}{Aggregated rates and timing with columns
#'       \code{method}, \code{null_model}, \code{n_sims},
#'       \code{reject_rate}, \code{mc_se}, \code{ci_low}, \code{ci_high},
#'       \code{mean_time_sec}.}
#'     \item{params}{Study parameters used.}
#'   }
#'
#' @details
#' Subject maps are simulated as smoothed Gaussian fields and standardized per
#' subject. A subject-level sign-flip null is created with
#' \code{\link{make_null_fun_subject_signflip}} and injected into
#' \code{\link{tfce_fwer}} / \code{\link{cluster_fdr}} using \code{null_fun}.
#'
#' @examples
#' \dontrun{
#' st <- run_resampling_fwer_study(
#'   n_sims = 20,
#'   dims = c(10, 10, 10),
#'   n_subj = 12,
#'   n_perm = 50,
#'   alpha = 0.1,
#'   seed = 1
#' )
#' st$summary
#' }
#'
#' @export
run_resampling_fwer_study <- function(n_sims = 50,
                                      dims = c(12, 12, 12),
                                      n_subj = 20,
                                      n_perm = 200,
                                      alpha = 0.05,
                                      cluster_thresh = 2.5,
                                      fwhm_vox = 2,
                                      tail = c("pos", "two"),
                                      two_sided_policy = c("BH_all", "split_q"),
                                      seed = 1,
                                      verbose = FALSE) {
  tail <- match.arg(tail)
  two_sided_policy <- match.arg(two_sided_policy)

  if (length(dims) != 3 || any(!is.finite(dims)) || any(dims <= 1)) {
    stop("dims must be length-3 with values > 1")
  }
  dims <- as.integer(dims)
  if (!is.numeric(n_sims) || length(n_sims) != 1L || n_sims < 1) {
    stop("n_sims must be >= 1")
  }
  n_sims <- as.integer(n_sims)
  if (!is.numeric(n_subj) || length(n_subj) != 1L || n_subj < 2) {
    stop("n_subj must be >= 2")
  }
  n_subj <- as.integer(n_subj)
  if (!is.numeric(n_perm) || length(n_perm) != 1L || n_perm < 1) {
    stop("n_perm must be >= 1")
  }
  n_perm <- as.integer(n_perm)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1)")
  }

  rows <- vector("list", 4L * n_sims)
  row_i <- 1L

  for (sim in seq_len(n_sims)) {
    if (isTRUE(verbose)) {
      message(sprintf("Simulation %d / %d", sim, n_sims))
    }

    sim_seed <- as.integer(seed + sim * 1000L)
    set.seed(sim_seed)

    subj <- .simulate_subject_null_vols(
      n_subj = n_subj,
      dims = dims,
      fwhm_vox = fwhm_vox
    )

    subj_null <- make_null_fun_subject_signflip(subj, seed = sim_seed + 1L)

    tfce_vox <- .study_eval(
      function() {
        tfce_fwer(
          subj_null$z_vol,
          mask = subj_null$mask,
          n_perm = n_perm,
          alpha = alpha,
          allow_voxel_signflip = TRUE,
          tail = tail,
          seed = sim_seed + 11L
        )
      }
    )
    rows[[row_i]] <- data.frame(
      sim = sim,
      method = "tfce_fwer",
      null_model = "voxel_signflip",
      reject = tfce_vox$reject,
      time_sec = tfce_vox$time_sec
    )
    row_i <- row_i + 1L

    tfce_subj <- .study_eval(
      function() {
        tfce_fwer(
          subj_null$z_vol,
          mask = subj_null$mask,
          n_perm = n_perm,
          alpha = alpha,
          null_fun = subj_null$null_fun,
          tail = tail,
          seed = sim_seed + 12L
        )
      }
    )
    rows[[row_i]] <- data.frame(
      sim = sim,
      method = "tfce_fwer",
      null_model = "subject_signflip",
      reject = tfce_subj$reject,
      time_sec = tfce_subj$time_sec
    )
    row_i <- row_i + 1L

    cfdr_vox <- .study_eval(
      function() {
        cluster_fdr(
          subj_null$z_vol,
          mask = subj_null$mask,
          cluster_thresh = cluster_thresh,
          q = alpha,
          p_method = "perm",
          n_perm = n_perm,
          allow_voxel_signflip = TRUE,
          tail = tail,
          two_sided_policy = two_sided_policy,
          seed = sim_seed + 21L
        )
      }
    )
    rows[[row_i]] <- data.frame(
      sim = sim,
      method = "cluster_fdr_perm",
      null_model = "voxel_signflip",
      reject = cfdr_vox$reject,
      time_sec = cfdr_vox$time_sec
    )
    row_i <- row_i + 1L

    cfdr_subj <- .study_eval(
      function() {
        cluster_fdr(
          subj_null$z_vol,
          mask = subj_null$mask,
          cluster_thresh = cluster_thresh,
          q = alpha,
          p_method = "perm",
          n_perm = n_perm,
          null_fun = subj_null$null_fun,
          tail = tail,
          two_sided_policy = two_sided_policy,
          seed = sim_seed + 22L
        )
      }
    )
    rows[[row_i]] <- data.frame(
      sim = sim,
      method = "cluster_fdr_perm",
      null_model = "subject_signflip",
      reject = cfdr_subj$reject,
      time_sec = cfdr_subj$time_sec
    )
    row_i <- row_i + 1L
  }

  draws <- do.call(rbind, rows)
  summary <- .study_summarize(draws)

  list(
    draws = draws,
    summary = summary,
    params = list(
      n_sims = n_sims,
      dims = dims,
      n_subj = n_subj,
      n_perm = n_perm,
      alpha = alpha,
      cluster_thresh = cluster_thresh,
      fwhm_vox = fwhm_vox,
      tail = tail,
      two_sided_policy = two_sided_policy,
      seed = seed
    )
  )
}

#' @keywords internal
.simulate_subject_null_vols <- function(n_subj, dims, fwhm_vox) {
  if (length(fwhm_vox) == 1L) fwhm_vox <- rep(fwhm_vox, 3)
  fwhm_vox <- as.numeric(fwhm_vox)
  if (length(fwhm_vox) != 3L || any(!is.finite(fwhm_vox)) || any(fwhm_vox < 0)) {
    stop("fwhm_vox must be length 1 or 3 with finite non-negative entries")
  }

  sigma_field <- .fwhm_to_sigma_vox(fwhm_vox)
  sigma_kernel <- sigma_field / sqrt(2)
  dims_i <- as.integer(dims)
  n_vox <- prod(dims_i)

  vols <- vector("list", n_subj)
  for (i in seq_len(n_subj)) {
    eps <- stats::rnorm(n_vox)
    sm <- gaussian_smooth3d_cpp(eps, dims_i, sigma_xyz = sigma_kernel)
    sm <- .standardize_mask(sm)
    vols[[i]] <- array(sm, dim = dims_i)
  }
  vols
}

#' @keywords internal
.study_eval <- function(fun) {
  t0 <- proc.time()[3]
  out <- fun()
  list(
    reject = isTRUE(any(as.logical(out$sig_mask))),
    time_sec = proc.time()[3] - t0
  )
}

#' @keywords internal
.study_summarize <- function(draws) {
  groups <- split(draws, interaction(draws$method, draws$null_model, drop = TRUE))
  out <- lapply(groups, function(d) {
    p <- mean(d$reject)
    se <- sqrt(p * (1 - p) / nrow(d))
    data.frame(
      method = d$method[1],
      null_model = d$null_model[1],
      n_sims = nrow(d),
      reject_rate = p,
      mc_se = se,
      ci_low = max(0, p - 1.96 * se),
      ci_high = min(1, p + 1.96 * se),
      mean_time_sec = mean(d$time_sec)
    )
  })
  out <- do.call(rbind, out)
  out <- out[order(out$method, out$null_model), , drop = FALSE]
  rownames(out) <- NULL
  out
}
