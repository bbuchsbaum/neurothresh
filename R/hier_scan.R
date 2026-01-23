#' Hierarchical LR-MFT Analysis
#'
#' Main entry point for hierarchical likelihood-ratio matched-filter
#' thresholding with stepdown inference.
#'
#' @param z_vol A NeuroVol object containing Z-scores (or other stat type)
#' @param prior_vol Optional NeuroVol of prior weights. If NULL, uniform
#'   priors are used.
#' @param parcels Optional ClusteredNeuroVol from neuroatlas for parcel-based
#'   initialization. If NULL, starts from whole brain with octree descent.
#' @param mask Optional LogicalNeuroVol defining analysis domain. If NULL,
#'   uses all finite voxels in z_vol.
#' @param alpha Significance level for FWER control (default 0.05)
#' @param kappa Temperature parameter(s) for soft-max. Can be a single value
#'   or vector of values to search over.
#' @param n_perm Number of permutations for null distribution (default 1000)
#' @param method Inference method: "stepdown" (Westfall-Young) or
#'   "hierarchical" (alpha-spending only)
#' @param gamma Fraction of alpha for testing at each hierarchical level
#' @param gamma_root Fraction of alpha for root/parcel level testing
#'   (used for stepdown only; for hierarchical mode this is used as \code{gamma})
#' @param min_voxels Minimum region size for octree descent
#' @param min_alpha Minimum alpha to continue hierarchical testing
#' @param seed Random seed for reproducibility
#' @param stat_type Input statistic type: "Z", "t", or "neglog10p"
#' @param df Degrees of freedom (required if stat_type = "t")
#' @param two_sided Logical, whether to perform two-sided inference
#' @param prior_eta Mixing weight for prior smoothing with uniform mass.
#'   Values in \[0, 1\], where 0 is uniform and 1 is the raw prior.
#' @param parallel Logical, whether to parallelize permutation loops using
#'   future.apply (optional).
#' @param n_workers Optional integer number of workers for future plan.
#'   Ignored unless \code{parallel = TRUE} and the future package is installed.
#'
#' @return An S3 object of class "neurothresh_result" containing:
#'   \describe{
#'     \item{significant_regions}{List of significant region index vectors}
#'     \item{tree}{Hierarchical tree structure with scores and p-values}
#'     \item{parcel_results}{Results from parcel-level testing (if used)}
#'     \item{method}{Method used}
#'     \item{params}{Parameters used}
#'   }
#'
#' @details
#' The analysis proceeds as follows:
#' 1. Canonicalize input to Z-scores
#' 2. If parcels provided, test parcels with Westfall-Young step-down
#' 3. For significant parcels (or whole brain), descend with octree
#' 4. At each level, test children with step-down, spend alpha budget
#' 5. Collect all significant regions
#'
#' @examples
#' \dontrun{
#' library(neuroim2)
#' library(neuroatlas)
#'
#' # Load statistical map
#' z_map <- read_vol("zstat1.nii.gz")
#'
#' # Run with default settings
#' result <- hier_scan(z_map, alpha = 0.05, n_perm = 1000)
#'
#' # With parcellation priors
#' atlas <- get_schaefer_atlas(parcels = "200", networks = "7")
#' result <- hier_scan(z_map, parcels = atlas$atlas)
#'
#' # Summary and visualization
#' summary(result)
#' plot(result, z_map)
#' }
#'
#' @export
hier_scan <- function(z_vol,
                      prior_vol = NULL,
                      parcels = NULL,
                      mask = NULL,
                      alpha = 0.05,
                      kappa = c(0.5, 1.0, 2.0),
                      n_perm = 1000,
                      method = c("stepdown", "hierarchical"),
                      gamma = 0.5,
                      gamma_root = 0.5,
                      min_voxels = 8,
                      min_alpha = 1e-6,
                      seed = NULL,
                      stat_type = "Z",
                      df = NULL,
                      two_sided = FALSE,
                      prior_eta = 0.9,
                      parallel = FALSE,
                      n_workers = NULL) {

  method <- match.arg(method)

  if (!is.null(seed)) set.seed(seed)

  parallel <- isTRUE(parallel)
  use_future <- FALSE
  if (parallel) {
    if (requireNamespace("future.apply", quietly = TRUE)) {
      use_future <- TRUE
    } else {
      warning("parallel=TRUE requested but future.apply is not installed; running serial.")
    }
  }
  if (use_future) {
    if (requireNamespace("future", quietly = TRUE)) {
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      workers <- if (is.null(n_workers)) future::availableCores() else n_workers
      future::plan(future::multisession, workers = workers)
    } else if (!is.null(n_workers)) {
      warning("n_workers ignored because future is not installed.")
    }
  }

  # --- (A) Canonicalize input ---
  if (stat_type != "Z") {
    canon <- canonicalize_stat(z_vol, stat_type = stat_type, df = df,
                               tail = if (two_sided) "two" else "pos")
    z_vol <- canon$Zeq
  }

  # --- (B) Set up mask ---
  if (is.null(mask)) {
    mask <- is.finite(as.numeric(z_vol))
  } else {
    mask <- as.logical(mask)
  }
  mask_idx <- which(mask)
  n_mask <- length(mask_idx)

  if (n_mask == 0) {
    stop("No valid voxels in mask")
  }

  # --- (C) Extract mask-space vectors ---
  z_vec <- as.numeric(z_vol)[mask_idx]
  if (two_sided) {
    z_vec <- abs(z_vec)
  }

  # Prior weights
  if (is.null(prior_vol)) {
    pi_vec <- rep(1.0, n_mask)
  } else {
    pi_vec <- as.numeric(prior_vol)[mask_idx]
  }
  pi_vec <- .prep_prior(pi_vec, eta = prior_eta)

  # Coordinates
  dims <- dim(z_vol)
  grid <- .index_to_grid(mask_idx, dims)
  x <- as.integer(grid[, 1])
  y <- as.integer(grid[, 2])
  z_coord <- as.integer(grid[, 3])

  # --- (D) Create permutation function ---
  perm_fun <- make_perm_fun(z_vec, n_perm, seed = seed)

  # --- (E) Parcel-level testing (if parcels provided) ---
  parcel_results <- NULL
  significant_regions <- list()

  if (!is.null(parcels)) {
    lab_vec <- as.integer(parcels)[mask_idx]
    lab_vec[is.na(lab_vec)] <- 0

    # Remove background (label 0)
    valid_labels <- lab_vec > 0

    if (sum(valid_labels) == 0) {
      warning("No valid parcel labels in mask")
    } else if (method == "stepdown") {
      # Test parcels
      parcel_results <- test_parcels(
        z_vec[valid_labels], pi_vec[valid_labels], lab_vec[valid_labels],
        n_perm = n_perm, kappa_grid = kappa,
        alpha = gamma_root * alpha, seed = seed,
        parallel = use_future
      )

      # Descend into significant parcels
      if (length(parcel_results$sig_parcels) > 0) {
        alpha_desc <- alpha - gamma_root * alpha
        parcel_mass <- parcel_results$den1[parcel_results$stepdown$rejected]
        w <- parcel_mass / sum(parcel_mass)

        for (k in seq_along(parcel_results$sig_parcels)) {
          p <- parcel_results$sig_parcels[k]
          p_str <- as.character(p)
          idx_parcel <- parcel_results$parcel_indices[[p_str]]

          # Map back to full mask-space
          idx_full <- which(lab_vec == p)

          if (length(idx_full) >= min_voxels) {
            bbox <- bbox_from_indices_cpp(idx_full, x, y, z_coord)

            hits <- hier_descend(
              idx_full, bbox,
              alpha_budget = alpha_desc * w[k],
              z_vec, pi_vec, x, y, z_coord,
              perm_fun, n_perm, kappa,
              gamma, min_alpha, min_voxels,
              parallel = use_future
            )

            significant_regions <- c(significant_regions, hits)
          }
        }
      }
    } else {
      parcel_ids <- sort(unique(lab_vec[valid_labels]))
      parcel_indices <- lapply(parcel_ids, function(p) which(lab_vec == p))
      names(parcel_indices) <- as.character(parcel_ids)
      parcel_mass <- vapply(parcel_indices, function(idx) sum(pi_vec[idx]),
                            numeric(1))
      parcel_mass <- parcel_mass / sum(parcel_mass)

      gamma_use <- gamma_root

      for (k in seq_along(parcel_ids)) {
        idx_parcel <- parcel_indices[[k]]
        if (length(idx_parcel) < min_voxels) next
        bbox <- bbox_from_indices_cpp(idx_parcel, x, y, z_coord)

        hits <- hier_descend_alpha(
          idx_parcel, bbox,
          alpha_budget = alpha * parcel_mass[k],
          z_vec, pi_vec, x, y, z_coord,
          perm_fun, n_perm, kappa,
          gamma = gamma_use, min_alpha, min_voxels,
          parallel = use_future
        )

        significant_regions <- c(significant_regions, hits)
      }
    }
  } else {
    # --- (F) Whole-brain octree descent ---
    idx_all <- seq_len(n_mask)
    bbox <- bbox_from_indices_cpp(idx_all, x, y, z_coord)

    if (method == "stepdown") {
      significant_regions <- hier_descend(
        idx_all, bbox,
        alpha_budget = alpha,
        z_vec, pi_vec, x, y, z_coord,
        perm_fun, n_perm, kappa,
        gamma, min_alpha, min_voxels,
        parallel = use_future
      )
    } else {
      significant_regions <- hier_descend_alpha(
        idx_all, bbox,
        alpha_budget = alpha,
        z_vec, pi_vec, x, y, z_coord,
        perm_fun, n_perm, kappa,
        gamma = gamma_root, min_alpha, min_voxels,
        parallel = use_future
      )
    }
  }

  # --- (G) Convert mask-space indices back to volume indices ---
  for (i in seq_along(significant_regions)) {
    significant_regions[[i]]$vol_indices <- mask_idx[significant_regions[[i]]$indices]
  }

  # --- (H) Build result object ---
  result <- list(
    significant_regions = significant_regions,
    parcel_results = parcel_results,
    n_significant = length(significant_regions),
    method = method,
    params = list(
      alpha = alpha,
      kappa = kappa,
      n_perm = n_perm,
      gamma = gamma,
      gamma_root = gamma_root,
      min_voxels = min_voxels,
      two_sided = two_sided,
      parallel = parallel,
      n_workers = n_workers
    ),
    mask_idx = mask_idx,
    dims = dims
  )

  class(result) <- "neurothresh_result"
  result
}

#' @keywords internal
.prep_prior <- function(pi_vec, eta = 0.9) {
  pi <- pi_vec
  pi[is.na(pi)] <- 0
  pi[pi < 0] <- 0

  s <- sum(pi)
  if (s <= 0) {
    stop("Prior has zero mass in mask")
  }
  pi <- pi / s

  n <- length(pi)
  eta <- max(0, min(1, eta))
  pi <- (1 - eta) * (1 / n) + eta * pi
  pi / sum(pi)
}

#' @keywords internal
.index_to_grid <- function(idx, dims) {
  # Convert 1-based linear indices to (i, j, k) coordinates
  idx0 <- idx - 1  # 0-based
  i <- (idx0 %% dims[1]) + 1
  j <- ((idx0 %/% dims[1]) %% dims[2]) + 1
  k <- (idx0 %/% (dims[1] * dims[2])) + 1
  cbind(i, j, k)
}

#' Print method for neurothresh_result
#'
#' @param x A neurothresh_result object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.neurothresh_result <- function(x, ...) {
  cat("neurothresh hierarchical analysis result\n")
  cat("----------------------------------------\n")
  cat("Method:", x$method, "\n")
  cat("Alpha:", x$params$alpha, "\n")
  cat("Permutations:", x$params$n_perm, "\n")
  cat("Kappa values:", paste(x$params$kappa, collapse = ", "), "\n")
  cat("Significant regions:", x$n_significant, "\n")
  if (!is.null(x$parcel_results)) {
    cat("Significant parcels:",
        length(x$parcel_results$sig_parcels), "\n")
  }
  invisible(x)
}

#' Summary method for neurothresh_result
#'
#' @param object A neurothresh_result object
#' @param ... Additional arguments (ignored)
#'
#' @export
summary.neurothresh_result <- function(object, ...) {
  cat("neurothresh Hierarchical Analysis Summary\n")
  cat("=========================================\n\n")

  cat("Parameters:\n")
  cat("  FWER alpha:", object$params$alpha, "\n")
  cat("  Permutations:", object$params$n_perm, "\n")
  cat("  Kappa grid:", paste(object$params$kappa, collapse = ", "), "\n")
  cat("  Method:", object$method, "\n")
  cat("  Two-sided:", object$params$two_sided, "\n")
  cat("\n")

  cat("Results:\n")
  cat("  Significant regions:", object$n_significant, "\n")

  if (!is.null(object$parcel_results)) {
    cat("  Significant parcels:",
        length(object$parcel_results$sig_parcels), "\n")
    if (length(object$parcel_results$sig_parcels) > 0) {
      cat("  Parcel IDs:", paste(object$parcel_results$sig_parcels,
                                 collapse = ", "), "\n")
    }
  }

  if (object$n_significant > 0) {
    cat("\n")
    cat("Top regions by score:\n")
    scores <- sapply(object$significant_regions, function(r) r$score)
    p_adj <- sapply(object$significant_regions, function(r) r$p_adj)
    sizes <- sapply(object$significant_regions, function(r) length(r$indices))
    ord <- order(scores, decreasing = TRUE)

    n_show <- min(10, length(ord))
    for (i in seq_len(n_show)) {
      idx <- ord[i]
      cat(sprintf("  %2d. Score=%.3f, p_adj=%.4f, size=%d voxels\n",
                  i, scores[idx], p_adj[idx], sizes[idx]))
    }
  }

  invisible(object)
}
