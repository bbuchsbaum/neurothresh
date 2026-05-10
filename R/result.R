#' Construct a neurothresh result
#'
#' Internal constructor used by every public inference method to return a
#' uniform result object. The S3 class is `"neurothresh_result"`; downstream
#' code should treat the object as opaque except for the documented slots.
#'
#' @param method Character; one of `"hier_scan"`, `"octree_scan_fwer"`,
#'   `"octree_scan_stepdown"`, `"tfce_fwer"`, `"cluster_fdr"`,
#'   `"rft_peak_fwer"`, `"rft_cluster_fwer"`, or `"lisa"`.
#' @param stat_vol The observed (canonicalized) statistic map.
#' @param mask Logical mask used.
#' @param threshold Critical value on the relevant statistic (numeric scalar).
#' @param reject Logical 3D rejection volume.
#' @param pvals Numeric 3D corrected p-value volume.
#' @param clusters Optional data frame of clusters (NULL for voxelwise-only
#'   methods).
#' @param null_max Optional vector of max-null statistics for permutation
#'   methods.
#' @param params List of method-specific parameters
#'   (alpha/q, tail, n_perm, kappa, H/E, smoothness, etc.).
#' @param call The captured `match.call()` from the public method.
#' @param extras Optional named list of method-specific extras (e.g., per-node
#'   tables for hier_scan, TFCE map for tfce_fwer). Stored under `$extras`.
#'
#' @return Object of class `"neurothresh_result"`.
#' @keywords internal
#' @aliases neurothresh_result
new_neurothresh_result <- function(method,
                                   stat_vol,
                                   mask,
                                   threshold,
                                   reject,
                                   pvals,
                                   clusters = NULL,
                                   null_max = NULL,
                                   params = list(),
                                   call = NULL,
                                   extras = list()) {
  if (!is.character(method) || length(method) != 1L) {
    stop("method must be a length-1 character")
  }
  if (!is.numeric(threshold) || length(threshold) != 1L) {
    stop("threshold must be a length-1 numeric")
  }

  structure(
    list(
      method = method,
      stat_vol = stat_vol,
      mask = mask,
      threshold = unname(threshold),
      reject = reject,
      pvals = pvals,
      clusters = clusters,
      null_max = null_max,
      params = params,
      call = call,
      extras = extras
    ),
    class = "neurothresh_result"
  )
}

#' Print a neurothresh result
#'
#' @param x A `neurothresh_result`.
#' @param ... Unused.
#' @export
print.neurothresh_result <- function(x, ...) {
  cat("<neurothresh_result>\n")
  cat("  method   :", x$method, "\n")
  if (!is.null(x$params$alpha)) cat("  alpha    :", x$params$alpha, "\n")
  if (!is.null(x$params$q))     cat("  q        :", x$params$q, "\n")
  if (!is.null(x$params$tail))  cat("  tail     :", x$params$tail, "\n")
  cat("  threshold:", format(x$threshold, digits = 4), "\n")

  n_reject <- if (is.null(x$reject)) NA_integer_ else sum(as.logical(x$reject), na.rm = TRUE)
  n_mask   <- if (is.null(x$mask))   NA_integer_ else sum(as.logical(x$mask),   na.rm = TRUE)
  if (!is.na(n_reject) && !is.na(n_mask)) {
    cat(sprintf("  rejected : %d / %d voxels (%.2f%%)\n",
                n_reject, n_mask, 100 * n_reject / max(n_mask, 1L)))
  }

  if (!is.null(x$clusters) && nrow(x$clusters) > 0L) {
    cat("  clusters :", nrow(x$clusters), "\n")
  }
  if (!is.null(x$null_max)) {
    cat("  n_perm   :", length(x$null_max), "\n")
  }
  invisible(x)
}

#' Summarize a neurothresh result
#'
#' @param object A `neurothresh_result`.
#' @param ... Unused.
#' @return The result object invisibly; prints a richer summary including
#'   cluster table head when present.
#' @export
summary.neurothresh_result <- function(object, ...) {
  print(object)
  if (!is.null(object$clusters) && nrow(object$clusters) > 0L) {
    cat("\nClusters (head):\n")
    print(utils::head(object$clusters, 10))
  }
  if (!is.null(object$params) && length(object$params) > 0L) {
    cat("\nParameters:\n")
    nm <- names(object$params)
    for (k in seq_along(nm)) {
      v <- object$params[[k]]
      vs <- if (is.numeric(v) && length(v) <= 5L) paste(format(v, digits = 4), collapse = ", ")
            else if (is.character(v) && length(v) <= 5L) paste(v, collapse = ", ")
            else paste0("<", class(v)[1L], "[", length(v), "]>")
      cat(sprintf("  %-12s : %s\n", nm[k], vs))
    }
  }
  invisible(object)
}

#' Coerce a neurothresh result's clusters to a data frame
#'
#' @param x A `neurothresh_result`.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Unused.
#' @return Data frame of clusters, or an empty data frame if none.
#' @export
as.data.frame.neurothresh_result <- function(x, row.names = NULL, optional = FALSE, ...) {
  if (is.null(x$clusters)) {
    return(data.frame())
  }
  as.data.frame(x$clusters)
}
