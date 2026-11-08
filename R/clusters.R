#' Label supra-threshold clusters
#'
#' Finds the connected components of the voxels whose statistic exceeds a
#' cluster-forming threshold and summarises each one by its extent, excess
#' mass and peak. The labelling runs in compiled code (union-find over the
#' neighbour graph of the mask) and is the building block for cluster-extent
#' and cluster-mass permutation inference.
#'
#' @param stat_vol 3D statistic map (array or `NeuroVol`).
#' @param threshold Cluster-forming threshold on the scale of `stat_vol`. For
#'   `tail = "neg"` clusters are formed where `-stat_vol >= threshold`.
#' @param mask Optional logical 3D mask. Defaults to finite voxels.
#' @param tail `"pos"` (default) or `"neg"`.
#' @param connectivity `6` (default), `18` or `26`.
#'
#' @return A list with
#'   \describe{
#'     \item{labels}{Integer array over the grid: cluster number, 0 outside
#'       clusters. Clusters are numbered in order of their smallest linear
#'       voxel index.}
#'     \item{table}{Data frame with one row per cluster: `cluster_id`,
#'       `size` (voxels), `mass` (sum of `stat - threshold` over the cluster,
#'       on the tested tail's scale), `peak_value` (signed statistic at the
#'       peak), `peak_index` (linear voxel index) and the 1-based grid
#'       coordinates `peak_i`, `peak_j`, `peak_k`.}
#'   }
#' @examples
#' z <- array(0, c(6, 6, 6))
#' z[2:3, 2:3, 2:3] <- 4
#' z[5, 5, 5] <- 5
#' label_clusters(z, threshold = 3)$table
#' @export
label_clusters <- function(stat_vol, threshold, mask = NULL,
                           tail = c("pos", "neg"), connectivity = 6L) {
  tail <- match.arg(tail)
  connectivity <- .check_connectivity(connectivity)
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold)) {
    stop("threshold must be one finite number", call. = FALSE)
  }
  dims <- dim(stat_vol)
  if (length(dims) != 3L) stop("stat_vol must be 3D", call. = FALSE)
  z <- as.numeric(stat_vol)
  if (is.null(mask)) mask <- is.finite(z)
  mask_idx <- which(as.logical(mask))
  sgn <- if (tail == "neg") -1 else 1
  labels <- integer(prod(dims))
  empty <- data.frame(
    cluster_id = integer(0), size = integer(0), mass = numeric(0),
    peak_value = numeric(0), peak_index = integer(0),
    peak_i = integer(0), peak_j = integer(0), peak_k = integer(0)
  )
  if (!length(mask_idx)) {
    return(list(labels = array(labels, dim = dims), table = empty))
  }
  res <- cluster_label_cpp(sgn * z[mask_idx], as.integer(dims),
                           as.integer(mask_idx), as.numeric(threshold),
                           connectivity)
  labels[mask_idx] <- res$labels
  if (!length(res$sizes)) {
    return(list(labels = array(labels, dim = dims), table = empty))
  }
  peak_index <- mask_idx[res$peak]
  ijk <- arrayInd(peak_index, dims)
  table <- data.frame(
    cluster_id = seq_along(res$sizes),
    size = res$sizes,
    mass = res$mass,
    peak_value = z[peak_index],
    peak_index = as.integer(peak_index),
    peak_i = as.integer(ijk[, 1]),
    peak_j = as.integer(ijk[, 2]),
    peak_k = as.integer(ijk[, 3])
  )
  list(labels = array(labels, dim = dims), table = table)
}
