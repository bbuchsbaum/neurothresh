#' R wrappers for octree operations
#'
#' These functions provide R interfaces to the C++ octree splitting
#' and scoring routines.
#'
#' @name octree-wrappers
NULL

#' Build octree split info for a parent node
#'
#' Computes child assignments and caches constants needed for one-pass
#' scoring of all children.
#'
#' @param idx_parent Integer vector of 1-based mask-space indices
#' @param bbox Integer vector of length 6: (x0, x1, y0, y1, z0, z1)
#' @param x Integer vector of x coordinates (mask-space, length N)
#' @param y Integer vector of y coordinates (mask-space, length N)
#' @param z Integer vector of z coordinates (mask-space, length N)
#' @param pi_vec Numeric vector of prior weights (mask-space, length N)
#' @param min_pi_mass Minimum prior mass to keep a child (default 1e-10)
#'
#' @return A list with components:
#'   \describe{
#'     \item{m}{Number of non-empty children (0 to 8)}
#'     \item{child_id}{Integer vector of child assignments (-1 or 0..m-1)}
#'     \item{w_parent}{Numeric vector of prior weights for idx_parent}
#'     \item{log_den1}{Numeric vector of log(sum(pi)) per child}
#'     \item{inv_sqrt_den2}{Numeric vector of 1/sqrt(sum(pi^2)) per child}
#'     \item{keep8}{Integer vector indicating which of 1..8 octants are kept}
#'   }
#'
#' @export
octree_split_info <- function(idx_parent, bbox, x, y, z, pi_vec,
                              min_pi_mass = 1e-10) {
  octree_split_info_cpp(idx_parent, bbox, x, y, z, pi_vec, min_pi_mass)
}

#' Score all children in one pass through parent voxels
#'
#' Computes the omnibus score (max of U_0 and best S_kappa) for all
#' children in a single O(n) pass through the parent voxels.
#'
#' @param Z_vec Numeric vector of Z-scores (mask-space, length N)
#' @param idx_parent Integer vector of parent mask-space indices
#' @param split List returned by \code{\link{octree_split_info}}
#' @param kappa_grid Numeric vector of positive kappa values
#' @param do_abs Logical, whether to take absolute value of Z (for two-sided)
#'
#' @return Numeric vector of length m with child scores
#'
#' @export
score_children_onepass <- function(Z_vec, idx_parent, split, kappa_grid,
                                   do_abs = FALSE) {
  score_children_onepass_cpp(
    Z_vec = Z_vec,
    idx_parent = idx_parent,
    w_parent = split$w_parent,
    child_id = split$child_id,
    log_den1 = split$log_den1,
    inv_sqrt_den2 = split$inv_sqrt_den2,
    kappa_grid_pos = kappa_grid,
    do_abs = do_abs
  )
}

#' Split parent indices into child index lists
#'
#' Partitions the parent voxel indices into lists for each child,
#' for recursive descent.
#'
#' @param idx_parent Integer vector of parent mask-space indices
#' @param split List returned by \code{\link{octree_split_info}}
#'
#' @return List of m IntegerVectors, one per child
#'
#' @export
split_child_indices <- function(idx_parent, split) {
  split_indices_cpp(idx_parent, split$child_id, split$m)
}

#' Compute bounding box from mask-space indices
#'
#' @param idx Integer vector of 1-based mask-space indices
#' @param x Integer vector of x coordinates (mask-space)
#' @param y Integer vector of y coordinates (mask-space)
#' @param z Integer vector of z coordinates (mask-space)
#'
#' @return Integer vector of length 6: (x0, x1, y0, y1, z0, z1)
#'
#' @export
bbox_from_indices <- function(idx, x, y, z) {
  bbox_from_indices_cpp(idx, x, y, z)
}

#' Split indices with bounding boxes
#'
#' Splits parent indices into children and computes bounding box
#' for each child in a single pass.
#'
#' @param idx_parent Integer vector of parent mask-space indices
#' @param split List returned by \code{\link{octree_split_info}}
#' @param x Integer vector of x coordinates
#' @param y Integer vector of y coordinates
#' @param z Integer vector of z coordinates
#'
#' @return List with components:
#'   \describe{
#'     \item{idx}{List of m IntegerVectors with child indices}
#'     \item{bbox}{Integer matrix (m x 6) with child bounding boxes}
#'   }
#'
#' @export
split_indices_bbox <- function(idx_parent, split, x, y, z) {
  split_indices_bbox_cpp(idx_parent, split$child_id, split$m, x, y, z)
}

#' Split indices with bounding boxes and prior masses
#'
#' Like \code{\link{split_indices_bbox}} but also computes sum(pi)
#' and sum(pi^2) for each child for alpha allocation.
#'
#' @param idx_parent Integer vector of parent mask-space indices
#' @param split List returned by \code{\link{octree_split_info}}
#' @param x Integer vector of x coordinates
#' @param y Integer vector of y coordinates
#' @param z Integer vector of z coordinates
#' @param pi_vec Numeric vector of prior weights
#'
#' @return List with components:
#'   \describe{
#'     \item{idx}{List of m IntegerVectors with child indices}
#'     \item{bbox}{Integer matrix (m x 6) with child bounding boxes}
#'     \item{den1}{Numeric vector of sum(pi) per child}
#'     \item{den2}{Numeric vector of sum(pi^2) per child}
#'   }
#'
#' @export
split_indices_bbox_mass <- function(idx_parent, split, x, y, z, pi_vec) {
  split_indices_bbox_mass_cpp(idx_parent, split$child_id, split$m, x, y, z, pi_vec)
}

#' Recursive octree split
#'
#' Recursively partitions a region into octree children until
#' regions are smaller than min_voxels.
#'
#' @param indices Integer vector of mask-space indices
#' @param x Integer vector of x coordinates (mask-space)
#' @param y Integer vector of y coordinates (mask-space)
#' @param z Integer vector of z coordinates (mask-space)
#' @param min_voxels Minimum region size to continue splitting
#'
#' @return List of index vectors for leaf regions
#'
#' @export
octree_split <- function(indices, x, y, z, min_voxels = 8) {
  if (length(indices) < min_voxels) {
    return(list(indices))
  }

  bbox <- bbox_from_indices(indices, x, y, z)

  # Check if we can split (bbox spans more than 1 voxel in at least one dim)
  x_span <- bbox[2] - bbox[1]
  y_span <- bbox[4] - bbox[3]
  z_span <- bbox[6] - bbox[5]


  if (x_span == 0 && y_span == 0 && z_span == 0) {
    return(list(indices))
  }

  # Compute midpoints
  xm <- floor((bbox[1] + bbox[2]) / 2)
  ym <- floor((bbox[3] + bbox[4]) / 2)
  zm <- floor((bbox[5] + bbox[6]) / 2)

  # Assign to octants
  xi <- x[indices]
  yi <- y[indices]
  zi <- z[indices]

  child_id <- 1L +
    (xi > xm) * 1L +
    (yi > ym) * 2L +
    (zi > zm) * 4L

  # Split into children
  children <- split(indices, child_id)

  # Recursively split children
  result <- list()
  for (child_idx in children) {
    if (length(child_idx) >= min_voxels) {
      result <- c(result, octree_split(child_idx, x, y, z, min_voxels))
    } else if (length(child_idx) > 0) {
      result <- c(result, list(child_idx))
    }
  }

  result
}
