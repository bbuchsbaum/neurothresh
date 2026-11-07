# Benchmark: compiled TFCE versus the pure-R reference.
#
# Rscript inst/bench/bench_tfce.R
#
# 1. A 40 x 48 x 40 grid (ellipsoid mask) where the R reference is still
#    feasible: both implementations, plus their maximum difference.
# 2. A 2 mm MNI-sized grid (91 x 109 x 91, ellipsoid mask of ~244k voxels):
#    compiled TFCE (one map) and a 64-draw batch through tfce_max_batch_cpp.
suppressPackageStartupMessages(library(neurothresh))
ns <- asNamespace("neurothresh")

ellipsoid_mask <- function(dims, frac = 0.85) {
  g <- expand.grid(i = seq_len(dims[1]), j = seq_len(dims[2]), k = seq_len(dims[3]))
  c0 <- (dims + 1) / 2
  r <- frac * (dims - 1) / 2
  inside <- ((g$i - c0[1]) / r[1])^2 + ((g$j - c0[2]) / r[2])^2 +
    ((g$k - c0[3]) / r[3])^2 <= 1
  array(inside, dims)
}

smooth_z <- function(dims, sd_vox = 2) {
  z <- ns$gaussian_smooth3d_cpp(rnorm(prod(dims)), as.integer(dims),
                                rep(sd_vox, 3), 3L)
  array(z / stats::sd(z), dims)
}

set.seed(1)
dims <- c(40L, 48L, 40L)
mask <- ellipsoid_mask(dims)
z <- smooth_z(dims)
cat(sprintf("Grid %s, %d in-mask voxels\n", paste(dims, collapse = "x"), sum(mask)))
t_r <- system.time(a <- ns$.tfce_from_vec_r(as.numeric(z), dims, mask, 2, 0.5, 0.1, "pos"))
t_c <- system.time(b <- ns$.tfce_from_vec(as.numeric(z), dims, mask, 2, 0.5, 0.1, "pos"))
cat(sprintf("  R reference: %.2f s   compiled: %.4f s   speed-up: %.0fx   max |diff|: %.2e\n",
            t_r[["elapsed"]], t_c[["elapsed"]], t_r[["elapsed"]] / max(t_c[["elapsed"]], 1e-4),
            max(abs(a - b))))

dims <- c(91L, 109L, 91L)
mask <- ellipsoid_mask(dims, frac = 0.81)
mask_idx <- which(mask)
cat(sprintf("Grid %s, %d in-mask voxels\n", paste(dims, collapse = "x"), length(mask_idx)))
z <- smooth_z(dims)
t1 <- system.time(for (i in 1:5) ns$tfce_masked_cpp(z[mask_idx], dims, mask_idx, 2, 0.5, 0.1, 6L, 0L))
cat(sprintf("  two-sided TFCE, one map: %.3f s\n", t1[["elapsed"]] / 5))
Z <- vapply(1:64, function(b) smooth_z(dims)[mask_idx], numeric(length(mask_idx)))
t2 <- system.time(ns$tfce_max_batch_cpp(Z, dims, mask_idx, 2, 0.5, 0.1, 6L, 0L))
cat(sprintf("  two-sided max-TFCE, 64-draw batch: %.2f s (%.3f s per draw)\n",
            t2[["elapsed"]], t2[["elapsed"]] / 64))
