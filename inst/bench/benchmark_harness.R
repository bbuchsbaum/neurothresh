# Benchmark harness for neurothresh (synthetic + method comparisons)
#
# Intended usage (from repo root):
#   Rscript inst/bench/run_benchmarks.R
#
# This script is deliberately dependency-light (base R + neurothresh).

insert_cube <- function(z, x0, x1, y0, y1, z0, z1, amp = 1.5) {
  z[x0:x1, y0:y1, z0:z1] <- z[x0:x1, y0:y1, z0:z1] + amp
  z
}

make_mask <- function(dims) {
  rep(TRUE, prod(dims))
}

run_once <- function(z, dims, fwhm_mm = 6, alpha = 0.05, n_perm = 200, seed = 1) {
  # Simple baselines; note RFT expects fwhm_mm.
  out <- list()

  t0 <- proc.time()[3]
  res_oct <- neurothresh::octree_scan_fwer(z, alpha = alpha, n_perm = n_perm, null = "mc_fwhm")
  out$octree_fwer <- list(rej = nrow(res_oct$sig_nodes) > 0, time = proc.time()[3] - t0)

  t0 <- proc.time()[3]
  res_sd <- neurothresh::octree_scan_stepdown(z, alpha = alpha, n_perm = n_perm, null = "mc_fwhm")
  out$octree_stepdown <- list(rej = nrow(res_sd$selected) > 0, time = proc.time()[3] - t0)

  t0 <- proc.time()[3]
  res_tf <- neurothresh::tfce_fwer(z, n_perm = n_perm, alpha = alpha, tail = "pos")
  out$tfce <- list(rej = any(res_tf$sig_mask), time = proc.time()[3] - t0)

  t0 <- proc.time()[3]
  res_rft <- neurothresh::rft_peak_fwer(z, fwhm_mm = fwhm_mm, alpha = alpha, tail = "pos")
  out$rft_peak <- list(rej = any(res_rft$sig_mask), time = proc.time()[3] - t0)

  out
}

simulate_null <- function(dims) {
  array(rnorm(prod(dims)), dim = dims)
}

simulate_blob <- function(dims, amp = 1.5) {
  z <- simulate_null(dims)
  insert_cube(z, 9, 12, 9, 12, 9, 12, amp = amp)
}

benchmark_suite <- function(n_sims = 20,
                            dims = c(20, 20, 20),
                            alpha = 0.05,
                            n_perm = 200,
                            fwhm_mm = 6,
                            seed = 1) {
  set.seed(seed)

  methods <- c("octree_fwer", "octree_stepdown", "tfce", "rft_peak")
  rows <- list()

  for (sim in seq_len(n_sims)) {
    z0 <- simulate_null(dims)
    r0 <- run_once(z0, dims, fwhm_mm = fwhm_mm, alpha = alpha, n_perm = n_perm, seed = seed + sim)
    for (m in methods) {
      rows[[length(rows) + 1]] <- data.frame(
        sim = sim,
        regime = "null",
        method = m,
        reject = r0[[m]]$rej,
        time_sec = r0[[m]]$time
      )
    }

    z1 <- simulate_blob(dims, amp = 1.5)
    r1 <- run_once(z1, dims, fwhm_mm = fwhm_mm, alpha = alpha, n_perm = n_perm, seed = seed + 1000 + sim)
    for (m in methods) {
      rows[[length(rows) + 1]] <- data.frame(
        sim = sim,
        regime = "blob",
        method = m,
        reject = r1[[m]]$rej,
        time_sec = r1[[m]]$time
      )
    }
  }

  do.call(rbind, rows)
}

