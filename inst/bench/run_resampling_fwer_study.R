#!/usr/bin/env Rscript

if (!requireNamespace("neurothresh", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Either install 'neurothresh' or install 'devtools' for load_all('.')")
  }
  suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
} else {
  suppressPackageStartupMessages(library(neurothresh))
}

args <- commandArgs(trailingOnly = TRUE)

# Defaults chosen to keep runtime practical while reducing Monte Carlo noise.
n_sims <- if (length(args) >= 1) as.integer(args[[1]]) else 40L
n_perm <- if (length(args) >= 2) as.integer(args[[2]]) else 100L

study <- run_resampling_fwer_study(
  n_sims = n_sims,
  dims = c(10, 10, 10),
  n_subj = 16,
  n_perm = n_perm,
  alpha = 0.05,
  cluster_thresh = 2.5,
  fwhm_vox = 2,
  tail = "pos",
  seed = 2026
)

cat("\nResampling FWER Study Summary\n")
print(study$summary)

if (length(args) >= 3) {
  out_csv <- args[[3]]
  utils::write.csv(study$draws, file = out_csv, row.names = FALSE)
  cat("\nSaved per-simulation draws to:", out_csv, "\n")
}
