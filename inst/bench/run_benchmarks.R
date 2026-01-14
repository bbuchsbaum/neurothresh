#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(neurothresh))

source("inst/bench/benchmark_harness.R")

res <- benchmark_suite(n_sims = 10, n_perm = 100, dims = c(20, 20, 20), alpha = 0.05, fwhm_mm = 6)

summary <- aggregate(cbind(reject, time_sec) ~ regime + method, data = res, FUN = mean)
print(summary)

