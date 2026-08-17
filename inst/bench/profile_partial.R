#!/usr/bin/env Rscript
# Profile exact (native vs lpSolve) and entropic partial paths separately.

args <- commandArgs(trailingOnly = TRUE)
out_csv <- if (length(args) >= 1L) args[[1]] else "inst/bench/results/profile_partial.csv"
seed <- if (length(args) >= 2L) as.integer(args[[2]]) else 20260816L
reps <- if (length(args) >= 3L) as.integer(args[[3]]) else 2L

if (!file.exists("DESCRIPTION")) {
  stop("Run inst/bench/profile_partial.R from the package root.", call. = FALSE)
}

suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
  library(rfugw)
})

source("inst/bench/protocol.R")
threads <- bench_pin_threads(1L)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
has_lp <- requireNamespace("lpSolve", quietly = TRUE)

median_ms <- function(fn, reps) {
  fn()
  times <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    fn()
    times[[i]] <- (proc.time()[["elapsed"]] - t0) * 1000
  }
  stats::median(times)
}

exact_sizes <- c(16L, 24L, 32L)
entropic_sizes <- c(48L, 64L, 96L)
rows <- list()
k <- 1L
add <- function(n, method, path, ms) {
  rows[[k]] <<- data.frame(n = n, method = method, path = path, median_ms = ms,
                           stringsAsFactors = FALSE)
  k <<- k + 1L
}

for (n in exact_sizes) {
  d <- bench_make_problem("fgw", n, seed)
  add(n, "partial_gw", "cpp_transport", median_ms(function() {
    partial_gromov_wasserstein(d$C1, d$C2, p = d$p, q = d$q, m = 0.7,
                               numItermax = 30L, tol = 1e-6, log = TRUE)
  }, reps))
  if (has_lp) {
    add(n, "partial_gw", "lp_matrix", median_ms(function() {
      partial_gromov_wasserstein(d$C1, d$C2, p = d$p, q = d$q, m = 0.7,
                                 numItermax = 30L, tol = 1e-6, log = TRUE,
                                 lp_solver = "lp_matrix")
    }, reps))
  }
}

for (n in entropic_sizes) {
  d <- bench_make_problem("fgw", n, seed)
  add(n, "entropic_partial_gw", "cpp_outer", median_ms(function() {
    entropic_partial_gromov_wasserstein(d$C1, d$C2, p = d$p, q = d$q, m = 0.7,
                                        reg = 0.2, numItermax = 40L, tol = 1e-6,
                                        log = TRUE)
  }, reps))
  add(n, "entropic_partial_gw", "r_outer", median_ms(function() {
    G0 <- (d$p %o% d$q) * (0.7 / (sum(d$p) * sum(d$q)))
    rfugw:::.partial_fgw_entropic_core(
      matrix(0, n, n), d$C1, d$C2, d$p, d$q, 0.7, 0.2, 1, G0,
      40L, 1e-6, TRUE, 200L, 1e-12, FALSE, 2L
    )
  }, reps))
  add(n, "entropic_partial_gw", "cpp_outer_zero_M", median_ms(function() {
    rfugw:::cpp_partial_fgw_entropic_square(
      matrix(0, n, n), d$C1, d$C2, d$p, d$q, 0.7, 0.2, 1, TRUE,
      matrix(numeric(0), 0, 0), 40L, 1e-6, 200L, 1e-12, 2L
    )
  }, reps))
}

df <- do.call(rbind, rows)
utils::write.csv(df, out_csv, row.names = FALSE)
cat("Wrote", out_csv, "threads=", threads, "\n")
print(df, row.names = FALSE)
