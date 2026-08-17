#!/usr/bin/env Rscript
# Profile UCOOT sample vs feature phases, and C++ vs R BCD.

args <- commandArgs(trailingOnly = TRUE)
out_csv <- if (length(args) >= 1L) args[[1]] else "inst/bench/results/profile_ucoot.csv"
seed <- if (length(args) >= 2L) as.integer(args[[2]]) else 20260816L
reps <- if (length(args) >= 3L) as.integer(args[[3]]) else 2L

if (!file.exists("DESCRIPTION")) {
  stop("Run inst/bench/profile_ucoot.R from the package root.", call. = FALSE)
}

suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
  library(rfugw)
})

source("inst/bench/protocol.R")
threads <- bench_pin_threads(1L)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

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

sizes <- c(16L, 24L, 32L)
rows <- list()
k <- 1L
add <- function(n, method, path, ms, feat_ms = NA_real_, samp_ms = NA_real_) {
  rows[[k]] <<- data.frame(
    n = n, method = method, path = path, median_ms = ms,
    feat_ms = feat_ms, samp_ms = samp_ms, stringsAsFactors = FALSE
  )
  k <<- k + 1L
}

for (n in sizes) {
  d <- bench_make_problem("ucoot", n, seed)
  add(n, "ucoot_independent", "cpp", median_ms(function() {
    unbalanced_co_optimal_transport(
      d$X, d$Y, wx_samp = d$wx_samp, wy_samp = d$wy_samp,
      wx_feat = d$wx_feat, wy_feat = d$wy_feat,
      max_iter = 30L, tol = 1e-6, max_iter_ot = 120L, log = TRUE
    )
  }, reps))
  timed <- unbalanced_co_optimal_transport(
    d$X, d$Y, wx_samp = d$wx_samp, wy_samp = d$wy_samp,
    wx_feat = d$wx_feat, wy_feat = d$wy_feat,
    max_iter = 30L, tol = 1e-6, max_iter_ot = 120L, log = TRUE
  )
  rows[[k - 1L]]$feat_ms <- timed$feat_ms %||% NA_real_
  rows[[k - 1L]]$samp_ms <- timed$samp_ms %||% NA_real_

  add(n, "ucoot_independent", "r_core", median_ms(function() {
    rfugw:::.ucoot_kl_r_core(
      d$X, d$Y, d$wx_samp, d$wx_feat, d$wy_samp, d$wy_feat,
      d$wx_samp %o% d$wy_samp, d$wx_feat %o% d$wy_feat,
      10, 8, 0.05, 0.03, NULL, NULL,
      d$wx_samp %o% d$wy_samp, d$wx_feat %o% d$wy_feat,
      "independent", TRUE, 30L, 1e-6, 120L, 1e-7, TRUE, FALSE
    )
  }, reps))
  add(n, "across_spaces_joint", "cpp", median_ms(function() {
    fused_unbalanced_across_spaces_divergence(
      d$X, d$Y, wx_samp = d$wx_samp, wy_samp = d$wy_samp,
      wx_feat = d$wx_feat, wy_feat = d$wy_feat,
      reg_type = "joint", max_iter = 30L, tol = 1e-6,
      max_iter_ot = 120L, log = TRUE
    )
  }, reps))
}

df <- do.call(rbind, rows)
utils::write.csv(df, out_csv, row.names = FALSE)
cat("Wrote", out_csv, "threads=", threads, "\n")
print(df, row.names = FALSE)
