#!/usr/bin/env Rscript
# Profile semi-relaxed GW/FGW hot paths at representative sizes.
# Writes inst/bench/results/profile_semirelaxed.csv and prints a summary.

args <- commandArgs(trailingOnly = TRUE)
out_csv <- if (length(args) >= 1L) args[[1]] else "inst/bench/results/profile_semirelaxed.csv"
seed <- if (length(args) >= 2L) as.integer(args[[2]]) else 20260816L
reps <- if (length(args) >= 3L) as.integer(args[[3]]) else 3L

if (!file.exists("DESCRIPTION")) {
  stop("Run inst/bench/profile_semirelaxed.R from the package root.", call. = FALSE)
}

suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
  library(rfugw)
})

source("inst/bench/protocol.R")
threads <- bench_pin_threads(1L)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

empty_M <- function() matrix(numeric(0), nrow = 0, ncol = 0)

median_ms <- function(fn, reps) {
  times <- numeric(reps)
  fn()
  for (i in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    fn()
    times[[i]] <- (proc.time()[["elapsed"]] - t0) * 1000
  }
  stats::median(times)
}

sizes <- c(32L, 64L, 96L)
rows <- list()
k <- 1L

add <- function(n, method, path, ms, extra = "") {
  rows[[k]] <<- data.frame(
    n = n, method = method, path = path, median_ms = ms, note = extra,
    stringsAsFactors = FALSE
  )
  k <<- k + 1L
}

for (n in sizes) {
  d <- bench_make_problem("fgw", n, seed)
  C1_as <- d$C1
  C1_as[1L, 2L] <- C1_as[1L, 2L] + 0.2
  p <- d$p
  M0 <- matrix(0, n, n)

  add(n, "srfgw", "r_core_symmetric", median_ms(function() {
    rfugw:::.semirelaxed_fgw_cg_core(
      d$M, d$C1, d$C2, p, 0.5, TRUE, NULL, 40L, 1e-6, 1e-6
    )
  }, reps))
  add(n, "srfgw", "cpp_exact_symmetric", median_ms(function() {
    rfugw:::cpp_semirelaxed_fgw_exact_square(
      d$M, d$C1, d$C2, p, 0.5, TRUE, empty_M(), 40L, 1e-6, 1e-6
    )
  }, reps))
  add(n, "srfgw", "cpp_fast_double", median_ms(function() {
    rfugw:::cpp_semirelaxed_fgw_cg_square_fast(
      d$M, d$C1, d$C2, p, 0.5, empty_M(), 40L, 1e-6, 1e-6, FALSE, FALSE
    )
  }, reps))
  add(n, "srfgw", "cpp_fast_mixed", median_ms(function() {
    rfugw:::cpp_semirelaxed_fgw_cg_square_fast(
      d$M, d$C1, d$C2, p, 0.5, empty_M(), 40L, 1e-6, 1e-6, FALSE, TRUE
    )
  }, reps))
  add(n, "srgw", "r_core_asymmetric", median_ms(function() {
    rfugw:::.semirelaxed_fgw_cg_core(
      M0, C1_as, d$C2, p, 1, FALSE, NULL, 40L, 1e-6, 1e-6
    )
  }, reps))
  add(n, "srgw", "cpp_exact_asymmetric", median_ms(function() {
    rfugw:::cpp_semirelaxed_fgw_exact_square(
      empty_M(), C1_as, d$C2, p, 1, FALSE, empty_M(), 40L, 1e-6, 1e-6
    )
  }, reps))
  add(n, "srgw", "cpp_fast_empty_M", median_ms(function() {
    rfugw:::cpp_semirelaxed_fgw_cg_square_fast(
      empty_M(), d$C1, d$C2, p, 1, empty_M(), 40L, 1e-6, 1e-6, FALSE, FALSE
    )
  }, reps))
  add(n, "srgw", "cpp_fast_zero_M", median_ms(function() {
    rfugw:::cpp_semirelaxed_fgw_cg_square_fast(
      M0, d$C1, d$C2, p, 1, empty_M(), 40L, 1e-6, 1e-6, FALSE, FALSE
    )
  }, reps))
}

df <- do.call(rbind, rows)
utils::write.csv(df, out_csv, row.names = FALSE)
cat("Wrote", out_csv, "threads=", threads, "\n")
print(df, row.names = FALSE)
