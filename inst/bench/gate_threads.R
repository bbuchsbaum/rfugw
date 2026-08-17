#!/usr/bin/env Rscript
# Nightly 1-vs-N thread equivalence smoke for the batched kernel.

args <- commandArgs(trailingOnly = TRUE)
out_csv <- if (length(args) >= 1L) args[[1]] else "inst/bench/results/current/threads.csv"
seed <- if (length(args) >= 2L) as.integer(args[[2]]) else 20260816L

if (!file.exists("DESCRIPTION")) {
  message("INFRA: run inst/bench/gate_threads.R from the package root.")
  quit(status = 2L)
}

suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
  library(rfugw)
})

make_set <- function(n, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * 3L), n, 3L)
  F <- matrix(rnorm(n * 2L), n, 2L)
  C <- as.matrix(dist(X)); C <- C / max(C)
  list(C = C, F = F, w = rep(1 / n, n))
}

s1 <- c(make_set(10L, seed + 1L), list(id = "s1"))
s2 <- c(make_set(9L, seed + 2L), list(id = "s2"))
s3 <- c(make_set(8L, seed + 3L), list(id = "s3"))
tpl <- c(make_set(8L, seed + 4L), list(id = "tpl"))
ctrl <- list(
  subjects = list(s1, s2, s3),
  template_mode = "fixed",
  template = tpl,
  method = "fgw_entropic",
  alpha = 0.5,
  epsilon = 0.08,
  max_iter = 40L,
  sinkhorn_max_iter = 120L,
  use_cpp_batch = TRUE
)

one <- tryCatch(do.call(multialign_fit, c(ctrl, list(n_threads = 1L))), error = function(e) e)
two <- tryCatch(do.call(multialign_fit, c(ctrl, list(n_threads = 2L))), error = function(e) e)
if (inherits(one, "error") || inherits(two, "error")) {
  message("SOLVER: threaded multialign failed: ",
          if (inherits(one, "error")) one$message else two$message)
  quit(status = 1L)
}

gap <- abs(one$objective_total - two$objective_total)
ok <- is.finite(gap) && gap <= 1e-8
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  data.frame(
    seed = seed,
    used_threads_1 = one$used_threads,
    used_threads_2 = two$used_threads,
    objective_1 = one$objective_total,
    objective_2 = two$objective_total,
    abs_gap = gap,
    ok = ok,
    stringsAsFactors = FALSE
  ),
  out_csv,
  row.names = FALSE
)
cat("Wrote ", out_csv, " gap=", gap, "\n", sep = "")
if (!ok) {
  message("SOLVER: 1-thread and 2-thread objectives differ by ", gap)
  quit(status = 1L)
}
