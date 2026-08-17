#!/usr/bin/env Rscript
# Peak-input and densification probe for batched multiset paths.

args <- commandArgs(trailingOnly = TRUE)
out_csv <- if (length(args) >= 1L) args[[1]] else "inst/bench/results/threading_memory.csv"
seed <- if (length(args) >= 2L) as.integer(args[[2]]) else 20260816L
n <- if (length(args) >= 3L) as.integer(args[[3]]) else 64L
S <- if (length(args) >= 4L) as.integer(args[[4]]) else 4L

if (!file.exists("DESCRIPTION")) {
  stop("Run inst/bench/profile_threading_memory.R from the package root.", call. = FALSE)
}

suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
  library(rfugw)
})

source("inst/bench/protocol.R")
threads <- bench_pin_threads(1L)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

make_set <- function(nn, d = 3L) {
  X <- matrix(rnorm(nn * d), nn, d)
  F <- matrix(rnorm(nn * d), nn, d)
  C <- as.matrix(dist(X)); C <- C / max(C)
  list(C = C, F = F, w = rep(1 / nn, nn))
}

sets <- lapply(seq_len(S), function(i) c(make_set(n), list(id = sprintf("s%d", i))))
tpl <- c(make_set(n), list(id = "tpl"))

C_bytes <- sum(vapply(sets, function(s) as.numeric(object.size(s$C)), numeric(1))) +
  as.numeric(object.size(tpl$C))
F_bytes <- sum(vapply(sets, function(s) as.numeric(object.size(s$F)), numeric(1))) +
  as.numeric(object.size(tpl$F))
knn <- rfugw:::.knn_sparsify_structure(sets[[1]]$C, k = 8L)
knn_bytes <- as.numeric(object.size(knn)) * (S + 1L)

M_list <- lapply(sets, function(s) {
  as.matrix(dist(rbind(s$F, tpl$F)))[seq_len(n), n + seq_len(n)]
})
M_bytes <- sum(vapply(M_list, function(M) as.numeric(object.size(M)), numeric(1)))
plan_bytes <- as.numeric(object.size(matrix(0, n, n))) * S

ctrl <- list(
  subjects = sets,
  template_mode = "fixed",
  template = tpl,
  method = "fgw_entropic",
  alpha = 0.5,
  epsilon = 0.08,
  max_iter = 20L,
  sinkhorn_max_iter = 80L,
  use_cpp_batch = TRUE
)

elapsed_ms <- function(fn) {
  t0 <- proc.time()[["elapsed"]]
  fn()
  (proc.time()[["elapsed"]] - t0) * 1000
}

invisible(do.call(multialign_fit, c(ctrl, list(use_cpp_feature_fused = TRUE, n_threads = 1L))))
ms_fused <- elapsed_ms(function() {
  do.call(multialign_fit, c(ctrl, list(use_cpp_feature_fused = TRUE, n_threads = 1L)))
})
ms_split <- elapsed_ms(function() {
  do.call(multialign_fit, c(ctrl, list(use_cpp_feature_fused = FALSE, n_threads = 1L)))
})

out <- data.frame(
  n = n,
  subjects = S,
  threads = threads,
  structure_bytes = C_bytes,
  feature_bytes = F_bytes,
  knn_filled_bytes = knn_bytes,
  M_list_bytes = M_bytes,
  plan_bytes = plan_bytes,
  fused_ms = ms_fused,
  split_ms = ms_split,
  stringsAsFactors = FALSE
)
utils::write.csv(out, out_csv, row.names = FALSE)
print(out)
cat("wrote ", out_csv, "\n", sep = "")
