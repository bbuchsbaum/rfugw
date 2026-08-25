#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !args[[1L]] %in% c("write", "read")) {
  stop("Usage: verify_sinkhorn_state_roundtrip.R write|read STATE.rds", call. = FALSE)
}
suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
  library(rfugw)
})

set.seed(7)
M <- matrix(runif(48), 8, 6)
p <- seq_len(8); p <- p / sum(p)
q <- rev(seq_len(6)); q <- q / sum(q)

if (identical(args[[1L]], "write")) {
  fit <- ot_sinkhorn(
    M, p, q, epsilon = 0.5, method = "auto",
    max_iter = 5000L, tol = 1e-9
  )
  stopifnot(fit$converged)
  saveRDS(
    list(state = fit$dual_state, plan = fit$plan, value = fit$ot_dist),
    args[[2L]]
  )
} else {
  saved <- readRDS(args[[2L]])
  fit <- ot_sinkhorn(
    M, p, q, epsilon = 0.5, method = "auto",
    max_iter = 5000L, tol = 1e-9, init_duals = saved$state
  )
  stopifnot(
    fit$converged,
    max(abs(fit$plan - saved$plan)) <= 2e-9,
    abs(fit$ot_dist - saved$value) <= 2e-9
  )
}
