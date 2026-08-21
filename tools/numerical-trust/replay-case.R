#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !file.exists(args[[1L]])) {
  stop(
    "Usage: Rscript tools/numerical-trust/replay-case.R CASE.rds",
    call. = FALSE
  )
}

if (file.exists("DESCRIPTION") && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(rfugw)
}

case <- readRDS(args[[1L]])
if (!identical(case$family, "transport_simplex")) {
  stop("Unsupported counterexample family: ", case$family, call. = FALSE)
}

out <- ot_emd(case$M, case$p, case$q, max_iter = 20000L, tol = 1e-11)
reduced <- case$M - outer(
  as.numeric(out$source_potential),
  as.numeric(out$target_potential),
  "+"
)
evidence <- list(
  family = case$family,
  id = case$id,
  seed = case$seed,
  mode = case$mode,
  dimensions = dim(case$M),
  termination_reason = out$termination_reason,
  lp_ok = out$lp_ok,
  row_residual = max(abs(rowSums(out$plan) - case$p)),
  col_residual = max(abs(colSums(out$plan) - case$q)),
  min_reduced_cost = min(reduced),
  primal_objective = sum(case$M * out$plan),
  dual_objective = sum(case$p * out$source_potential) +
    sum(case$q * out$target_potential),
  failed_checks_at_capture = case$failed_checks
)
print(evidence)

ok <- identical(out$termination_reason, "optimal") && isTRUE(out$lp_ok) &&
  evidence$row_residual <= out$feasibility_tolerance &&
  evidence$col_residual <= out$feasibility_tolerance &&
  evidence$min_reduced_cost >= -out$reduced_cost_tolerance &&
  abs(evidence$primal_objective - evidence$dual_objective) <=
    out$duality_gap_tolerance
if (!ok) quit(status = 1L)
