#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_csv <- if (length(args)) args[[1]] else "inst/bench/linear-foundations-baseline.csv"
reps <- if (length(args) >= 2L) as.integer(args[[2]]) else 3L
seed <- if (length(args) >= 3L) as.integer(args[[3]]) else 20260820L

suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
  library(rfugw)
})
source("inst/bench/protocol.R")
bench_pin_threads(1L)

certify <- function(out) {
  ok <- identical(out$status, "converged") &&
    isTRUE(out$converged) && isTRUE(out$feasible) &&
    isTRUE(out$objective_consistent) &&
    isTRUE(out$objective_components_consistent)
  if (!is.null(out$inner_converged) && length(out$inner_converged) &&
      !is.na(out$inner_converged[[1]])) {
    ok <- ok && isTRUE(out$inner_converged)
  }
  if (!ok) stop("Foundation benchmark result lacks its formulation certificate.", call. = FALSE)
  invisible(TRUE)
}

measure <- function(name, n, prepare, solve) {
  warm <- prepare()
  warm_out <- solve(warm)
  certify(warm_out)
  setup <- solve_ms <- numeric(reps)
  last <- warm_out
  for (i in seq_len(reps)) {
    p <- bench_time_ms(prepare)
    s <- bench_time_ms(function() solve(p$value))
    certify(s$value)
    setup[[i]] <- p$ms
    solve_ms[[i]] <- s$ms
    last <- s$value
  }
  mem_solve <- mem_e2e <- NA_real_
  if (requireNamespace("bench", quietly = TRUE)) {
    prepared_for_memory <- prepare()
    mem_solve <- as.numeric(bench::mark(
      solve(prepared_for_memory), iterations = 1L, check = FALSE, memory = TRUE
    )$mem_alloc[[1]])
    mem_e2e <- as.numeric(bench::mark({
      d <- prepare(); solve(d)
    }, iterations = 1L, check = FALSE, memory = TRUE)$mem_alloc[[1]])
  }
  info <- Sys.info()
  data.frame(
    method = name,
    n = n,
    evidence_channel = "local",
    certified = TRUE,
    status = last$status,
    objective = rfugw_value(last),
    residual = last$residual,
    feasibility_residual = last$feasibility_residual,
    objective_residual = last$objective_residual,
    mass_residual = last$mass_residual %||% NA_real_,
    setup_ms = median(setup),
    solve_ms = median(solve_ms),
    e2e_ms = median(setup + solve_ms),
    mem_solve_bytes = mem_solve,
    mem_e2e_bytes = mem_e2e,
    requested_precision = last$requested_precision %||% NA_character_,
    effective_precision = last$effective_precision %||% NA_character_,
    requested_threads = last$requested_threads %||% 1L,
    used_threads = last$used_threads %||% 1L,
    commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE)[1],
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    sysname = unname(info[["sysname"]]),
    machine = unname(info[["machine"]]),
    omp_num_threads = Sys.getenv("OMP_NUM_THREADS", ""),
    blas_threads = Sys.getenv("VECLIB_MAXIMUM_THREADS", ""),
    seed = seed,
    reps = reps,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    stringsAsFactors = FALSE
  )
}

rows <- list()
for (n in c(16L, 32L)) {
  partial_prepare <- function() {
    d <- bench_make_problem("linear", n, seed)
    list(M = d$M, p = d$p, q = d$q)
  }
  rows[[length(rows) + 1L]] <- measure(
    "ot_partial_emd", n, partial_prepare,
    function(d) ot_partial_emd(d$M, d$p, d$q, mass = 0.7)
  )
  rows[[length(rows) + 1L]] <- measure(
    "ot_sinkhorn_unbalanced_log", n, partial_prepare,
    function(d) ot_sinkhorn_unbalanced(
      d$M, d$p, d$q, epsilon = 0.08, rho = 2, method = "log",
      max_iter = 2000L, tol = 1e-6
    )
  )
}

out <- do.call(rbind, rows)
utils::write.csv(out, out_csv, row.names = FALSE)
print(out)
cat("Wrote", out_csv, "\n")
