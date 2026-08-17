#!/usr/bin/env Rscript
# Canonical baseline generator. See inst/bench/PROTOCOL.md.

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1L) as.integer(args[[1]]) else 3L
seed <- if (length(args) >= 2L) as.integer(args[[2]]) else 20260816L
out_dir <- if (length(args) >= 3L) args[[3]] else "inst/bench/results/current"
threads <- if (length(args) >= 4L) as.integer(args[[4]]) else 1L
warmup <- 1L

if (!file.exists("DESCRIPTION")) {
  stop("Run inst/bench/run_protocol.R from the package root.", call. = FALSE)
}

suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
  library(rfugw)
  library(jsonlite)
})

source("inst/bench/protocol.R")

profile <- if (nzchar(Sys.getenv("RFUGW_FAST_FLAGS"))) "fast" else "conservative"
threads <- bench_pin_threads(threads)
meta <- bench_capture_env(seed, threads, warmup, reps, profile)

bench_archive_current(dirname(out_dir))

sizes_linear <- c(16L, 24L)
sizes_fgw <- c(12L, 16L)

runs <- list()
quality <- list()
k <- 1L

add_row <- function(suite, method, n, timing, extra = list()) {
  runs[[k]] <<- data.frame(
    suite = suite,
    method = method,
    n = n,
    valid = timing$valid,
    reject_reason = timing$reject_reason,
    prepare_ms = timing$prepare_ms,
    solve_ms = timing$solve_ms,
    e2e_ms = timing$e2e_ms,
    mem_solve_bytes = timing$mem_solve_bytes,
    mem_e2e_bytes = timing$mem_e2e_bytes,
    stringsAsFactors = FALSE
  )
  quality[[k]] <<- data.frame(
    suite = suite,
    method = method,
    n = n,
    valid = timing$valid,
    reject_reason = timing$reject_reason,
    status = timing$result$status %||% NA_character_,
    residual = timing$result$residual %||% NA_real_,
    value = tryCatch(rfugw_value(timing$result), error = function(e) NA_real_),
    stringsAsFactors = FALSE
  )
  k <<- k + 1L
}

for (n in sizes_linear) {
  prob <- bench_make_problem("linear", n, seed)
  stop <- bench_read_thresholds("ot_sinkhorn")$stop
  add_row("linear", "ot_sinkhorn", n, bench_run_split(
    prepare_fn = function() bench_make_problem("linear", n, seed),
    solve_fn = function(d) {
      ot_sinkhorn(d$M, d$p, d$q, epsilon = stop$epsilon,
                  max_iter = as.integer(stop$max_iter), tol = stop$tol)
    },
    warmup = warmup, reps = reps, method = "ot_sinkhorn",
    problem_for_quality = prob
  ))
  stop_e <- bench_read_thresholds("ot_emd")$stop
  add_row("linear", "ot_emd", n, bench_run_split(
    prepare_fn = function() bench_make_problem("linear", n, seed),
    solve_fn = function(d) {
      ot_emd(d$M, d$p, d$q, max_iter = as.integer(stop_e$max_iter), tol = stop_e$tol)
    },
    warmup = warmup, reps = reps, method = "ot_emd",
    problem_for_quality = prob
  ))
}

for (n in sizes_fgw) {
  prob <- bench_make_problem("fgw", n, seed)
  stop <- bench_read_thresholds("fgw_entropic")$stop
  add_row("fgw", "fgw_entropic", n, bench_run_split(
    prepare_fn = function() bench_make_problem("fgw", n, seed),
    solve_fn = function(d) {
      fgw_entropic(d$M, d$C1, d$C2, p = d$p, q = d$q,
                   alpha = stop$alpha, epsilon = stop$epsilon,
                   max_iter = as.integer(stop$max_iter), tol = stop$tol)
    },
    warmup = warmup, reps = reps, method = "fgw_entropic",
    problem_for_quality = prob
  ))
}

sizes_sr <- c(12L, 16L)
for (n in sizes_sr) {
  prob <- bench_make_problem("fgw", n, seed)
  stop_sr <- bench_read_thresholds("semirelaxed_gromov_wasserstein")$stop
  add_row("semirelaxed", "semirelaxed_gromov_wasserstein", n, bench_run_split(
    prepare_fn = function() bench_make_problem("fgw", n, seed),
    solve_fn = function(d) {
      semirelaxed_gromov_wasserstein(
        d$C1, d$C2, p = d$p,
        max_iter = as.integer(stop_sr$max_iter),
        tol_rel = stop_sr$tol_rel, tol_abs = stop_sr$tol_abs
      )
    },
    warmup = warmup, reps = reps, method = "semirelaxed_gromov_wasserstein",
    problem_for_quality = prob
  ))
  stop_srf <- bench_read_thresholds("semirelaxed_fused_gromov_wasserstein")$stop
  add_row("semirelaxed", "semirelaxed_fused_gromov_wasserstein", n, bench_run_split(
    prepare_fn = function() bench_make_problem("fgw", n, seed),
    solve_fn = function(d) {
      semirelaxed_fused_gromov_wasserstein(
        d$M, d$C1, d$C2, p = d$p, alpha = stop_srf$alpha,
        max_iter = as.integer(stop_srf$max_iter),
        tol_rel = stop_srf$tol_rel, tol_abs = stop_srf$tol_abs
      )
    },
    warmup = warmup, reps = reps, method = "semirelaxed_fused_gromov_wasserstein",
    problem_for_quality = prob
  ))
  stop_esr <- bench_read_thresholds("entropic_semirelaxed_gromov_wasserstein")$stop
  add_row("semirelaxed", "entropic_semirelaxed_gromov_wasserstein", n, bench_run_split(
    prepare_fn = function() bench_make_problem("fgw", n, seed),
    solve_fn = function(d) {
      entropic_semirelaxed_gromov_wasserstein(
        d$C1, d$C2, p = d$p, epsilon = stop_esr$epsilon,
        max_iter = as.integer(stop_esr$max_iter), tol = stop_esr$tol,
        check_every = as.integer(stop_esr$check_every %||% 10L)
      )
    },
    warmup = warmup, reps = reps, method = "entropic_semirelaxed_gromov_wasserstein",
    problem_for_quality = prob
  ))
}

sizes_partial <- c(8L, 10L)
for (n in sizes_partial) {
  prob <- bench_make_problem("fgw", n, seed)
  stop_pgw <- bench_read_thresholds("partial_gromov_wasserstein")$stop
  add_row("partial", "partial_gromov_wasserstein", n, bench_run_split(
    prepare_fn = function() bench_make_problem("fgw", n, seed),
    solve_fn = function(d) {
      partial_gromov_wasserstein(
        d$C1, d$C2, p = d$p, q = d$q, m = stop_pgw$m,
        numItermax = as.integer(stop_pgw$max_iter), tol = stop_pgw$tol,
        log = TRUE
      )
    },
    warmup = warmup, reps = reps, method = "partial_gromov_wasserstein",
    problem_for_quality = prob
  ))
  stop_pfgw <- bench_read_thresholds("partial_fused_gromov_wasserstein")$stop
  add_row("partial", "partial_fused_gromov_wasserstein", n, bench_run_split(
    prepare_fn = function() bench_make_problem("fgw", n, seed),
    solve_fn = function(d) {
      partial_fused_gromov_wasserstein(
        d$M, d$C1, d$C2, p = d$p, q = d$q, m = stop_pfgw$m,
        alpha = stop_pfgw$alpha,
        numItermax = as.integer(stop_pfgw$max_iter), tol = stop_pfgw$tol,
        log = TRUE
      )
    },
    warmup = warmup, reps = reps, method = "partial_fused_gromov_wasserstein",
    problem_for_quality = prob
  ))
  stop_ep <- bench_read_thresholds("entropic_partial_gromov_wasserstein")$stop
  add_row("partial", "entropic_partial_gromov_wasserstein", n, bench_run_split(
    prepare_fn = function() bench_make_problem("fgw", n, seed),
    solve_fn = function(d) {
      entropic_partial_gromov_wasserstein(
        d$C1, d$C2, p = d$p, q = d$q, m = stop_ep$m, reg = stop_ep$reg,
        numItermax = as.integer(stop_ep$max_iter), tol = stop_ep$tol,
        log = TRUE
      )
    },
    warmup = warmup, reps = reps, method = "entropic_partial_gromov_wasserstein",
    problem_for_quality = prob
  ))
}

sizes_ucoot <- c(8L, 12L)
for (n in sizes_ucoot) {
  prob <- bench_make_problem("ucoot", n, seed)
  stop_u <- bench_read_thresholds("unbalanced_co_optimal_transport")$stop
  add_row("ucoot", "unbalanced_co_optimal_transport", n, bench_run_split(
    prepare_fn = function() bench_make_problem("ucoot", n, seed),
    solve_fn = function(d) {
      unbalanced_co_optimal_transport(
        d$X, d$Y,
        wx_samp = d$wx_samp, wy_samp = d$wy_samp,
        wx_feat = d$wx_feat, wy_feat = d$wy_feat,
        reg_marginals = unlist(stop_u$reg_marginals),
        epsilon = unlist(stop_u$epsilon),
        max_iter = as.integer(stop_u$max_iter), tol = stop_u$tol,
        max_iter_ot = as.integer(stop_u$max_iter_ot), tol_ot = stop_u$tol_ot,
        log = TRUE
      )
    },
    warmup = warmup, reps = reps, method = "unbalanced_co_optimal_transport",
    problem_for_quality = prob
  ))
  stop_j <- bench_read_thresholds("fused_unbalanced_across_spaces_divergence")$stop
  add_row("ucoot", "fused_unbalanced_across_spaces_divergence", n, bench_run_split(
    prepare_fn = function() bench_make_problem("ucoot", n, seed),
    solve_fn = function(d) {
      fused_unbalanced_across_spaces_divergence(
        d$X, d$Y,
        wx_samp = d$wx_samp, wy_samp = d$wy_samp,
        wx_feat = d$wx_feat, wy_feat = d$wy_feat,
        reg_marginals = unlist(stop_j$reg_marginals),
        epsilon = unlist(stop_j$epsilon),
        reg_type = "joint",
        max_iter = as.integer(stop_j$max_iter), tol = stop_j$tol,
        max_iter_ot = as.integer(stop_j$max_iter_ot), tol_ot = stop_j$tol_ot,
        log = TRUE
      )
    },
    warmup = warmup, reps = reps, method = "fused_unbalanced_across_spaces_divergence",
    problem_for_quality = prob
  ))
}

run_df <- do.call(rbind, runs)
qual_df <- do.call(rbind, quality)
bench_write_baseline(meta, run_df, qual_df, out_dir)

cat("Wrote", out_dir, "\n")
print(run_df)
if (!all(run_df$valid)) {
  stop("Invalid-quality rows cannot count as baselines:\n",
       paste(run_df$reject_reason[!run_df$valid], collapse = "\n"),
       call. = FALSE)
}
if (identical(profile, "fast")) {
  message("profile=fast: recorded, but not a certified 0.1 baseline.")
}
