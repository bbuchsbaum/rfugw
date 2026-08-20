#!/usr/bin/env Rscript
# Canonical baseline generator. See inst/bench/PROTOCOL.md.

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1L) as.integer(args[[1]]) else 3L
seed <- if (length(args) >= 2L) as.integer(args[[2]]) else 20260816L
out_dir <- if (length(args) >= 3L) args[[3]] else "inst/bench/results/current"
threads <- if (length(args) >= 4L) as.integer(args[[4]]) else 1L
scale <- Sys.getenv("RFUGW_PROTOCOL_SCALE", unset = "full")
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
bench_validate_threshold_history()
suites <- bench_parse_suites(if (length(args) >= 5L) args[[5]] else "")

profile <- if (nzchar(Sys.getenv("RFUGW_FAST_FLAGS"))) "fast" else "conservative"
threads <- bench_pin_threads(threads)
meta <- bench_capture_env(seed, threads, warmup, reps, profile)
meta$suites <- paste(suites, collapse = ",")
meta$scale <- scale
meta$thresholds_md5 <- unname(tools::md5sum(bench_thresholds_path())[[1]])
meta$threshold_policy_version <- bench_read_thresholds()$`_policy`$schema_version

bench_archive_current(dirname(out_dir))

if (identical(scale, "pr")) {
  sizes_linear <- 16L
  sizes_fgw <- 12L
  sizes_fugw <- 8L
  sizes_sr <- 12L
  sizes_partial <- 8L
  sizes_ucoot <- 8L
  sizes_sampled <- 12L
} else if (identical(scale, "nightly")) {
  sizes_linear <- c(16L, 24L)
  sizes_fgw <- c(16L, 24L)
  sizes_fugw <- c(8L, 12L)
  sizes_sr <- c(12L, 16L)
  sizes_partial <- c(8L, 10L)
  sizes_ucoot <- c(8L, 12L)
  sizes_sampled <- 12L
} else {
  sizes_linear <- c(16L, 24L)
  sizes_fgw <- c(12L, 16L)
  sizes_fugw <- c(8L, 12L)
  sizes_sr <- c(12L, 16L)
  sizes_partial <- c(8L, 10L)
  sizes_ucoot <- c(8L, 12L)
  sizes_sampled <- 12L
}

runs <- list()
quality <- list()
k <- 1L

add_row <- function(suite, method, n, timing, extra = list()) {
  detail <- modifyList(list(
    budget_p = NA_integer_, budget_q = NA_integer_,
    reference_value = NA_real_, approximation_error_abs = NA_real_,
    approximation_error_rel = NA_real_
  ), extra)
  result <- timing$result
  runs[[k]] <<- data.frame(
    suite = suite,
    method = method,
    n = n,
    valid = timing$valid,
    certified = timing$certified,
    comparison_eligible = timing$comparison_eligible,
    performance_regression_eligible = timing$performance_regression_eligible,
    evidence_class = timing$evidence_class,
    reject_reason = timing$reject_reason,
    prepare_ms = timing$prepare_ms,
    setup_ms = timing$prepare_ms,
    solve_ms = timing$solve_ms,
    e2e_ms = timing$e2e_ms,
    mem_solve_bytes = timing$mem_solve_bytes,
    mem_e2e_bytes = timing$mem_e2e_bytes,
    requested_precision = result$requested_precision %||% NA_character_,
    effective_precision = result$effective_precision %||% NA_character_,
    compute_precision = result$compute_precision %||% NA_character_,
    requested_threads = result$requested_threads %||% threads,
    used_threads = result$used_threads %||% 1L,
    budget_p = detail$budget_p,
    budget_q = detail$budget_q,
    stringsAsFactors = FALSE
  )
  quality[[k]] <<- data.frame(
    suite = suite,
    method = method,
    n = n,
    valid = timing$valid,
    certified = timing$certified,
    comparison_eligible = timing$comparison_eligible,
    performance_regression_eligible = timing$performance_regression_eligible,
    evidence_class = timing$evidence_class,
    reject_reason = timing$reject_reason,
    status = result$status %||% NA_character_,
    converged = result$converged %||% FALSE,
    feasible = result$feasible %||% FALSE,
    objective_consistent = result$objective_consistent %||% FALSE,
    objective_components_consistent = result$objective_components_consistent %||% FALSE,
    inner_converged = result$inner_converged %||% NA,
    residual = result$residual %||% NA_real_,
    feasibility_residual = result$feasibility_residual %||% NA_real_,
    objective_residual = result$objective_residual %||% NA_real_,
    row_residual = result$row_residual %||% NA_real_,
    col_residual = result$col_residual %||% NA_real_,
    mass_residual = result$mass_residual %||% NA_real_,
    value = tryCatch(
      rfugw_value(result),
      error = function(e) result$gw_dist_estimated %||% NA_real_
    ),
    reference_value = detail$reference_value,
    approximation_error_abs = detail$approximation_error_abs,
    approximation_error_rel = detail$approximation_error_rel,
    budget_p = detail$budget_p,
    budget_q = detail$budget_q,
    stringsAsFactors = FALSE
  )
  k <<- k + 1L
}

if ("linear" %in% suites) for (n in sizes_linear) {
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

if ("fgw" %in% suites) for (n in sizes_fgw) {
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

if ("fugw" %in% suites) for (n in sizes_fugw) {
  prob <- bench_make_problem("fgw", n, seed)
  stop_f <- bench_read_thresholds("fugw_kl")$performance_stop
  add_row("fugw", "fugw_kl", n, bench_run_split(
    prepare_fn = function() bench_make_problem("fgw", n, seed),
    solve_fn = function(d) {
      fugw_kl(
        Cx = d$C1, Cy = d$C2, wx = d$p, wy = d$q, M = d$M,
        alpha = stop_f$alpha, epsilon = stop_f$epsilon,
        reg_marginals = unlist(stop_f$reg_marginals),
        max_iter = as.integer(stop_f$max_iter), tol = stop_f$tol,
        max_iter_ot = as.integer(stop_f$max_iter_ot), tol_ot = stop_f$tol_ot
      )
    },
    warmup = warmup, reps = reps, method = "fugw_kl",
    problem_for_quality = prob,
    evidence_class = "fixed_budget_performance"
  ))
}

if ("semirelaxed" %in% suites) for (n in sizes_sr) {
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

if ("partial" %in% suites) for (n in sizes_partial) {
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

if ("ucoot" %in% suites) for (n in sizes_ucoot) {
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

if ("sampled" %in% suites) for (n in sizes_sampled) {
  prob <- bench_make_problem("fgw", n, seed)
  stop_s <- bench_read_thresholds("sampled_gromov_wasserstein")$stop
  exact <- gromov_wasserstein(
    prob$C1, prob$C2, p = prob$p, q = prob$q,
    max_iter = 100L, tol_rel = 1e-8, tol_abs = 1e-8
  )
  reference_value <- rfugw_value(exact)
  budgets <- list(c(2L, 1L), c(4L, 2L), c(8L, 4L))
  for (budget in budgets) {
    timing <- bench_run_split(
      prepare_fn = function() bench_make_problem("fgw", n, seed),
      solve_fn = function(d) {
        sampled_gromov_wasserstein(
          d$C1, d$C2, p = d$p, q = d$q,
          nb_samples_grad = budget,
          epsilon = stop_s$epsilon,
          max_iter = as.integer(stop_s$max_iter),
          random_state = as.integer(seed),
          log = TRUE
        )
      },
      warmup = warmup, reps = reps, method = "sampled_gromov_wasserstein",
      problem_for_quality = prob,
      evidence_class = "fixed_budget_performance"
    )
    approx_value <- if (is.matrix(timing$result$plan)) {
      ot_gw_square(prob$C1, prob$C2, timing$result)
    } else {
      NA_real_
    }
    abs_error <- abs(approx_value - reference_value)
    add_row(
      "sampled", "sampled_gromov_wasserstein", n, timing,
      extra = list(
        budget_p = budget[[1]], budget_q = budget[[2]],
        reference_value = reference_value,
        approximation_error_abs = abs_error,
        approximation_error_rel = abs_error / max(abs(reference_value), .Machine$double.eps)
      )
    )
  }
}

run_df <- do.call(rbind, runs)
qual_df <- do.call(rbind, quality)
bench_write_baseline(meta, run_df, qual_df, out_dir)

cat("Wrote", out_dir, "\n")
print(run_df)
bad_cert <- run_df$evidence_class == "certified_comparison" & !run_df$comparison_eligible
bad_perf <- run_df$evidence_class == "fixed_budget_performance" &
  !run_df$performance_regression_eligible
if (any(bad_cert) || any(bad_perf)) {
  stop("Ineligible benchmark rows cannot count in their evidence class:\n",
       paste(run_df$reject_reason[bad_cert | bad_perf], collapse = "\n"),
       call. = FALSE)
}
if (identical(profile, "fast")) {
  message("profile=fast: recorded, but not a certified 0.1 baseline.")
}
