#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) >= 1L) as.integer(args[[1L]]) else 1000L
repetitions <- if (length(args) >= 2L) as.integer(args[[2L]]) else 3L
label <- if (length(args) >= 3L) args[[3L]] else "profile"
output <- if (length(args) >= 4L) {
  args[[4L]]
} else file.path("inst", "bench", "results", paste0(
  "moment-fugw-kernel-profile-", label, "-n", n, ".csv"
))

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Source profiling requires pkgload.", call. = FALSE)
}
pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
source(file.path("inst", "bench", "protocol.R"), local = .GlobalEnv)
source(
  file.path("inst", "bench", "benchmark_moment_fugw.R"),
  local = .GlobalEnv
)
bench_pin_threads(1L)

internal <- function(name) getFromNamespace(name, "rfugw")
profile_call <- function(component, fn, reps = repetitions) {
  elapsed <- numeric(reps)
  for (index in seq_len(reps)) {
    gc(FALSE)
    started <- proc.time()[["elapsed"]]
    value <- fn()
    elapsed[[index]] <- proc.time()[["elapsed"]] - started
    if (is.null(value)) stop("Profiled call returned NULL.", call. = FALSE)
  }
  data.frame(
    label = label,
    n_source = n,
    n_target = n,
    component = component,
    repetitions = reps,
    median_seconds = stats::median(elapsed),
    minimum_seconds = min(elapsed),
    maximum_seconds = max(elapsed),
    stringsAsFactors = FALSE
  )
}

config <- bench_moment_fugw_profile("full")
problem <- bench_make_moment_fugw_problem(
  n,
  structure_rank = config$structure_rank,
  feature_dimension = config$feature_dimension,
  seed = config$seed
)
fit <- bench_solve_moment_fugw(
  problem,
  max_iter = 1L,
  max_iter_ot = 1L,
  tol = config$tol,
  tol_ot = config$tol_ot,
  epsilon = config$epsilon,
  reg_marginals = config$reg_marginals,
  block_size = config$block_size
)
quality <- bench_check_quality(
  "fugw_factorized", fit,
  list(expected_work = list(outer_iterations = 1L)),
  evidence_class = "fixed_work_scaling"
)
if (!isTRUE(quality$scaling_eligible)) {
  stop("Profile fixture failed fixed-work quality: ", quality$reject_reason,
       call. = FALSE)
}

sample_plan <- fit$implicit_plans$sample
feature_plan <- fit$implicit_plans$feature
sample_state <- internal(".factorized_plan_state")(sample_plan)
feature_state <- internal(".factorized_plan_state")(feature_plan)
source_factors <- internal(".cost_factor_pair")(problem$Cx)
target_factors <- internal(".cost_factor_pair")(problem$Cy)
feature_terms <- internal(".cost_native_terms")(problem$M)

native_ti <- function() {
  terms <- sample_state$cost
  native <- internal("cpp_ot_sinkhorn_unbalanced_ti_factorized")(
    terms$row, terms$column, terms$left, terms$right,
    sample_state$source_measure, sample_state$target_measure,
    sample_state$epsilon,
    config$reg_marginals[[1L]], config$reg_marginals[[2L]],
    1L, 1e-30, config$block_size,
    sample_state$source_bar, sample_state$target_bar
  )
  native$source_bar
}

rows <- list(
  profile_call("whole_fixed_work_solver", function() {
    bench_solve_moment_fugw(
      problem, 1L, 1L, config$tol, config$tol_ot, config$epsilon,
      config$reg_marginals, config$block_size
    )
  }),
  profile_call("ti_native_sweeps", native_ti),
  profile_call("plan_statistics", function() {
    internal(".factorized_plan_stats")(sample_state)
  }),
  profile_call("moment_reduction", function() {
    internal(".factorized_structure_moments")(
      sample_plan, source_factors$right, target_factors$right
    )
  }),
  profile_call("plan_difference", function() {
    internal(".factorized_plan_difference")(
      sample_plan, feature_plan, config$block_size
    )
  }),
  profile_call("objective_recomputation", function() {
    internal(".fugw_factorized_objective")(
      sample_plan, feature_plan, problem$M, fit$moments$feature,
      source_factors$left, target_factors$left,
      problem$wx, problem$wy, fit$reg_marginals, fit$regularization,
      fit$feature_weight, fit$structure_weight
    )
  }),
  profile_call("final_two_sided_kkt", function() {
    internal(".fugw_final_block_audits")(
      sample_plan, feature_plan, fit$moments$sample, fit$moments$feature,
      source_factors$left, target_factors$left, feature_terms,
      problem$wx, problem$wy, fit$reg_marginals, fit$regularization,
      fit$feature_weight, fit$structure_weight, config$tol_ot,
      config$block_size
    )
  }),
  profile_call("apply_batch_10", function() {
    transport_plan_apply(
      sample_plan, bench_moment_fugw_map_batch(n, 10L, config$seed)
    )
  }),
  profile_call("adjoint_batch_10", function() {
    transport_plan_adjoint(
      sample_plan, bench_moment_fugw_map_batch(n, 10L, config$seed)
    )
  })
)
rows <- do.call(rbind, rows)
rows$share_of_whole_solver <- rows$median_seconds /
  rows$median_seconds[rows$component == "whole_fixed_work_solver"]
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(rows, output, row.names = FALSE)
print(rows)
cat("Profile artifact:", normalizePath(output), "\n")
