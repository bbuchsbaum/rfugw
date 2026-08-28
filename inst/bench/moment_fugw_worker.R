#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: moment_fugw_worker.R SPEC.json RESULT.json", call. = FALSE)
}
`%||%` <- function(x, y) if (is.null(x)) y else x
spec_path <- normalizePath(args[[1L]], mustWork = TRUE)
result_path <- args[[2L]]

write_result <- function(value) {
  dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    value, result_path, auto_unbox = TRUE, pretty = TRUE, null = "null",
    na = "null"
  )
}

tryCatch({
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The worker requires jsonlite.", call. = FALSE)
  }
  spec <- jsonlite::fromJSON(spec_path, simplifyVector = TRUE)
  package_mode <- as.character(spec$package_mode %||% "installed")[[1L]]
  if (identical(package_mode, "source")) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("Source mode requires pkgload.", call. = FALSE)
    }
    package_root <- normalizePath(spec$package_root, mustWork = TRUE)
    pkgload::load_all(package_root, quiet = TRUE, export_all = FALSE)
  } else {
    suppressPackageStartupMessages(library(rfugw))
  }

  protocol_path <- normalizePath(spec$protocol_path, mustWork = TRUE)
  benchmark_path <- normalizePath(spec$benchmark_path, mustWork = TRUE)
  source(protocol_path, local = .GlobalEnv)
  source(benchmark_path, local = .GlobalEnv)
  bench_pin_threads(as.integer(spec$threads %||% 1L))

  action <- as.character(spec$action)[[1L]]
  if (identical(action, "baseline")) {
    write_result(list(
      valid = TRUE,
      action = action,
      package_mode = package_mode,
      elapsed_seconds = 0,
      r_visible_allocation_bytes = NA_real_
    ))
    quit(save = "no", status = 0L)
  }

  if (action %in% c("solve", "allocation")) {
    setup_start <- proc.time()[["elapsed"]]
    problem <- bench_make_moment_fugw_problem(
      n = as.integer(spec$n),
      structure_rank = as.integer(spec$structure_rank),
      feature_dimension = as.integer(spec$feature_dimension),
      seed = as.integer(spec$seed),
      guard_dense = TRUE
    )
    setup_seconds <- proc.time()[["elapsed"]] - setup_start
    solve_call <- function() bench_solve_moment_fugw(
        problem,
        max_iter = as.integer(spec$max_iter),
        max_iter_ot = as.integer(spec$max_iter_ot),
        tol = as.numeric(spec$tol),
        tol_ot = as.numeric(spec$tol_ot),
        epsilon = as.numeric(spec$epsilon),
        reg_marginals = as.numeric(spec$reg_marginals),
        block_size = as.integer(spec$block_size)
      )
    solve_start <- proc.time()[["elapsed"]]
    allocation <- if (identical(action, "allocation")) {
      bench_profile_r_visible_allocation(solve_call)
    } else {
      list(value = solve_call(), bytes = NA_real_, events = NA_integer_)
    }
    fit <- allocation$value
    solve_seconds <- proc.time()[["elapsed"]] - solve_start
    evidence_class <- as.character(spec$evidence_class)[[1L]]
    quality <- bench_check_quality(
      "fugw_factorized", fit,
      list(expected_work = list(
        outer_iterations = as.integer(spec$max_iter),
        inner_iterations = as.integer(spec$max_iter_ot)
      )),
      evidence_class = evidence_class
    )
    sentinel <- bench_moment_fugw_dense_sentinel(fit, problem)
    valid <- isTRUE(quality$valid) && isTRUE(sentinel$passed)
    reject <- c(
      quality$reject_reason %||% "",
      if (!isTRUE(sentinel$passed)) sentinel$reasons else character()
    )
    reject <- reject[nzchar(reject)]
    fit_path <- as.character(spec$fit_path %||% "")[[1L]]
    if (identical(action, "solve") && nzchar(fit_path)) {
      dir.create(dirname(fit_path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(list(
        fit = fit,
        problem = problem,
        source_quality = quality,
        source_spec = spec
      ), fit_path, compress = FALSE)
    }
    write_result(list(
      valid = valid,
      action = "solve",
      method = "fugw_factorized",
      evidence_class = evidence_class,
      n_source = problem$n_source,
      n_target = problem$n_target,
      structure_rank = problem$structure_rank,
      feature_dimension = problem$feature_dimension,
      status = fit$status,
      stationarity_classification =
        fit$certificate$outer_stationarity$classification,
      certified = isTRUE(quality$certified),
      comparison_eligible = isTRUE(quality$comparison_eligible),
      timing_eligible = isTRUE(quality$timing_eligible),
      scaling_eligible = isTRUE(quality$scaling_eligible),
      memory_contract_eligible = identical(action, "allocation") && valid,
      threshold_update_eligible = isTRUE(quality$threshold_update_eligible),
      contract_status = quality$contract_status,
      reject_reason = paste(unique(reject), collapse = ";"),
      limitation_reason = quality$limitation_reason,
      objective = fit$fugw_cost,
      iterations = fit$iterations,
      inner_iterations = fit$inner_iterations,
      residual = fit$residual,
      inner_residual = fit$inner_residual,
      setup_seconds = setup_seconds,
      operation_seconds = solve_seconds,
      end_to_end_seconds = setup_seconds + solve_seconds,
      result_object_bytes = as.numeric(object.size(fit)),
      r_visible_allocation_bytes = allocation$bytes,
      r_visible_allocation_events = allocation$events,
      dense_sentinel = sentinel,
      fit_path = if (nzchar(fit_path)) fit_path else NULL
    ))
    quit(save = "no", status = 0L)
  }

  if (!action %in% c("objective", "apply", "adjoint")) {
    stop("Unknown worker action: ", action, call. = FALSE)
  }
  fit_path <- normalizePath(spec$fit_path, mustWork = TRUE)
  bundle <- readRDS(fit_path)
  fit <- bundle$fit
  problem <- bundle$problem
  source_evidence_class <- as.character(
    spec$source_evidence_class %||% "operator_throughput"
  )[[1L]]
  quality <- bench_check_quality(
    "fugw_factorized", fit,
    if (identical(source_evidence_class, "fixed_work_scaling")) {
      list(expected_work = list(
        outer_iterations = as.integer(bundle$source_spec$max_iter),
        inner_iterations = as.integer(bundle$source_spec$max_iter_ot)
      ))
    } else list(),
    evidence_class = source_evidence_class
  )
  batch <- as.integer(spec$batch %||% 1L)
  run_operation <- function() {
    if (identical(action, "objective")) {
      value <- bench_recompute_moment_fugw_objective(bundle)
      valid <- all(is.finite(unlist(value))) &&
        isTRUE(all.equal(
          as.numeric(value$fugw_cost), as.numeric(fit$fugw_cost),
          tolerance = max(1e-10, 100 * fit$objective_tolerance)
        ))
      return(list(
        value = value, valid = valid, rows = 1L,
        columns = length(value), checksum = sum(unlist(value))
      ))
    }
    if (identical(action, "apply")) {
      input <- bench_moment_fugw_map_batch(
        problem$n_target, batch, as.integer(spec$seed)
      )
      value <- transport_plan_apply(fit$plans$sample, input)
      valid <- is.matrix(value) && all(is.finite(value)) &&
        identical(dim(value), c(problem$n_source, batch))
    } else {
      input <- bench_moment_fugw_map_batch(
        problem$n_source, batch, as.integer(spec$seed)
      )
      value <- transport_plan_adjoint(fit$plans$sample, input)
      valid <- is.matrix(value) && all(is.finite(value)) &&
        identical(dim(value), c(problem$n_target, batch))
    }
    list(
      value = value, valid = valid, rows = nrow(value),
      columns = ncol(value), checksum = sum(value)
    )
  }
  operation_start <- proc.time()[["elapsed"]]
  allocation <- if (isTRUE(spec$measure_r_visible %||% FALSE)) {
    bench_profile_r_visible_allocation(run_operation)
  } else {
    list(value = run_operation(), bytes = NA_real_, events = NA_integer_)
  }
  operation <- allocation$value
  operation_seconds <- proc.time()[["elapsed"]] - operation_start
  operation_valid <- operation$valid
  output_rows <- operation$rows
  output_columns <- operation$columns
  checksum <- operation$checksum
  sentinel <- bench_moment_fugw_dense_sentinel(fit, problem)
  source_eligible <- if (identical(
      source_evidence_class, "fixed_work_scaling")) {
    isTRUE(quality$scaling_eligible)
  } else {
    isTRUE(quality$operator_eligible)
  }
  valid <- source_eligible && operation_valid &&
    isTRUE(sentinel$passed)
  reject <- c(
    quality$reject_reason %||% "",
    if (!operation_valid) paste0(action, "_invalid_output") else character(),
    if (!isTRUE(sentinel$passed)) sentinel$reasons else character()
  )
  reject <- reject[nzchar(reject)]
  write_result(list(
    valid = valid,
    action = action,
    method = "fugw_factorized",
    evidence_class = if (identical(
        source_evidence_class, "fixed_work_scaling")) {
      "matrix_free_memory_contract"
    } else "operator_throughput",
    n_source = problem$n_source,
    n_target = problem$n_target,
    structure_rank = problem$structure_rank,
    feature_dimension = problem$feature_dimension,
    batch = batch,
    source_status = fit$status,
    source_certified = isTRUE(quality$certified),
    timing_eligible = isTRUE(quality$timing_eligible),
    operator_eligible = isTRUE(quality$operator_eligible),
    memory_contract_eligible = identical(
      source_evidence_class, "fixed_work_scaling"
    ) && isTRUE(quality$scaling_eligible),
    threshold_update_eligible = isTRUE(quality$threshold_update_eligible),
    contract_status = quality$contract_status,
    reject_reason = paste(unique(reject), collapse = ";"),
    operation_seconds = operation_seconds,
    output_rows = output_rows,
    output_columns = output_columns,
    checksum = checksum,
    r_visible_allocation_bytes = allocation$bytes,
    r_visible_allocation_events = allocation$events,
    dense_sentinel = sentinel
  ))
}, error = function(error) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    write_result(list(
      valid = FALSE,
      infrastructure_error = TRUE,
      error = conditionMessage(error)
    ))
  }
  message(conditionMessage(error))
  quit(save = "no", status = 1L)
})
