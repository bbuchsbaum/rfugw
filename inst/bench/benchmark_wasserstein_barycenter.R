#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rfugw))

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1L]] else tempfile(fileext = ".csv")

profile_call <- function(code) {
  profile <- tempfile(fileext = ".Rprofmem")
  gc()
  Rprofmem(profile)
  elapsed <- system.time(value <- force(code))[["elapsed"]]
  Rprofmem(NULL)
  lines <- readLines(profile, warn = FALSE)
  bytes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  list(
    value = value,
    elapsed_seconds = elapsed,
    allocated_bytes = sum(bytes[is.finite(bytes)])
  )
}

support <- seq(0, 1, length.out = 5L)
support_cost <- outer(support, support, function(x, y) (x - y)^2)
costs <- list(support_cost, support_cost, support_cost)
measures <- list(
  c(0.1, 0.2, 0.3, 0.2, 0.2),
  c(0.2, 0.1, 0.2, 0.3, 0.2),
  c(0.3, 0.2, 0.1, 0.2, 0.2)
)
coefficients <- c(0.2, 0.3, 0.5)

exact_profile <- profile_call(ot_barycenter_weights(
  costs, measures, coefficients, mode = "exact", tol = 1e-10
))
regularized_profile <- profile_call(ot_barycenter_weights(
  costs, measures, coefficients,
  mode = "regularized",
  support_cost = support_cost,
  epsilon = 0.1,
  max_iter = 300L,
  tol = 2e-5,
  sinkhorn_method = "auto",
  sinkhorn_max_iter = 5000L,
  sinkhorn_tol = 1e-8
))
warm_profile <- profile_call(ot_barycenter_weights(
  costs, measures, coefficients,
  mode = "regularized",
  support_cost = support_cost,
  epsilon = 0.1,
  init_state = regularized_profile$value,
  max_iter = 300L,
  tol = 2e-5,
  sinkhorn_method = "auto",
  sinkhorn_max_iter = 5000L,
  sinkhorn_tol = 1e-8
))

profiles <- list(exact_profile, regularized_profile, warm_profile)
labels <- c("exact_joint_lp", "regularized_cold", "regularized_warm")
result <- do.call(rbind, Map(function(profile, label) {
  fit <- profile$value
  data.frame(
    mode = label,
    status = fit$status,
    converged = fit$converged,
    objective = fit$barycenter_objective,
    simplex_residual = fit$simplex_residual,
    kkt_residual = fit$kkt_residual,
    kkt_tolerance = fit$kkt_tolerance,
    outer_iterations = fit$iterations,
    component_solves = fit$component_solves,
    component_iterations = fit$total_component_iterations,
    warm_reuse_count = fit$warm_start_reuse_count %||% 0L,
    elapsed_seconds = profile$elapsed_seconds,
    allocated_bytes = profile$allocated_bytes,
    stringsAsFactors = FALSE
  )
}, profiles, labels))

if (!all(result$converged) || any(result$kkt_residual > result$kkt_tolerance)) {
  stop("Barycenter benchmark quality gate failed.", call. = FALSE)
}
if (warm_profile$value$total_component_iterations >
    regularized_profile$value$total_component_iterations ||
    max(abs(warm_profile$value$weights -
      regularized_profile$value$weights)) > 1e-12) {
  stop("Barycenter warm-state work/parity gate failed.", call. = FALSE)
}
utils::write.csv(result, output, row.names = FALSE)
print(result)
cat(sprintf("wrote %s\n", normalizePath(output, mustWork = FALSE)))
