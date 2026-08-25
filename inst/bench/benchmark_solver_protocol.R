#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rfugw))

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1L]] else tempfile(fileext = ".csv")

set.seed(20260824)
ns <- 40L
nt <- 32L
source_points <- matrix(rnorm(ns * 4L), ns, 4L)
target_points <- matrix(rnorm(nt * 4L), nt, 4L)
cost <- outer(rowSums(source_points^2), rowSums(target_points^2), "+") -
  2 * tcrossprod(source_points, target_points)
cost[cost < 0 & cost > -1e-12] <- 0
source <- seq_len(ns) / sum(seq_len(ns))
target <- rev(seq_len(nt)) / sum(seq_len(nt))
tolerances <- c(1e-5, 1e-7, 1e-9)

state <- NULL
rows <- vector("list", length(tolerances))
for (i in seq_along(tolerances)) {
  problem <- transport_problem_sinkhorn(
    cost,
    source,
    target,
    epsilon = 0.65,
    method = "log",
    max_iter = 10000L,
    tol = tolerances[[i]],
    mass_policy = "probability"
  )
  warm_elapsed <- system.time({
    warm <- transport_solve(problem, init_state = state)
  })[["elapsed"]]
  cold_elapsed <- system.time({
    cold <- transport_solve(problem)
  })[["elapsed"]]
  if (!isTRUE(warm$converged) || !isTRUE(cold$converged)) {
    stop("Protocol benchmark requires both paths to converge.", call. = FALSE)
  }
  rows[[i]] <- data.frame(
    tolerance = tolerances[[i]],
    warm_iterations = warm$iterations,
    cold_iterations = cold$iterations,
    warm_elapsed_seconds = warm_elapsed,
    cold_elapsed_seconds = cold_elapsed,
    max_plan_difference = max(abs(warm$plan - cold$plan)),
    regularized_objective_difference = abs(
      warm$regularized_objective - cold$regularized_objective
    ),
    warm_residual = warm$residual,
    cold_residual = cold$residual,
    warm_state_accepted = warm$warm_start_accepted,
    stringsAsFactors = FALSE
  )
  state <- rfugw_state(warm)
}

result <- do.call(rbind, rows)
if (sum(result$warm_iterations) > sum(result$cold_iterations)) {
  stop("Warm continuation used more total iterations than cold restarts.",
       call. = FALSE)
}
if (max(result$max_plan_difference) > 2e-8 ||
    max(result$regularized_objective_difference) > 2e-8) {
  stop("Warm and cold answers differ beyond the benchmark quality gate.",
       call. = FALSE)
}
utils::write.csv(result, output, row.names = FALSE)
print(result)
cat(sprintf(
  "total iterations: warm=%d cold=%d\n",
  sum(result$warm_iterations), sum(result$cold_iterations)
))
cat(sprintf("wrote %s\n", normalizePath(output, mustWork = FALSE)))
