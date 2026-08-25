#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
  library(rfugw)
})

set.seed(42)
X <- matrix(rnorm(60), 20, 3)
Y <- matrix(rnorm(45), 15, 3)
Z <- rbind(X, Y)
M <- as.matrix(stats::dist(Z))[seq_len(nrow(X)), nrow(X) + seq_len(nrow(Y))]^2
p <- rep(1 / nrow(X), nrow(X))
q <- rep(1 / nrow(Y), nrow(Y))
epsilon <- c(1, 0.6, 0.35, 0.2)

state <- NULL
rows <- lapply(epsilon, function(reg) {
  warm <- ot_sinkhorn(
    M, p, q, epsilon = reg, method = "scaling",
    max_iter = 5000L, tol = 1e-9, init_duals = state
  )
  cold <- ot_sinkhorn(
    M, p, q, epsilon = reg, method = "scaling",
    max_iter = 5000L, tol = 1e-9
  )
  stopifnot(warm$converged, cold$converged)
  state <<- warm$dual_state
  data.frame(
    method = "scaling",
    epsilon = reg,
    cold_iterations = cold$iterations,
    warm_iterations = warm$iterations,
    max_plan_error = max(abs(warm$plan - cold$plan)),
    objective_error = abs(warm$ot_dist - cold$ot_dist)
  )
})

result <- do.call(rbind, rows)
stopifnot(
  sum(result$warm_iterations) <= sum(result$cold_iterations),
  max(result$max_plan_error) <= 2e-8,
  max(result$objective_error) <= 2e-8
)
print(result)
