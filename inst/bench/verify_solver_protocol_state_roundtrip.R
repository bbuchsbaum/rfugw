#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rfugw))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--consume")) {
  bundle <- readRDS(args[[2L]])
  fit <- transport_solve(bundle$problem, bundle$state)
  stopifnot(
    isTRUE(fit$converged),
    isTRUE(fit$warm_start_accepted),
    identical(rfugw_status(fit), "converged"),
    max(abs(rfugw_plan(fit) - bundle$plan)) <= 2e-10
  )
  cat("fresh installed-session protocol state: PASS\n")
  quit(status = 0L)
}

set.seed(1204)
source_points <- matrix(rnorm(42), 14L, 3L)
target_points <- matrix(rnorm(30), 10L, 3L)
cost <- outer(rowSums(source_points^2), rowSums(target_points^2), "+") -
  2 * tcrossprod(source_points, target_points)
problem <- transport_problem_sinkhorn(
  cost,
  epsilon = 0.7,
  method = "auto",
  max_iter = 5000L,
  tol = 1e-10,
  mass_policy = "probability"
)
fit <- transport_solve(problem)
stopifnot(isTRUE(fit$converged))
bundle_path <- tempfile(fileext = ".rds")
saveRDS(
  list(problem = problem, state = rfugw_state(fit), plan = rfugw_plan(fit)),
  bundle_path
)
script <- system.file(
  "bench", "verify_solver_protocol_state_roundtrip.R", package = "rfugw"
)
if (!nzchar(script)) {
  script <- normalizePath(
    "inst/bench/verify_solver_protocol_state_roundtrip.R", mustWork = TRUE
  )
}
status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c("--vanilla", shQuote(script), "--consume", shQuote(bundle_path))
)
if (!identical(status, 0L)) {
  stop("Fresh-session state consumer failed.", call. = FALSE)
}
cat("solver protocol state round-trip: PASS\n")
