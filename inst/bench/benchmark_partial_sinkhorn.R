#!/usr/bin/env Rscript

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools to run this source-tree benchmark.", call. = FALSE)
}
devtools::load_all(quiet = TRUE)

profile_call <- function(call) {
  profile <- tempfile("rfugw-partial-sinkhorn-", fileext = ".mem")
  on.exit(unlink(profile), add = TRUE)
  gc()
  Rprofmem(profile)
  elapsed <- system.time(value <- force(call))[["elapsed"]]
  Rprofmem(NULL)
  allocations <- readLines(profile, warn = FALSE)
  bytes <- suppressWarnings(as.numeric(sub(" .*", "", allocations)))
  list(
    value = value,
    elapsed_ms = 1000 * elapsed,
    allocation_count = sum(is.finite(bytes)),
    allocation_bytes = sum(bytes[is.finite(bytes)])
  )
}

set.seed(20260824L)
M <- matrix(runif(120), 12L, 10L)
p <- runif(12); p <- p / sum(p)
q <- runif(10); q <- q / sum(q)
mass <- 0.75
epsilon <- 0.08

run_case <- function(label, method, cost = M, state = NULL) {
  measured <- profile_call(ot_partial_sinkhorn(
    cost, p, q, mass = mass, epsilon = epsilon,
    method = method, max_iter = 20000L, tol = 1e-9,
    check_every = 10L, init_state = state
  ))
  fit <- measured$value
  stopifnot(
    isTRUE(fit$converged), isTRUE(fit$feasible),
    isTRUE(fit$objective_consistent),
    abs(fit$duality_gap) <= fit$duality_gap_tolerance
  )
  data.frame(
    case = label,
    requested_method = method,
    effective_method = fit$effective_sinkhorn_method,
    warm_started = fit$warm_started,
    iterations = fit$iterations,
    residual = fit$residual,
    feasibility_residual = fit$feasibility_residual,
    objective_residual = fit$objective_residual,
    duality_gap = fit$duality_gap,
    dynamic_range = fit$sinkhorn_dynamic_range,
    scaling_threshold = fit$sinkhorn_scaling_threshold,
    elapsed_ms = measured$elapsed_ms,
    allocation_count = measured$allocation_count,
    allocation_bytes = measured$allocation_bytes,
    status = fit$status,
    certified = fit$converged,
    stringsAsFactors = FALSE
  )
}

moderate_scaling <- run_case("moderate_scaling", "scaling")
moderate_log <- run_case("moderate_log", "log")
cold_log <- ot_partial_sinkhorn(
  M, p, q, mass = mass, epsilon = epsilon,
  method = "log", max_iter = 20000L, tol = 1e-9, check_every = 10L
)
warm_log <- run_case("moderate_log_warm", "log", state = cold_log)
adversarial_cost <- M
adversarial_cost[1, 1] <- 100
adversarial_log <- run_case("adversarial_auto", "auto", adversarial_cost)

result <- rbind(moderate_scaling, moderate_log, warm_log, adversarial_log)
result$seed <- 20260824L
result$r_version <- as.character(getRversion())
result$sysname <- Sys.info()[["sysname"]]
result$machine <- Sys.info()[["machine"]]
result$commit <- tryCatch(
  system2("git", c("rev-parse", "--short=12", "HEAD"), stdout = TRUE),
  error = function(e) "unknown"
)

output <- Sys.getenv(
  "RFUGW_PARTIAL_SINKHORN_BENCH_OUTPUT",
  unset = file.path("inst", "bench", "partial-sinkhorn-baseline.csv")
)
utils::write.csv(result, output, row.names = FALSE)
print(result)
