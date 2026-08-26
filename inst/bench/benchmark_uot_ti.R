#!/usr/bin/env Rscript

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install `devtools` to run the source-tree benchmark.", call. = FALSE)
}
devtools::load_all(quiet = TRUE)

peak_rss_call <- function(expression, environment = parent.frame()) {
  expression <- substitute(expression)
  started <- proc.time()[["elapsed"]]
  rss_kib <- function(pid) {
    rss <- suppressWarnings(system2(
      "/bin/ps", c("-o", "rss=", "-p", as.integer(pid)),
      stdout = TRUE, stderr = FALSE
    ))
    rss <- suppressWarnings(as.numeric(trimws(rss)))
    if (length(rss) && is.finite(rss[[1L]])) rss[[1L]] else NA_real_
  }
  if (.Platform$OS.type == "unix" &&
      requireNamespace("parallel", quietly = TRUE)) {
    job <- parallel::mcparallel(
      eval(expression, envir = environment),
      mc.set.seed = FALSE,
      silent = TRUE
    )
    peak_kib <- rss_kib(Sys.getpid())
    if (!is.finite(peak_kib)) peak_kib <- 0
    collected <- NULL
    repeat {
      rss <- rss_kib(job$pid)
      if (is.finite(rss)) {
        peak_kib <- max(peak_kib, rss, na.rm = TRUE)
      }
      collected <- parallel::mccollect(job, wait = FALSE)
      if (!is.null(collected)) break
      Sys.sleep(.005)
    }
    value <- collected[[1L]]
    parent_rss <- rss_kib(Sys.getpid())
    if (is.finite(parent_rss)) peak_kib <- max(peak_kib, parent_rss)
    return(list(
      value = value,
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      peak_rss_bytes = if (peak_kib > 0) peak_kib * 1024 else NA_real_
    ))
  }
  value <- eval(expression, envir = environment)
  list(
    value = value,
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    peak_rss_bytes = NA_real_
  )
}

n <- 300L
k <- 8L
epsilon <- .3
rho <- c(2, 3)
p <- seq_len(n)
p <- p / sum(p)
q <- rev(seq_len(n))
q <- 1.4 * q / sum(q)
source <- rep(seq_len(n), each = k)
offset <- rep(0:(k - 1L), times = n)
target <- ((source - 1L + offset) %% n) + 1L
distance <- pmin(abs(source - target), n - abs(source - target)) / n
edges <- data.frame(source = source, target = target, cost = distance^2)

dense_cost <- matrix(25, n, n)
dense_cost[cbind(edges$source, edges$target)] <- edges$cost

dense_profile <- peak_rss_call(
  ot_sinkhorn_unbalanced_ti(
    dense_cost, p, q, epsilon, rho,
    max_iter = 5000L, tol = 1e-9, plan = "operator"
  )
)
sparse_profile <- peak_rss_call(
  ot_sinkhorn_unbalanced_ti(
    edges, p, q, epsilon, rho,
    n_source = n, n_target = n,
    max_iter = 5000L, tol = 1e-9, plan = "operator"
  )
)
dense_fit <- dense_profile$value
sparse_fit <- sparse_profile$value
if (!isTRUE(dense_fit$converged) || !isTRUE(sparse_fit$converged)) {
  stop("Benchmark solve was not certified.", call. = FALSE)
}

signal <- cbind(sin(seq_len(n) / 17), cos(seq_len(n) / 23))
dense_apply <- peak_rss_call(
  transport_plan_adjoint(rfugw_plan(dense_fit), signal)
)
sparse_apply <- peak_rss_call(
  transport_plan_adjoint(rfugw_plan(sparse_fit), signal)
)
dense_materialize <- peak_rss_call(
  transport_plan_materialize(rfugw_plan(dense_fit))
)
sparse_materialize <- peak_rss_call(
  transport_plan_materialize(rfugw_plan(sparse_fit))
)

dense_plan <- dense_materialize$value
sparse_plan <- sparse_materialize$value
plan_error <- max(abs(dense_plan - sparse_plan))
objective_error <- abs(
  dense_fit$regularized_objective - sparse_fit$regularized_objective
)

row <- function(label, fit, profile, apply_profile, materialize_profile, input) {
  support_working_bytes <- 24 * fit$support_size +
    80 * (n + n)
  solve_problem_bytes <- as.numeric(object.size(input)) + support_working_bytes
  data.frame(
    backend = label,
    n_source = n,
    n_target = n,
    support_size = fit$support_size,
    support_density = fit$support_density,
    setup_seconds = fit$runtime_provenance$timing$setup_seconds,
    solve_seconds = fit$runtime_provenance$timing$solve_seconds,
    certificate_seconds = fit$runtime_provenance$timing$certificate_seconds,
    end_to_end_seconds = profile$elapsed_seconds,
    apply_seconds = apply_profile$elapsed_seconds,
    materialize_seconds = materialize_profile$elapsed_seconds,
    solve_peak_rss_bytes = profile$peak_rss_bytes,
    apply_peak_rss_bytes = apply_profile$peak_rss_bytes,
    materialize_peak_rss_bytes = materialize_profile$peak_rss_bytes,
    estimated_solve_problem_bytes = solve_problem_bytes,
    estimated_apply_problem_bytes = solve_problem_bytes +
      as.numeric(object.size(signal)) + 8 * n * ncol(signal),
    estimated_materialize_problem_bytes = solve_problem_bytes + 8 * n * n,
    iterations = fit$iterations,
    kkt_residual = fit$kkt_residual,
    primal_dual_gap = fit$primal_dual_gap,
    transported_mass = fit$transported_mass,
    dense_plan_max_abs_error = plan_error,
    dense_objective_abs_error = objective_error,
    stringsAsFactors = FALSE
  )
}

result <- rbind(
  row("dense_large_forbidden_cost", dense_fit, dense_profile,
      dense_apply, dense_materialize, dense_cost),
  row("sparse_exact_support", sparse_fit, sparse_profile,
      sparse_apply, sparse_materialize, edges)
)
output <- file.path("inst", "bench", "uot-ti-baseline.csv")
write.csv(result, output, row.names = FALSE)
print(result)
