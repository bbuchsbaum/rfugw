# DKGE-shaped client fixture. Cost construction and intensive/extensive field
# semantics stay in DKGE; this adapter consumes only the exported rfugw
# protocol and returns its own client-shaped diagnostics.
dkge_rfugw_sinkhorn <- function(cost, source_reliability, target_reliability,
                                epsilon, tol, state = NULL) {
  problem <- rfugw::transport_problem_sinkhorn(
    cost,
    p = source_reliability,
    q = target_reliability,
    epsilon = epsilon,
    method = "log",
    max_iter = 10000L,
    tol = tol,
    mass_policy = "probability"
  )
  fit <- rfugw::transport_solve(problem, init_state = state)
  residuals <- rfugw::rfugw_residuals(fit)
  list(
    plan = rfugw::rfugw_plan(fit),
    state = rfugw::rfugw_state(fit),
    diagnostics = list(
      converged = identical(rfugw::rfugw_status(fit), "converged"),
      iterations = fit$iterations,
      marginal_error = max(residuals$row_residual, residuals$col_residual),
      warm_started = isTRUE(fit$warm_started),
      warm_start_accepted = isTRUE(fit$warm_start_accepted),
      rejection = fit$warm_start_rejection
    ),
    provenance = rfugw::rfugw_provenance(fit)
  )
}
