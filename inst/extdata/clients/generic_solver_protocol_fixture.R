# Installed-package client conformance fixture. It intentionally uses only
# exported rfugw functions and does not inspect solver-specific list fields.
generic_transport_client <- function(cost, source, target, epsilon, tol,
                                     state = NULL) {
  problem <- rfugw::transport_problem_sinkhorn(
    cost,
    p = source,
    q = target,
    epsilon = epsilon,
    method = "auto",
    max_iter = 10000L,
    tol = tol,
    mass_policy = "probability"
  )
  result <- rfugw::transport_solve(problem, init_state = state)
  list(
    result = result,
    plan = rfugw::rfugw_plan(result),
    state = rfugw::rfugw_state(result),
    masses = rfugw::rfugw_masses(result),
    provenance = rfugw::rfugw_provenance(result),
    status = rfugw::rfugw_status(result),
    residuals = rfugw::rfugw_residuals(result)
  )
}

generic_exact_transport_client <- function(cost, source, target) {
  problem <- rfugw::transport_problem_emd(
    cost,
    p = source,
    q = target,
    mass_policy = "probability"
  )
  result <- rfugw::transport_solve(problem)
  list(
    result = result,
    plan = rfugw::rfugw_plan(result),
    value = rfugw::rfugw_value(result),
    provenance = rfugw::rfugw_provenance(result)
  )
}
