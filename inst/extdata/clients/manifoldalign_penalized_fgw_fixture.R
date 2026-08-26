# Minimal downstream conformance fixture for re-enabling a manifold-alignment
# lambda mode. The downstream caller consumes a plan only after rfugw certifies
# variable-mass feasibility, objective consistency, all nested direction solves,
# and Frank-Wolfe stationarity.
manifoldalign_accept_penalized_fgw_result <- function(fit) {
  certified <- isTRUE(fit$converged) &&
    identical(rfugw::rfugw_status(fit), "converged") &&
    isTRUE(fit$feasible) &&
    isTRUE(fit$objective_consistent) &&
    isTRUE(fit$objective_components_consistent) &&
    isTRUE(fit$inner_converged) &&
    isTRUE(fit$stationarity_consistent) &&
    is.finite(fit$frank_wolfe_gap) &&
    fit$frank_wolfe_gap <= fit$frank_wolfe_gap_tolerance
  if (!certified) {
    stop("rfugw did not certify the penalized partial FGW result.", call. = FALSE)
  }
  fit
}

manifoldalign_penalized_fgw_step <- function(
    feature_cost, source_structure, target_structure, lambda,
    source_weights = NULL, target_weights = NULL, alpha = 0.5) {
  fit <- rfugw::penalized_partial_fused_gromov_wasserstein(
    M = feature_cost,
    C1 = source_structure,
    C2 = target_structure,
    discard_penalty = lambda,
    p = source_weights,
    q = target_weights,
    alpha = alpha
  )
  fit <- manifoldalign_accept_penalized_fgw_result(fit)
  list(
    coupling = rfugw::rfugw_plan(fit),
    transported_mass = fit$transported_mass,
    result = fit
  )
}
