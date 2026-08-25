# Temporary downstream conformance fixture for migrating manifoldalign's
# generic TI-Sinkhorn, coupling-free map, and sparse coupling extraction into
# rfugw. The downstream package retains cost construction, fit orchestration,
# and prediction policy.
manifoldalign_rfugw_ti_accept <- function(fit) {
  certified <- isTRUE(fit$converged) &&
    identical(rfugw::rfugw_status(fit), "converged") &&
    isTRUE(fit$feasible) &&
    isTRUE(fit$objective_consistent) &&
    isTRUE(fit$objective_components_consistent) &&
    isTRUE(fit$inner_converged) &&
    isTRUE(fit$kkt_consistent) &&
    isTRUE(fit$fixed_point_consistent) &&
    isTRUE(fit$gauge_invariant) &&
    isTRUE(fit$uot_certificate$objective_consistent)
  if (!certified) {
    stop("rfugw did not certify the translation-invariant KL-UOT result.",
         call. = FALSE)
  }
  fit
}

manifoldalign_rfugw_uot_ti_sinkhorn_kl <- function(
    cost, alpha, beta, epsilon, rho1, rho2,
    max_iter = 2000L, tol = 1e-6,
    plan = c("operator", "sparse")) {
  plan <- match.arg(plan)
  fit <- rfugw::ot_sinkhorn_unbalanced_ti(
    cost = cost,
    p = alpha,
    q = beta,
    epsilon = epsilon,
    rho = c(rho1, rho2),
    max_iter = max_iter,
    tol = tol,
    plan = plan
  )
  manifoldalign_rfugw_ti_accept(fit)
}

manifoldalign_rfugw_uot_apply_map <- function(fit, signal, delta = 1e-8) {
  fit <- manifoldalign_rfugw_ti_accept(fit)
  stopifnot(is.numeric(delta), length(delta) == 1L, is.finite(delta), delta >= 0)
  plan <- rfugw::rfugw_plan(fit)
  target_mass <- rfugw::transport_plan_mass(plan, "target")
  if (is.matrix(signal)) {
    if (ncol(signal) != rfugw::transport_plan_shape(plan)[["source"]]) {
      stop("`signal` must have one column per source node.", call. = FALSE)
    }
    numerator <- t(rfugw::transport_plan_adjoint(plan, t(signal)))
    return(sweep(numerator, 2L, target_mass + delta, "/"))
  }
  if (!is.numeric(signal) ||
      length(signal) != rfugw::transport_plan_shape(plan)[["source"]]) {
    stop("`signal` must have one value per source node.", call. = FALSE)
  }
  as.numeric(rfugw::transport_plan_adjoint(plan, signal) /
    (target_mass + delta))
}

manifoldalign_rfugw_uot_extract_coupling <- function(fit) {
  fit <- manifoldalign_rfugw_ti_accept(fit)
  representation <- rfugw::transport_plan_representation(rfugw::rfugw_plan(fit))
  if (!identical(representation$representation, "edge_list")) {
    stop(
      "Refit with `plan = \"sparse\"` to request explicit coupling extraction.",
      call. = FALSE
    )
  }
  rfugw::rfugw_plan(fit)
}
