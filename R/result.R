#' Extract the coupling from an rfugw result
#'
#' @param x An `rfugw_result`, transport plan, or a list with `plan` /
#'   `pi_samp`.
#' @param materialize If `TRUE`, explicitly convert a transport-plan
#'   representation to a dense matrix. Dense legacy results are unchanged.
#' @return The stored coupling matrix/representation, or a dense matrix when
#'   explicitly requested.
#' @export
rfugw_plan <- function(x, materialize = FALSE) {
  if (!is.logical(materialize) || length(materialize) != 1L ||
      is.na(materialize)) {
    stop("`materialize` must be TRUE or FALSE.", call. = FALSE)
  }
  if (.is_transport_plan(x)) {
    return(if (materialize) transport_plan_materialize(x) else x)
  }
  if (!is.null(x$plan)) {
    plan <- x$plan
    return(if (materialize && .is_transport_plan(plan)) {
      transport_plan_materialize(plan)
    } else {
      plan
    })
  }
  if (!is.null(x$pi_samp)) {
    return(x$pi_samp)
  }
  stop("`x` does not contain a coupling.", call. = FALSE)
}

#' Extract the documented objective value
#'
#' @param x An `rfugw_result` or compatible list.
#' @return Numeric scalar objective.
#' @export
rfugw_value <- function(x) {
  if (!is.null(x$barycenter_objective)) {
    return(as.numeric(x$barycenter_objective)[1])
  }
  if (identical(x$formulation, "ot_sinkhorn_unbalanced_ti") &&
      !is.null(x$regularized_objective)) {
    return(as.numeric(x$regularized_objective)[1])
  }
  if (!is.null(x$partial_sinkhorn_objective)) {
    return(as.numeric(x$partial_sinkhorn_objective)[1])
  }
  if (!is.null(x$penalized_partial_objective)) {
    return(as.numeric(x$penalized_partial_objective)[1])
  }
  if (!is.null(x$sinkhorn_divergence)) {
    return(as.numeric(x$sinkhorn_divergence)[1])
  }
  if (!is.null(x$wasserstein_value)) {
    return(as.numeric(x$wasserstein_value)[1])
  }
  for (nm in c("ot_dist", "fgw_dist", "gw_dist", "fugw_cost", "ucoot_cost",
               "srgw_dist", "srfgw_dist", "partial_gw_dist", "partial_fgw_dist")) {
    if (!is.null(x[[nm]])) {
      return(as.numeric(x[[nm]])[1])
    }
  }
  stop("`x` does not contain a documented objective field.", call. = FALSE)
}

#' Extract solver status
#'
#' @param x An `rfugw_result`.
#' @return Character status string.
#' @export
rfugw_status <- function(x) {
  if (is.null(x$status)) {
    stop("`x` does not contain `status`.", call. = FALSE)
  }
  x$status
}

#' Extract residual diagnostics
#'
#' @param x An `rfugw_result`.
#' @return A list with stopping, marginal, feasibility, objective, and nested
#'   solver certificate fields.
#' @export
rfugw_residuals <- function(x) {
  plan_representation <- tryCatch(
    transport_plan_representation(rfugw_plan(x)),
    error = function(e) NULL
  )
  list(
    residual = x$residual %||% x$error,
    row_residual = x$row_residual,
    col_residual = x$col_residual,
    mass = x$mass,
    mass_residual = x$mass_residual,
    mass_target = x$mass_target,
    mass_certified = x$mass_certified,
    mass_certification = x$mass_certification,
    transported_mass = x$transported_mass,
    original_source_mass = x$original_source_mass,
    original_target_mass = x$original_target_mass,
    effective_source_mass = x$effective_source_mass,
    effective_target_mass = x$effective_target_mass,
    normalization = x$measure_normalization,
    source_marginal_kl = x$source_marginal_kl,
    target_marginal_kl = x$target_marginal_kl,
    plan_product_kl = x$plan_product_kl,
    regularized_objective = x$regularized_objective,
    uot_certificate = x$uot_certificate,
    fixed_point_residual = x$fixed_point_residual,
    fixed_point_tolerance = x$fixed_point_tolerance,
    fixed_point_consistent = x$fixed_point_consistent,
    source_kkt_residual = x$source_kkt_residual,
    target_kkt_residual = x$target_kkt_residual,
    kkt_residual = x$kkt_residual,
    kkt_tolerance = x$kkt_tolerance,
    kkt_consistent = x$kkt_consistent,
    primal_dual_gap = x$primal_dual_gap,
    gauge_residual = x$gauge_residual,
    gauge_invariant = x$gauge_invariant,
    support_certificate = x$support_certificate,
    wasserstein_power = x$wasserstein_power,
    input_cost_power = x$input_cost_power,
    effective_cost_power = x$effective_cost_power,
    value_kind = x$value_kind,
    value_certified = x$value_certified,
    value_certification = x$value_certification,
    sinkhorn_divergence = x$sinkhorn_divergence,
    component_status = x$component_status,
    component_residuals = x$component_residuals,
    component_converged = x$component_converged,
    discarded_source_mass = x$discarded_source_mass,
    discarded_target_mass = x$discarded_target_mass,
    discard_penalty_term = x$discard_penalty_term,
    entropy_minus_one = x$entropy_minus_one,
    weighted_entropy_minus_one = x$weighted_entropy_minus_one,
    stationarity_residual = x$stationarity_residual,
    complementarity_residual = x$complementarity_residual,
    dual_feasibility_residual = x$dual_feasibility_residual,
    kkt_tolerance = x$kkt_tolerance,
    frank_wolfe_gap = x$frank_wolfe_gap,
    frank_wolfe_gap_tolerance = x$frank_wolfe_gap_tolerance,
    stationarity_consistent = x$stationarity_consistent,
    line_search_residual = x$line_search_residual,
    line_search_tolerance = x$line_search_tolerance,
    line_search_consistent = x$line_search_consistent,
    plan_representation = x$plan_representation %||% plan_representation,
    pruning_lost_mass = x$pruning_lost_mass,
    certificate_invalidated_by_pruning =
      x$certificate_invalidated_by_pruning,
    feasibility = x$feasibility,
    feasibility_residual = x$feasibility_residual,
    feasibility_tolerance = x$feasibility_tolerance,
    feasible = x$feasible,
    inner_residual = x$inner_residual,
    max_inner_residual = x$max_inner_residual,
    inner_converged = x$inner_converged,
    inner_status = x$inner_status,
    objective_recomputed = x$objective_recomputed,
    objective_residual = x$objective_residual,
    objective_tolerance = x$objective_tolerance,
    objective_consistent = x$objective_consistent,
    objective_components_consistent = x$objective_components_consistent,
    runtime_provenance = x$runtime_provenance
  )
}

#' Print an rfugw result
#'
#' @param x An `rfugw_result`.
#' @param ... Ignored.
#' @export
print.rfugw_result <- function(x, ...) {
  value <- tryCatch(rfugw_value(x), error = function(e) NA_real_)
  cat("<rfugw_result>\n")
  cat(sprintf("  formulation: %s\n", x$formulation %||% "unknown"))
  cat(sprintf("  backend:     %s\n", x$backend %||% "unknown"))
  cat(sprintf("  status:      %s\n", x$status %||% "unknown"))
  representation <- tryCatch(
    transport_plan_representation(rfugw_plan(x))$representation,
    error = function(e) NULL
  )
  if (!is.null(representation) &&
      !identical(representation, "dense_materialized")) {
    cat(sprintf("  plan:        %s\n", representation))
  }
  cat(sprintf("  value:       %s\n", format(value, digits = 6)))
  cat(sprintf("  iterations:  %s / %s\n", x$iterations %||% NA, x$max_iter %||% NA))
  cat(sprintf("  residual:    %s\n", format(x$residual %||% x$error, digits = 4)))
  if (!is.null(x$row_residual) && is.finite(x$row_residual)) {
    cat(sprintf("  row/col res: %s / %s\n",
                format(x$row_residual, digits = 4),
                format(x$col_residual, digits = 4)))
  }
  if (!is.null(x$regularization) && is.finite(x$regularization)) {
    cat(sprintf("  regularizer: %s\n", format(x$regularization, digits = 4)))
  }
  invisible(x)
}

#' Summarize an rfugw result
#'
#' @param object An `rfugw_result`.
#' @param ... Ignored.
#' @export
summary.rfugw_result <- function(object, ...) {
  out <- list(
    formulation = object$formulation,
    backend = object$backend,
    status = object$status,
    converged = object$converged,
    value = tryCatch(rfugw_value(object), error = function(e) NA_real_),
    iterations = object$iterations,
    max_iter = object$max_iter,
    residuals = rfugw_residuals(object),
    plan_dim = dim(rfugw_plan(object))
  )
  class(out) <- "summary.rfugw_result"
  out
}

#' Print an rfugw result summary
#'
#' @param x A `summary.rfugw_result`.
#' @param ... Ignored.
#' @export
print.summary.rfugw_result <- function(x, ...) {
  cat("rfugw result summary\n")
  cat(sprintf("  %s via %s: %s\n", x$formulation, x$backend, x$status))
  cat(sprintf("  value=%s  iterations=%s/%s  plan=%s x %s\n",
              format(x$value, digits = 6),
              x$iterations, x$max_iter,
              x$plan_dim[1], x$plan_dim[2]))
  invisible(x)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
