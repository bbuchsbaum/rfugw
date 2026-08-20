#' Extract the coupling from an rfugw result
#'
#' @param x An `rfugw_result` or a list with `plan` / `pi_samp`.
#' @return The coupling matrix.
#' @export
rfugw_plan <- function(x) {
  if (!is.null(x$plan)) {
    return(x$plan)
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
  list(
    residual = x$residual %||% x$error,
    row_residual = x$row_residual,
    col_residual = x$col_residual,
    mass = x$mass,
    mass_residual = x$mass_residual,
    mass_target = x$mass_target,
    mass_certified = x$mass_certified,
    mass_certification = x$mass_certification,
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
