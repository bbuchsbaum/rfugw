.prepare_linear_ot <- function(M, p, q) {
  M <- .validate_finite_matrix(M, "M")
  ns <- nrow(M)
  nt <- ncol(M)
  if (is.null(p)) p <- rep(1 / ns, ns)
  if (is.null(q)) q <- rep(1 / nt, nt)
  p <- .assert_prob(p, ns, "p")
  q <- .assert_prob(q, nt, "q")
  list(M = unname(M), p = unname(p), q = unname(q))
}

#' Balanced entropic optimal transport
#'
#' Scaling or log-domain Sinkhorn for a linear cost. Reported `ot_dist` is
#' the unregularized `<M, plan>` cost.
#'
#' @param M Cost matrix (`ns x nt`).
#' @param p Source weights (default uniform). Renormalized to sum 1.
#' @param q Target weights (default uniform). Renormalized to sum 1.
#' @param epsilon Positive entropic regularization.
#' @param method `"scaling"`, `"log"`, or `"auto"`. Auto uses scaling only
#'   when the maximum scaled exponent magnitude and span are at most 500;
#'   otherwise it selects the genuine log-domain implementation.
#' @param max_iter Maximum Sinkhorn iterations.
#' @param tol Marginal residual tolerance.
#' @return An `rfugw_result` with `plan`, `ot_dist`, `status`, and residuals.
#' @examples
#' M <- matrix(c(0, 1, 1, 0), 2, 2)
#' out <- ot_sinkhorn(M, epsilon = 0.1)
#' out$status
#' rfugw_value(out)
#' @export
ot_sinkhorn <- function(
    M,
    p = NULL,
    q = NULL,
    epsilon = 0.05,
    method = c("scaling", "log", "auto"),
    max_iter = 1000L,
    tol = 1e-9) {
  requested_method <- match.arg(method)
  dat <- .prepare_linear_ot(M, p, q)
  epsilon <- .validate_positive_scalar(epsilon, "epsilon")
  max_iter <- .validate_count(max_iter, "max_iter")
  tol <- .validate_positive_scalar(tol, "tol")
  dispatch <- .select_sinkhorn_method(
    requested_method, dat$M, epsilon, precision = "double",
    context = "Balanced Sinkhorn"
  )
  method <- dispatch$effective
  out <- cpp_ot_sinkhorn(
    M = dat$M,
    p = dat$p,
    q = dat$q,
    epsilon = epsilon,
    max_iter = max_iter,
    tol = tol,
    use_log = identical(method, "log")
  )
  out$ot_dist <- as.numeric(out$ot_dist)
  out$formulation <- "ot_sinkhorn"
  out$regularization <- epsilon
  out$backend <- if (identical(method, "log")) "cpp_log" else "cpp_scaling"
  out$requested_sinkhorn_method <- dispatch$requested
  out$effective_sinkhorn_method <- dispatch$effective
  out$sinkhorn_backend_transition <- dispatch$transition
  out$sinkhorn_dispatch_reason <- dispatch$reason
  out$sinkhorn_dynamic_range <- dispatch$metric
  out$sinkhorn_scaling_threshold <- dispatch$threshold
  residual <- if (!is.null(out$error)) out$error else Inf
  ans <- .attach_solver_diagnostics(
    out,
    residual = residual,
    converged = is.finite(residual) && residual <= tol,
    iterations = out$iterations,
    max_iter = max_iter,
    p = dat$p,
    q = dat$q,
    plan = out$plan,
    feasibility = "balanced",
    feasibility_tol = tol,
    objective_recomputed = ot_linear_cost(dat$M, out$plan)
  )
  ans$termination_reason <- .termination_reason_from_result(ans, max_iter)
  ans
}

#' Exact balanced linear transport
#'
#' Network-simplex / assignment backend. Reported `ot_dist` is `<M, plan>`.
#' An optimal result is accepted only when marginal feasibility, nonbasic
#' reduced costs, and the primal-dual gap pass the reported certificate
#' tolerances. Dual potentials and all certificate components are returned.
#'
#' @inheritParams ot_sinkhorn
#' @param max_iter Maximum simplex iterations.
#' @param tol Optimality tolerance.
#' @return An `rfugw_result` with `plan`, `ot_dist`, `status`, exact termination
#'   reason, source and target dual potentials, and optimality-certificate
#'   diagnostics.
#' @examples
#' M <- matrix(c(0, 2, 2, 0), 2, 2)
#' out <- ot_emd(M)
#' out$converged
#' ot_validate_plan(out, c(0.5, 0.5), c(0.5, 0.5))
#' @export
ot_emd <- function(
    M,
    p = NULL,
    q = NULL,
    max_iter = 20000L,
    tol = 1e-12) {
  dat <- .prepare_linear_ot(M, p, q)
  max_iter <- .validate_count(max_iter, "max_iter")
  tol <- .validate_positive_scalar(tol, "tol")
  out <- cpp_ot_emd(
    M = dat$M,
    p = dat$p,
    q = dat$q,
    max_iter = max_iter,
    tol = tol
  )
  out$formulation <- "ot_emd"
  out$regularization <- 0
  out$backend <- "cpp_transport"
  residual <- out$error
  ans <- .attach_solver_diagnostics(
    out,
    residual = residual,
    converged = isTRUE(out$lp_ok),
    iterations = out$iterations,
    max_iter = max_iter,
    p = dat$p,
    q = dat$q,
    plan = out$plan,
    feasibility = "balanced",
    feasibility_tol = out$feasibility_tolerance,
    objective_recomputed = ot_linear_cost(dat$M, out$plan),
    lp_ok = isTRUE(out$lp_ok)
  )
  if (!isTRUE(out$lp_ok)) {
    ans$status <- out$termination_reason
    ans$converged <- FALSE
  }
  ans
}

#' Exact partial linear optimal transport
#'
#' Solves the nonnegative-cost partial transport problem
#' `min <M, G>` subject to `rowSums(G) <= p`, `colSums(G) <= q`, and
#' `sum(G) = mass`. Weights are normalized to probability vectors, so `mass`
#' is in `[0, 1]`. The implementation reduces the problem to the certified
#' exact transport primitive with one dummy source and target; it does not add
#' a second simplex implementation.
#'
#' @inheritParams ot_emd
#' @param mass Transported probability mass in `[0, 1]`.
#' @return An `rfugw_result` with `plan`, `ot_dist`, `partial_ot_dist`, partial
#'   feasibility and mass certificates, and the exact certificate for the
#'   equivalent augmented transport problem. Augmented dual potentials include
#'   the final dummy coordinate.
#' @examples
#' M <- matrix(c(0, 3, 2, 0), 2, 2)
#' out <- ot_partial_emd(M, mass = 0.5)
#' out$status
#' sum(rfugw_plan(out))
#' @export
ot_partial_emd <- function(
    M,
    p = NULL,
    q = NULL,
    mass = 1,
    max_iter = 20000L,
    tol = 1e-12) {
  dat <- .prepare_linear_ot(M, p, q)
  if (any(dat$M < 0)) {
    stop("`M` must be nonnegative for exact partial transport.", call. = FALSE)
  }
  if (!is.numeric(mass) || length(mass) != 1L ||
      !is.finite(mass) || mass < 0 || mass > 1) {
    stop("`mass` must be one finite number in [0, 1].", call. = FALSE)
  }
  max_iter <- .validate_count(max_iter, "max_iter")
  tol <- .validate_positive_scalar(tol, "tol")

  ns <- nrow(dat$M)
  nt <- ncol(dat$M)
  penalty <- 2 * max(dat$M, 0) + 1
  augmented_cost <- rbind(
    cbind(dat$M, rep(0, ns)),
    c(rep(0, nt), penalty)
  )
  augmented_p <- c(dat$p, 1 - mass)
  augmented_q <- c(dat$q, 1 - mass)
  native <- cpp_ot_emd(
    M = augmented_cost,
    p = augmented_p,
    q = augmented_q,
    max_iter = max_iter,
    tol = tol
  )
  plan <- native$plan[seq_len(ns), seq_len(nt), drop = FALSE]
  value <- ot_linear_cost(dat$M, plan)
  native$plan <- plan
  native$ot_dist <- value
  native$partial_ot_dist <- value
  native$formulation <- "ot_partial_emd"
  native$regularization <- 0
  native$backend <- "cpp_transport_dummy_reduction"
  native$transported_mass_target <- mass
  native$transported_mass_defaulted <- FALSE
  native$dummy_penalty <- penalty
  native$augmented_shape <- c(ns + 1L, nt + 1L)
  residual <- max(
    native$error,
    max(c(rowSums(plan) - dat$p, 0)),
    max(c(colSums(plan) - dat$q, 0)),
    abs(sum(plan) - mass)
  )
  ans <- .attach_solver_diagnostics(
    native,
    residual = residual,
    converged = isTRUE(native$lp_ok),
    iterations = native$iterations,
    max_iter = max_iter,
    p = dat$p,
    q = dat$q,
    plan = plan,
    feasibility = "partial",
    feasibility_tol = native$feasibility_tolerance,
    mass_target = mass,
    objective_recomputed = value,
    lp_ok = isTRUE(native$lp_ok)
  )
  if (!isTRUE(native$lp_ok)) {
    ans$status <- native$termination_reason
    ans$converged <- FALSE
  }
  ans$termination_reason <- native$termination_reason
  ans
}

#' KL-unbalanced entropic optimal transport
#'
#' @inheritParams ot_sinkhorn
#' @param rho Marginal KL penalties, length 1 or 2.
#' @param init_plan Optional warm start.
#' @param method `"scaling"`, genuine `"log"`, or `"auto"`. Auto uses the
#'   documented dynamic-range criterion from `ot_sinkhorn()`.
#' @return An `rfugw_result` with `plan`, `ot_dist`, `mass`, and status.
#' @examples
#' M <- matrix(c(0, 1, 1, 0), 2, 2)
#' out <- ot_sinkhorn_unbalanced(M, epsilon = 0.1, rho = 2)
#' out$mass
#' @export
ot_sinkhorn_unbalanced <- function(
    M,
    p = NULL,
    q = NULL,
    epsilon = 0.05,
    rho = 10,
    max_iter = 500L,
    tol = 1e-7,
    init_plan = NULL,
    method = c("scaling", "log", "auto")) {
  requested_method <- match.arg(method)
  dat <- .prepare_linear_ot(M, p, q)
  epsilon <- .validate_positive_scalar(epsilon, "epsilon")
  max_iter <- .validate_count(max_iter, "max_iter")
  tol <- .validate_positive_scalar(tol, "tol")
  if (length(rho) == 1L) rho <- c(rho, rho)
  if (length(rho) != 2L || any(!is.finite(rho)) || any(rho <= 0)) {
    stop("`rho` must be one or two finite positive numbers.", call. = FALSE)
  }
  init_plan <- .validate_optional_init_plan(init_plan, length(dat$p), length(dat$q), "init_plan")
  dispatch <- .select_sinkhorn_method(
    requested_method, dat$M, epsilon, precision = "double",
    context = "Unbalanced Sinkhorn"
  )
  out <- if (identical(dispatch$effective, "log")) {
    .sinkhorn_unbalanced_log(
      M = dat$M, a = dat$p, b = dat$q, epsilon = epsilon, rho = rho,
      max_iter = max_iter, tol = tol, init_plan = init_plan
    )
  } else {
    cpp_ot_sinkhorn_unbalanced(
      M = dat$M,
      a = dat$p,
      b = dat$q,
      epsilon = epsilon,
      rho1 = rho[[1]],
      rho2 = rho[[2]],
      max_iter = max_iter,
      tol = tol,
      init_plan = init_plan
    )
  }
  out$formulation <- "ot_sinkhorn_unbalanced"
  out$regularization <- epsilon
  out$backend <- if (identical(dispatch$effective, "log")) "r_log" else "cpp_scaling"
  out$requested_sinkhorn_method <- dispatch$requested
  out$effective_sinkhorn_method <- dispatch$effective
  out$sinkhorn_backend_transition <- dispatch$transition
  out$sinkhorn_dispatch_reason <- dispatch$reason
  out$sinkhorn_dynamic_range <- dispatch$metric
  out$sinkhorn_scaling_threshold <- dispatch$threshold
  residual <- if (!is.null(out$error)) out$error else Inf
  ans <- .attach_solver_diagnostics(
    out,
    residual = residual,
    converged = is.finite(residual) && residual <= tol,
    iterations = out$iterations,
    max_iter = max_iter,
    plan = out$plan,
    inner_residual = out$inner_residual,
    max_inner_residual = out$max_inner_residual,
    inner_iterations = out$inner_iterations,
    inner_converged = out$inner_converged,
    inner_status = out$inner_status,
    feasibility = "unbalanced",
    feasibility_tol = tol,
    objective_recomputed = ot_linear_cost(dat$M, out$plan)
  )
  ans$termination_reason <- .termination_reason_from_result(ans, max_iter)
  ans
}
