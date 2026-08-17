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
#' @param method `"scaling"` or `"log"`.
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
    method = c("scaling", "log"),
    max_iter = 1000L,
    tol = 1e-9) {
  method <- match.arg(method)
  dat <- .prepare_linear_ot(M, p, q)
  epsilon <- .validate_positive_scalar(epsilon, "epsilon")
  max_iter <- .validate_count(max_iter, "max_iter")
  tol <- .validate_positive_scalar(tol, "tol")
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
  residual <- if (!is.null(out$error)) out$error else Inf
  .attach_solver_diagnostics(
    out,
    residual = residual,
    converged = is.finite(residual) && residual <= tol,
    iterations = out$iterations,
    max_iter = max_iter,
    p = dat$p,
    q = dat$q,
    plan = out$plan
  )
}

#' Exact balanced linear transport
#'
#' Network-simplex / assignment backend. Reported `ot_dist` is `<M, plan>`.
#'
#' @inheritParams ot_sinkhorn
#' @param max_iter Maximum simplex iterations.
#' @param tol Optimality tolerance.
#' @return An `rfugw_result` with `plan`, `ot_dist`, `status`, and residuals.
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
  residual <- if (isTRUE(out$lp_ok)) 0 else Inf
  .attach_solver_diagnostics(
    out,
    residual = residual,
    converged = isTRUE(out$lp_ok),
    iterations = out$iterations,
    max_iter = max_iter,
    p = dat$p,
    q = dat$q,
    plan = out$plan,
    lp_ok = isTRUE(out$lp_ok)
  )
}

#' KL-unbalanced entropic optimal transport
#'
#' @inheritParams ot_sinkhorn
#' @param rho Marginal KL penalties, length 1 or 2.
#' @param init_plan Optional warm start.
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
    init_plan = NULL) {
  dat <- .prepare_linear_ot(M, p, q)
  epsilon <- .validate_positive_scalar(epsilon, "epsilon")
  max_iter <- .validate_count(max_iter, "max_iter")
  tol <- .validate_positive_scalar(tol, "tol")
  if (length(rho) == 1L) rho <- c(rho, rho)
  if (length(rho) != 2L || any(!is.finite(rho)) || any(rho <= 0)) {
    stop("`rho` must be one or two finite positive numbers.", call. = FALSE)
  }
  init_plan <- .validate_optional_init_plan(init_plan, length(dat$p), length(dat$q), "init_plan")
  out <- cpp_ot_sinkhorn_unbalanced(
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
  out$formulation <- "ot_sinkhorn_unbalanced"
  out$regularization <- epsilon
  out$backend <- "cpp"
  residual <- if (!is.null(out$error)) out$error else Inf
  .attach_solver_diagnostics(
    out,
    residual = residual,
    converged = is.finite(residual) && residual <= tol,
    iterations = out$iterations,
    max_iter = max_iter,
    plan = out$plan
  )
}
