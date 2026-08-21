.as_plan <- function(plan, name = "plan") {
  if (is.list(plan)) {
    plan <- rfugw_plan(plan)
  }
  .assert_matrix(plan, name)
  if (nrow(plan) == 0L || ncol(plan) == 0L) {
    stop(sprintf("`%s` must be nonempty.", name), call. = FALSE)
  }
  if (any(!is.finite(plan))) {
    stop(sprintf("`%s` must be finite.", name), call. = FALSE)
  }
  if (any(plan < 0)) {
    stop(sprintf("`%s` must be nonnegative.", name), call. = FALSE)
  }
  plan
}

#' Validate a transport plan
#'
#' Checks shape, nonnegativity, total mass, and optional marginal constraints.
#' Accepts a raw matrix or an `rfugw_result`.
#'
#' @param plan Coupling matrix or result object.
#' @param p Optional source weights.
#' @param q Optional target weights.
#' @param mass Optional required total mass.
#' @param marginals `"balanced"` (row/col match `p`/`q`), `"partial"`
#'   (row/col do not exceed `p`/`q` and total mass matches), or `"relaxed"`
#'   (nonnegativity and finiteness only).
#' @param tol Residual tolerance.
#' @return The plan, invisibly, after validation.
#' @examples
#' M <- matrix(c(0, 1, 1, 0), 2, 2)
#' out <- ot_emd(M)
#' ot_validate_plan(out, c(0.5, 0.5), c(0.5, 0.5))
#' ot_linear_cost(M, out)
#' ot_entropy(out)
#' @export
ot_validate_plan <- function(plan,
                             p = NULL,
                             q = NULL,
                             mass = NULL,
                             marginals = c("balanced", "partial", "relaxed"),
                             tol = 1e-8) {
  marginals <- match.arg(marginals)
  plan <- .as_plan(plan)
  if (!is.null(p)) {
    if (!is.numeric(p) || length(p) != nrow(plan)) {
      stop("`p` must be numeric with length nrow(plan).", call. = FALSE)
    }
  }
  if (!is.null(q)) {
    if (!is.numeric(q) || length(q) != ncol(plan)) {
      stop("`q` must be numeric with length ncol(plan).", call. = FALSE)
    }
  }
  if (!is.null(mass)) {
    if (length(mass) != 1L || !is.finite(mass) || mass < 0) {
      stop("`mass` must be a finite nonnegative scalar.", call. = FALSE)
    }
    if (abs(sum(plan) - mass) > tol) {
      stop(sprintf("`plan` total mass must equal `mass` within %g.", tol), call. = FALSE)
    }
  }
  if (identical(marginals, "balanced")) {
    if (is.null(p) || is.null(q)) {
      stop("`marginals = \"balanced\"` requires `p` and `q`.", call. = FALSE)
    }
    if (max(abs(rowSums(plan) - p)) > tol || max(abs(colSums(plan) - q)) > tol) {
      stop("`plan` does not satisfy the balanced marginals.", call. = FALSE)
    }
  } else if (identical(marginals, "partial")) {
    if (is.null(p) || is.null(q)) {
      stop("`marginals = \"partial\"` requires `p` and `q`.", call. = FALSE)
    }
    if (any(rowSums(plan) - p > tol) || any(colSums(plan) - q > tol)) {
      stop("`plan` exceeds the partial marginal upper bounds.", call. = FALSE)
    }
  }
  invisible(plan)
}

#' Linear transport cost
#'
#' @param M Cost matrix.
#' @param plan Coupling or `rfugw_result`.
#' @return Scalar `<M, plan>`.
#' @export
ot_linear_cost <- function(M, plan) {
  plan <- .as_plan(plan)
  M <- .validate_finite_matrix(M, "M")
  if (!all(dim(M) == dim(plan))) {
    stop("`M` and `plan` must have the same shape.", call. = FALSE)
  }
  sum(M * plan)
}

#' Entropic term of a plan
#'
#' Computes `sum(plan * log(plan))` with the convention `0 log 0 = 0`.
#'
#' @param plan Coupling or `rfugw_result`.
#' @return Numeric scalar.
#' @export
ot_entropy <- function(plan) {
  plan <- .as_plan(plan)
  z <- plan[plan > 0]
  sum(z * log(z))
}

#' Generalized KL divergence of a plan from a product reference
#'
#' @param plan Coupling or `rfugw_result`.
#' @param p Finite nonnegative source reference weights.
#' @param q Finite nonnegative target reference weights.
#' @return Generalized `KL(plan || p %o% q)`, including the mass correction
#'   `-sum(plan) + sum(p %o% q)`. Positive plan mass outside zero reference
#'   support returns `Inf`; zero-over-zero contributes zero.
#' @export
ot_kl <- function(plan, p, q) {
  plan <- .as_plan(plan)
  if (!is.numeric(p) || length(p) != nrow(plan)) {
    stop("`p` must be numeric with length nrow(plan).", call. = FALSE)
  }
  if (!is.numeric(q) || length(q) != ncol(plan)) {
    stop("`q` must be numeric with length ncol(plan).", call. = FALSE)
  }
  if (any(!is.finite(p))) {
    stop("`p` must be finite.", call. = FALSE)
  }
  if (any(!is.finite(q))) {
    stop("`q` must be finite.", call. = FALSE)
  }
  if (any(p < 0)) {
    stop("`p` must be nonnegative.", call. = FALSE)
  }
  if (any(q < 0)) {
    stop("`q` must be nonnegative.", call. = FALSE)
  }
  ref <- tcrossprod(p, q)
  if (any(!is.finite(ref))) {
    stop("The product reference `p %o% q` must be finite.", call. = FALSE)
  }
  if (any(plan > 0 & ref == 0)) {
    return(Inf)
  }
  z <- plan > 0
  sum(plan[z] * log(plan[z] / ref[z])) - sum(plan) + sum(ref)
}

#' Square-loss Gromov-Wasserstein objective
#'
#' Independent evaluator of `sum_ijkl (C1_ij - C2_kl)^2 plan_ik plan_jl`.
#' Uses the standard factored expansion; tests compare it to the O(n^4) form.
#'
#' @param C1 Source structure matrix.
#' @param C2 Target structure matrix.
#' @param plan Coupling or `rfugw_result`.
#' @param symmetric `NULL` auto-detects; `TRUE` requires symmetry.
#' @return Numeric scalar.
#' @export
ot_gw_square <- function(C1, C2, plan, symmetric = NULL) {
  plan <- .as_plan(plan)
  C1 <- .validate_finite_matrix(C1, "C1", square = TRUE)
  C2 <- .validate_finite_matrix(C2, "C2", square = TRUE)
  if (nrow(C1) != nrow(plan) || nrow(C2) != ncol(plan)) {
    stop("`plan` must have shape nrow(C1) x nrow(C2).", call. = FALSE)
  }
  symmetric <- .resolve_symmetric(symmetric, C1, C2)
  p <- rowSums(plan)
  q <- colSums(plan)
  ns <- nrow(C1)
  nt <- nrow(C2)
  constC <- tcrossprod(as.vector((C1^2) %*% p), rep(1, nt)) +
    tcrossprod(rep(1, ns), as.vector(q %*% t(C2^2)))
  tens <- constC - C1 %*% plan %*% t(2 * C2)
  if (symmetric) {
    return(sum(tens * plan))
  }
  constCt <- tcrossprod(as.vector((t(C1)^2) %*% p), rep(1, nt)) +
    tcrossprod(rep(1, ns), as.vector(q %*% t(t(C2)^2)))
  tenst <- constCt - t(C1) %*% plan %*% t(2 * t(C2))
  0.5 * (sum(tens * plan) + sum(tenst * plan))
}

#' Square-loss fused Gromov-Wasserstein objective
#'
#' Unregularized FGW value: `(1 - alpha) * <M, plan> + alpha * GW(C1, C2, plan)`.
#'
#' @inheritParams ot_gw_square
#' @param M Feature cost matrix.
#' @param alpha Feature/structure trade-off in `[0, 1]`.
#' @return Numeric scalar.
#' @export
ot_fgw_square <- function(M, C1, C2, plan, alpha = 0.5, symmetric = NULL) {
  alpha <- .validate_alpha(alpha)
  plan <- .as_plan(plan)
  M <- .validate_finite_matrix(M, "M")
  if (!all(dim(M) == dim(plan))) {
    stop("`M` and `plan` must have the same shape.", call. = FALSE)
  }
  (1 - alpha) * ot_linear_cost(M, plan) + alpha * ot_gw_square(C1, C2, plan, symmetric)
}

#' Barycentric projection of a coupling
#'
#' Maps source points through `plan` onto the target support, or the reverse.
#' Rows (or columns) with zero transported mass return `NaN` by default.
#'
#' @param plan Coupling or `rfugw_result`.
#' @param points Point matrix on the destination support (`nt x d` for
#'   `source_to_target`, `ns x d` for `target_to_source`).
#' @param orientation `"source_to_target"` or `"target_to_source"`.
#' @param zero_mass `"nan"` or `"zero"` for empty-mass rows.
#' @return Projected point matrix.
#' @examples
#' plan <- diag(c(1, 0))
#' points <- matrix(c(0, 1, 2, 3), 2, 2, byrow = TRUE)
#' ot_barycentric_project(plan, points)
#' @export
ot_barycentric_project <- function(plan,
                                   points,
                                   orientation = c("source_to_target", "target_to_source"),
                                   zero_mass = c("nan", "zero")) {
  orientation <- match.arg(orientation)
  zero_mass <- match.arg(zero_mass)
  plan <- .as_plan(plan)
  points <- .validate_finite_matrix(points, "points")
  if (identical(orientation, "source_to_target")) {
    if (nrow(points) != ncol(plan)) {
      stop("`points` must have nrow equal to ncol(plan) for source_to_target.", call. = FALSE)
    }
    mass <- rowSums(plan)
    proj <- plan %*% points
  } else {
    if (nrow(points) != nrow(plan)) {
      stop("`points` must have nrow equal to nrow(plan) for target_to_source.", call. = FALSE)
    }
    mass <- colSums(plan)
    proj <- t(plan) %*% points
  }
  out <- matrix(if (identical(zero_mass, "nan")) NaN else 0, nrow = length(mass), ncol = ncol(points))
  keep <- mass > 0
  if (any(keep)) {
    out[keep, ] <- sweep(proj[keep, , drop = FALSE], 1L, mass[keep], "/")
  }
  out
}
