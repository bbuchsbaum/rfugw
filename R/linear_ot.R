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

.validate_finite_measure <- function(x, n, name, default) {
  if (is.null(x)) x <- default
  if (!is.numeric(x) || length(x) != n) {
    stop(sprintf("`%s` must be numeric with length %d.", name, n), call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(sprintf("`%s` must be finite.", name), call. = FALSE)
  }
  if (any(x < 0)) {
    stop(sprintf("`%s` must be nonnegative.", name), call. = FALSE)
  }
  total <- sum(x)
  if (!is.finite(total)) {
    stop(sprintf("`%s` must have finite total mass.", name), call. = FALSE)
  }
  unname(as.numeric(x))
}

.prepare_unbalanced_linear_ot <- function(M, p, q, normalization) {
  M <- .validate_finite_matrix(M, "M")
  ns <- nrow(M)
  nt <- ncol(M)
  original_p <- .validate_finite_measure(p, ns, "p", rep(1 / ns, ns))
  original_q <- .validate_finite_measure(q, nt, "q", rep(1 / nt, nt))
  original_mass <- c(source = sum(original_p), target = sum(original_q))

  scale <- switch(
    normalization,
    none = c(source = 1, target = 1),
    separate = {
      if (any(original_mass <= 0)) {
        stop(
          paste0(
            "`normalization = \"separate\"` requires positive source and ",
            "target mass; use `normalization = \"none\"` for zero measures."
          ),
          call. = FALSE
        )
      }
      1 / original_mass
    },
    joint = {
      total <- sum(original_mass)
      if (total <= 0) {
        stop(
          "`normalization = \"joint\"` requires positive combined mass; use `normalization = \"none\"` for two zero measures.",
          call. = FALSE
        )
      }
      rep(2 / total, 2L)
    }
  )
  names(scale) <- c("source", "target")
  effective_p <- original_p * scale[["source"]]
  effective_q <- original_q * scale[["target"]]
  if (any(!is.finite(tcrossprod(effective_p, effective_q)))) {
    stop(
      "The effective product reference `p %o% q` must be finite; rescale the input measures explicitly.",
      call. = FALSE
    )
  }

  list(
    M = unname(M),
    p = effective_p,
    q = effective_q,
    original_p = original_p,
    original_q = original_q,
    original_mass = original_mass,
    effective_mass = c(source = sum(effective_p), target = sum(effective_q)),
    scale = scale,
    normalization = normalization
  )
}

.generalized_kl_vector <- function(x, reference) {
  if (any(x > 0 & reference == 0)) return(Inf)
  positive <- x > 0
  sum(x[positive] * log(x[positive] / reference[positive])) -
    sum(x) + sum(reference)
}

.attach_unbalanced_measure_contract <- function(out, dat, rho, epsilon, tol) {
  plan <- out$plan
  source_marginal <- rowSums(plan)
  target_marginal <- colSums(plan)
  source_kl <- .generalized_kl_vector(source_marginal, dat$p)
  target_kl <- .generalized_kl_vector(target_marginal, dat$q)
  plan_kl <- ot_kl(plan, dat$p, dat$q)
  transport_cost <- ot_linear_cost(dat$M, plan)
  regularized_objective <- transport_cost +
    rho[[1L]] * source_kl + rho[[2L]] * target_kl + epsilon * plan_kl

  out$normalization <- dat$normalization
  out$measure_normalization <- dat$normalization
  out$source_measure_original <- dat$original_p
  out$target_measure_original <- dat$original_q
  out$source_measure_effective <- dat$p
  out$target_measure_effective <- dat$q
  out$original_source_mass <- unname(dat$original_mass[["source"]])
  out$original_target_mass <- unname(dat$original_mass[["target"]])
  out$effective_source_mass <- unname(dat$effective_mass[["source"]])
  out$effective_target_mass <- unname(dat$effective_mass[["target"]])
  out$source_measure_scale <- unname(dat$scale[["source"]])
  out$target_measure_scale <- unname(dat$scale[["target"]])
  out$source_marginal_mass <- sum(source_marginal)
  out$target_marginal_mass <- sum(target_marginal)
  out$transported_mass <- sum(plan)
  out$source_marginal_kl <- source_kl
  out$target_marginal_kl <- target_kl
  out$plan_product_kl <- plan_kl
  out$regularized_objective <- regularized_objective
  out$uot_objective_components <- list(
    transport = transport_cost,
    source_marginal_kl = source_kl,
    target_marginal_kl = target_kl,
    plan_product_kl = plan_kl,
    weighted_source_marginal_kl = rho[[1L]] * source_kl,
    weighted_target_marginal_kl = rho[[2L]] * target_kl,
    weighted_plan_product_kl = epsilon * plan_kl,
    total = regularized_objective
  )
  out$uot_certificate <- list(
    fixed_point_residual = out$residual,
    tolerance = tol,
    converged = isTRUE(out$converged),
    finite_plan = all(is.finite(plan)),
    nonnegative_plan = all(plan >= 0)
  )
  out
}

.canonicalize_balanced_duals <- function(source, target, p) {
  shift <- sum(p * source) / sum(p)
  list(source = source - shift, target = target + shift)
}

.log_sum_exp <- function(x) {
  maximum <- max(x)
  if (!is.finite(maximum)) return(maximum)
  maximum + log(sum(exp(x - maximum)))
}

.attach_balanced_regularized_contract <- function(out, dat, epsilon, tol) {
  plan <- out$plan
  transport <- ot_linear_cost(dat$M, plan)
  product_kl <- ot_kl(plan, dat$p, dat$q)
  primal <- transport + epsilon * product_kl
  entropy_minus_one <- transport + epsilon * (ot_entropy(plan) - sum(plan))

  active_source <- dat$p > 0
  active_target <- dat$q > 0
  source_kl_potential <- numeric(length(dat$p))
  target_kl_potential <- numeric(length(dat$q))
  source_kl_potential[active_source] <-
    out$source_potential[active_source] - epsilon * log(dat$p[active_source])
  target_kl_potential[active_target] <-
    out$target_potential[active_target] - epsilon * log(dat$q[active_target])

  log_reference <- outer(
    log(dat$p[active_source]), log(dat$q[active_target]), "+"
  )
  dual_exponent <- (
    outer(
      source_kl_potential[active_source],
      target_kl_potential[active_target],
      "+"
    ) - dat$M[active_source, active_target, drop = FALSE]
  ) / epsilon
  exponential_mass <- exp(.log_sum_exp(log_reference + dual_exponent))
  reference_mass <- sum(dat$p) * sum(dat$q)
  dual <- sum(dat$p * source_kl_potential) +
    sum(dat$q * target_kl_potential) -
    epsilon * (exponential_mass - reference_mass)
  gap <- primal - dual
  gap_tolerance <- max(
    1e-10,
    20 * tol * (1 + abs(primal) + abs(dual))
  )
  dual_consistent <- is.finite(gap) && gap >= -gap_tolerance &&
    abs(gap) <= gap_tolerance

  out$entropy_convention <- "epsilon_kl_plan_to_product_measure"
  out$entropy_reference <- "source_weights_tensor_target_weights"
  out$transport_objective <- transport
  out$product_measure_kl <- product_kl
  out$regularized_objective <- primal
  out$entropy_minus_one_objective <- entropy_minus_one
  out$product_reference_offset <- primal - entropy_minus_one
  out$regularized_source_potential <- source_kl_potential
  out$regularized_target_potential <- target_kl_potential
  out$regularized_dual_exponential_mass <- exponential_mass
  out$regularized_dual_objective <- dual
  out$regularized_duality_gap <- gap
  out$regularized_duality_gap_tolerance <- gap_tolerance
  out$regularized_dual_consistent <- dual_consistent
  out$regularized_objective_components <- list(
    transport = transport,
    product_measure_kl = product_kl,
    weighted_product_measure_kl = epsilon * product_kl,
    total = primal
  )
  out$regularized_certificate <- list(
    convention = out$entropy_convention,
    primal = primal,
    dual = dual,
    gap = gap,
    tolerance = gap_tolerance,
    consistent = dual_consistent
  )
  if (isTRUE(out$converged) && !dual_consistent) {
    out$converged <- FALSE
    out$status <- "objective_mismatch"
    out$termination_reason <- "regularized_primal_dual_gap"
    out$warning_payload <- list(
      code = "regularized_primal_dual_gap",
      message = sprintf(
        "Regularized primal-dual gap %s exceeds tolerance %s.",
        format(gap, digits = 6), format(gap_tolerance, digits = 6)
      )
    )
  }
  out
}

.validate_balanced_init_duals <- function(init_duals, p, q) {
  if (is.null(init_duals)) return(NULL)
  if (inherits(init_duals, "rfugw_result")) {
    init_duals <- init_duals$dual_state %||% list(
      source = init_duals$source_potential,
      target = init_duals$target_potential
    )
  }
  if (!is.list(init_duals) || is.null(init_duals$source) ||
      is.null(init_duals$target)) {
    stop(
      "`init_duals` must be a result or a list with `source` and `target`.",
      call. = FALSE
    )
  }
  source <- init_duals$source
  target <- init_duals$target
  if (!is.numeric(source) || length(source) != length(p)) {
    stop(
      sprintf("`init_duals$source` must be numeric with length %d.", length(p)),
      call. = FALSE
    )
  }
  if (!is.numeric(target) || length(target) != length(q)) {
    stop(
      sprintf("`init_duals$target` must be numeric with length %d.", length(q)),
      call. = FALSE
    )
  }
  if (any(!is.finite(source)) || any(!is.finite(target))) {
    stop("`init_duals` potentials must be finite.", call. = FALSE)
  }
  if (!is.null(init_duals$source_support) &&
      !identical(as.logical(init_duals$source_support), p > 0)) {
    stop(
      "`init_duals` source support is incompatible with the requested weights.",
      call. = FALSE
    )
  }
  if (!is.null(init_duals$target_support) &&
      !identical(as.logical(init_duals$target_support), q > 0)) {
    stop(
      "`init_duals` target support is incompatible with the requested weights.",
      call. = FALSE
    )
  }
  .canonicalize_balanced_duals(
    unname(as.numeric(source)), unname(as.numeric(target)), p
  )
}

.balanced_duals_from_plan <- function(init_plan, M, p, q, epsilon) {
  if (is.null(init_plan)) return(NULL)
  init_plan <- .validate_optional_init_plan(
    init_plan, length(p), length(q), "init_plan"
  )
  if (!length(init_plan)) return(NULL)
  active <- outer(p > 0, q > 0, "&")
  if (any(init_plan[!active] > 0)) {
    stop(
      "`init_plan` must be zero outside the positive source/target support.",
      call. = FALSE
    )
  }
  if (any(init_plan[active] <= 0)) {
    stop(
      "`init_plan` must be strictly positive on the active entropic support.",
      call. = FALSE
    )
  }

  active_i <- which(p > 0)
  active_j <- which(q > 0)
  pa <- p[active_i] / sum(p[active_i])
  qa <- q[active_j] / sum(q[active_j])
  log_factor <- epsilon * log(init_plan[active_i, active_j, drop = FALSE]) +
    M[active_i, active_j, drop = FALSE]
  source_active <- as.numeric(log_factor %*% qa)
  grand <- sum(pa * source_active)
  source_active <- source_active - grand
  target_active <- as.numeric(crossprod(pa, log_factor))

  source <- numeric(length(p))
  target <- numeric(length(q))
  source[active_i] <- source_active
  target[active_j] <- target_active
  .canonicalize_balanced_duals(source, target, p)
}

#' Balanced entropic optimal transport
#'
#' Scaling or log-domain Sinkhorn for a linear cost. Reported `ot_dist` is
#' the unregularized `<M, plan>` cost. `regularized_objective` uses the explicit
#' convention `<M, plan> + epsilon * KL(plan || p %o% q)`; its independently
#' derived dual and gap are returned without changing the legacy field.
#'
#' @param M Cost matrix (`ns x nt`).
#' @param p Source weights (default uniform). Renormalized to sum 1.
#' @param q Target weights (default uniform). Renormalized to sum 1.
#' @param epsilon Positive entropic regularization.
#' @param method `"scaling"`, `"log"`, or `"auto"`. Auto uses scaling only
#'   when the maximum scaled exponent magnitude and span are at most 500;
#'   otherwise it selects the genuine log-domain implementation.
#' @param init_plan Optional strictly positive entropic plan warm start. It must
#'   have the requested shape and be zero outside any zero-weight support.
#' @param init_duals Optional reusable dual state, either a prior
#'   `rfugw_result` or a list with numeric `source` and `target` potentials.
#'   When both initialization forms are supplied, `init_duals` takes precedence
#'   after both inputs are validated.
#' @param max_iter Maximum Sinkhorn iterations.
#' @param tol Marginal residual tolerance.
#' @return An `rfugw_result` with `plan`, `ot_dist`, `status`, residuals,
#'   canonical `source_potential` / `target_potential`, and reusable
#'   `dual_state`. The gauge is fixed by a zero source-weighted mean. The
#'   regularized primal, dual, product-reference KL, constant offset, and gap
#'   certificate are also returned.
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
    tol = 1e-9,
    init_plan = NULL,
    init_duals = NULL) {
  requested_method <- match.arg(method)
  dat <- .prepare_linear_ot(M, p, q)
  epsilon <- .validate_positive_scalar(epsilon, "epsilon")
  max_iter <- .validate_count(max_iter, "max_iter")
  tol <- .validate_positive_scalar(tol, "tol")
  plan_duals <- .balanced_duals_from_plan(
    init_plan, dat$M, dat$p, dat$q, epsilon
  )
  explicit_duals <- .validate_balanced_init_duals(init_duals, dat$p, dat$q)
  initial_duals <- explicit_duals %||% plan_duals
  dispatch <- .select_sinkhorn_method(
    requested_method, dat$M, epsilon, precision = "double",
    context = "Balanced Sinkhorn"
  )
  method <- dispatch$effective
  initialization <- if (!is.null(explicit_duals) && !is.null(plan_duals)) {
    "duals_over_plan"
  } else if (!is.null(explicit_duals)) {
    "duals"
  } else if (!is.null(plan_duals)) {
    "plan"
  } else {
    "cold"
  }
  out <- cpp_ot_sinkhorn(
    M = dat$M,
    p = dat$p,
    q = dat$q,
    epsilon = epsilon,
    max_iter = max_iter,
    tol = tol,
    use_log = identical(method, "log"),
    init_source_potential = initial_duals$source %||% numeric(),
    init_target_potential = initial_duals$target %||% numeric()
  )
  active_support <- outer(dat$p > 0, dat$q > 0, "&")
  out$plan[!active_support] <- 0
  out$ot_dist <- ot_linear_cost(dat$M, out$plan)
  out$source_potential <- as.numeric(out$source_potential)
  out$target_potential <- as.numeric(out$target_potential)
  out$formulation <- "ot_sinkhorn"
  out$regularization <- epsilon
  out$backend <- if (identical(method, "log")) "cpp_log" else "cpp_scaling"
  out$requested_sinkhorn_method <- dispatch$requested
  out$effective_sinkhorn_method <- dispatch$effective
  out$sinkhorn_backend_transition <- dispatch$transition
  out$sinkhorn_dispatch_reason <- dispatch$reason
  out$sinkhorn_dynamic_range <- dispatch$metric
  out$sinkhorn_scaling_threshold <- dispatch$threshold
  out$initialization <- initialization
  out$initialization_precedence <- "init_duals_over_init_plan"
  out$dual_state <- list(
    source = out$source_potential,
    target = out$target_potential,
    gauge = out$potential_gauge,
    epsilon = epsilon,
    method = method,
    source_support = dat$p > 0,
    target_support = dat$q > 0
  )
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
  .attach_balanced_regularized_contract(ans, dat, epsilon, tol)
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

.as_wasserstein_support <- function(x, name) {
  if (is.numeric(x) && is.null(dim(x))) {
    x <- matrix(x, ncol = 1L)
  }
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) < 1L || ncol(x) < 1L) {
    stop(
      sprintf("`%s` must be a nonempty numeric vector or matrix.", name),
      call. = FALSE
    )
  }
  if (any(!is.finite(x))) {
    stop(sprintf("`%s` must be finite.", name), call. = FALSE)
  }
  unname(x)
}

.wasserstein_support_cost <- function(source, target, metric, power) {
  if (ncol(source) != ncol(target)) {
    stop("`source` and `target` must have the same feature dimension.", call. = FALSE)
  }
  out <- matrix(0, nrow(source), nrow(target))
  for (i in seq_len(nrow(source))) {
    delta <- sweep(target, 2L, source[i, ], "-")
    out[i, ] <- if (identical(metric, "euclidean")) {
      sqrt(rowSums(delta^2))
    } else {
      rowSums(abs(delta))
    }
  }
  out^power
}

.prepare_wasserstein_problem <- function(
    source,
    target,
    cost,
    p,
    cost_power,
    metric,
    source_weights,
    target_weights) {
  p <- .validate_positive_scalar(p, "p")
  if (is.null(cost)) {
    if (is.null(source) || is.null(target)) {
      stop(
        "Supply both `source` and `target`, or supply `cost` with `cost_power`.",
        call. = FALSE
      )
    }
    if (!is.null(cost_power)) {
      stop("`cost_power` is only used with a supplied `cost` matrix.", call. = FALSE)
    }
    source <- .as_wasserstein_support(source, "source")
    target <- .as_wasserstein_support(target, "target")
    metric <- match.arg(metric, c("euclidean", "manhattan"))
    input_cost <- .wasserstein_support_cost(source, target, metric, 1)
    effective_cost <- input_cost^p
    input_cost_power <- 1
    input_mode <- "raw_supports"
    metric_certified <- TRUE
  } else {
    if (!is.null(source) || !is.null(target)) {
      stop("Supply either raw supports or `cost`, not both.", call. = FALSE)
    }
    if (is.null(cost_power)) {
      stop(
        "A supplied `cost` matrix requires an explicit positive `cost_power`.",
        call. = FALSE
      )
    }
    input_cost_power <- .validate_positive_scalar(cost_power, "cost_power")
    input_cost <- .validate_finite_matrix(cost, "cost")
    if (any(input_cost < 0)) {
      stop("`cost` must be nonnegative.", call. = FALSE)
    }
    effective_cost <- input_cost^(p / input_cost_power)
    metric <- "user_supplied"
    input_mode <- "cost_matrix"
    metric_certified <- FALSE
  }

  ns <- nrow(effective_cost)
  nt <- ncol(effective_cost)
  original_source <- .validate_finite_measure(
    source_weights, ns, "source_weights", rep(1 / ns, ns)
  )
  original_target <- .validate_finite_measure(
    target_weights, nt, "target_weights", rep(1 / nt, nt)
  )
  source_mass <- sum(original_source)
  target_mass <- sum(original_target)
  if (source_mass <= 0 || target_mass <= 0) {
    stop(
      "Classical Wasserstein helpers require positive source and target mass.",
      call. = FALSE
    )
  }

  list(
    cost = effective_cost,
    source_weights = original_source / source_mass,
    target_weights = original_target / target_mass,
    source_weights_original = original_source,
    target_weights_original = original_target,
    source_mass_original = source_mass,
    target_mass_original = target_mass,
    p = p,
    input_cost_power = input_cost_power,
    input_cost_scale = max(input_cost),
    effective_cost_scale = max(effective_cost),
    input_mode = input_mode,
    metric = metric,
    metric_certified = metric_certified
  )
}

.ot_wasserstein <- function(
    value,
    source,
    target,
    p,
    source_weights,
    target_weights,
    cost,
    cost_power,
    metric,
    solver,
    epsilon,
    sinkhorn_method,
    max_iter,
    tol,
    init_plan,
    init_duals,
    allow_uncertified) {
  solver <- match.arg(solver, c("exact", "sinkhorn"))
  if (length(allow_uncertified) != 1L || is.na(allow_uncertified) ||
      !is.logical(allow_uncertified)) {
    stop("`allow_uncertified` must be TRUE or FALSE.", call. = FALSE)
  }
  dat <- .prepare_wasserstein_problem(
    source = source,
    target = target,
    cost = cost,
    p = p,
    cost_power = cost_power,
    metric = metric,
    source_weights = source_weights,
    target_weights = target_weights
  )

  if (identical(solver, "exact")) {
    if (!is.null(epsilon)) {
      stop("`epsilon` is only operational for `solver = \"sinkhorn\"`.", call. = FALSE)
    }
    if (!is.null(sinkhorn_method)) {
      stop(
        "`sinkhorn_method` is only operational for `solver = \"sinkhorn\"`.",
        call. = FALSE
      )
    }
    if (!is.null(init_plan) || !is.null(init_duals)) {
      stop(
        "Warm starts are only supported for `solver = \"sinkhorn\"`.",
        call. = FALSE
      )
    }
    max_iter <- max_iter %||% 20000L
    tol <- tol %||% 1e-12
    out <- ot_emd(
      dat$cost, dat$source_weights, dat$target_weights,
      max_iter = max_iter, tol = tol
    )
    estimate_kind <- "exact"
  } else {
    epsilon <- epsilon %||% 0.05
    sinkhorn_method <- sinkhorn_method %||% "auto"
    max_iter <- max_iter %||% 1000L
    tol <- tol %||% 1e-9
    out <- ot_sinkhorn(
      dat$cost, dat$source_weights, dat$target_weights,
      epsilon = epsilon, method = sinkhorn_method,
      max_iter = max_iter, tol = tol,
      init_plan = init_plan, init_duals = init_duals
    )
    estimate_kind <- "entropic_plan"
  }

  if (!isTRUE(out$converged) && !allow_uncertified) {
    stop(
      sprintf(
        paste0(
          "Wasserstein %s is uncertified because the %s solve ended with ",
          "status `%s` (residual %s). Set `allow_uncertified = TRUE` only ",
          "to inspect the labeled plan value."
        ),
        value, solver, out$status %||% "unknown",
        format(out$residual %||% out$error, digits = 6)
      ),
      call. = FALSE
    )
  }

  base_formulation <- out$formulation
  p_cost <- max(0, as.numeric(out$ot_dist))
  p_distance <- p_cost^(1 / dat$p)
  out$formulation <- paste0("wasserstein_p_", value)
  out$base_formulation <- base_formulation
  out$wasserstein_value <- if (identical(value, "cost")) p_cost else p_distance
  out$wasserstein_p_cost <- p_cost
  out$wasserstein_p_distance <- p_distance
  out$wasserstein_power <- dat$p
  out$value_root <- if (identical(value, "cost")) 1 else 1 / dat$p
  out$estimate_kind <- estimate_kind
  out$value_kind <- paste(estimate_kind, paste0("p_", value), sep = "_")
  out$value_certified <- isTRUE(out$converged)
  out$value_certification <- if (isTRUE(out$converged) && identical(solver, "exact")) {
    "certified_exact_transport"
  } else if (isTRUE(out$converged)) {
    "certified_feasible_entropic_plan_value_not_exact_wasserstein"
  } else {
    "uncertified_plan_value_requested_explicitly"
  }
  out$input_mode <- dat$input_mode
  out$metric <- dat$metric
  out$metric_certified <- dat$metric_certified
  out$input_cost_power <- dat$input_cost_power
  out$effective_cost_power <- dat$p
  out$input_cost_scale <- dat$input_cost_scale
  out$effective_cost_scale <- dat$effective_cost_scale
  out$measure_normalization <- "separate_probability"
  out$source_measure_original <- dat$source_weights_original
  out$target_measure_original <- dat$target_weights_original
  out$source_measure_effective <- dat$source_weights
  out$target_measure_effective <- dat$target_weights
  out$original_source_mass <- dat$source_mass_original
  out$original_target_mass <- dat$target_mass_original
  out$effective_source_mass <- 1
  out$effective_target_mass <- 1
  out
}

#' Wasserstein p-cost with explicit cost-power semantics
#'
#' Solves classical balanced transport with ground cost `d^p` and returns a
#' certificate-rich result. Supply either raw supports or a nonnegative cost
#' matrix whose existing power is declared explicitly. A Sinkhorn result is
#' labeled as the cost of its certified entropic plan, not as exact
#' Wasserstein cost and not as a regularized objective.
#'
#' @param source Source support as a numeric vector or rows-by-features matrix.
#' @param target Target support with the same feature dimension.
#' @param p Positive Wasserstein ground-metric power.
#' @param source_weights Nonnegative source weights; normalized to probability.
#' @param target_weights Nonnegative target weights; normalized to probability.
#' @param cost Optional nonnegative source-by-target cost matrix. When supplied,
#'   `source` and `target` must be `NULL` and `cost_power` is required.
#' @param cost_power Positive power already represented by `cost`. The effective
#'   solver cost is `cost^(p / cost_power)`.
#' @param metric Raw-support metric, `"euclidean"` or `"manhattan"`.
#' @param solver `"exact"` for certified EMD or `"sinkhorn"` for the distinctly
#'   labeled cost of a certified entropic plan.
#' @param epsilon Entropic regularization for the Sinkhorn solver only.
#' @param sinkhorn_method `"auto"`, `"log"`, or `"scaling"` for Sinkhorn only.
#' @param max_iter Optional solver iteration limit; backend defaults apply.
#' @param tol Optional solver tolerance; backend defaults apply.
#' @param init_plan Optional Sinkhorn plan warm start.
#' @param init_duals Optional Sinkhorn dual warm state.
#' @param allow_uncertified If `FALSE` (default), fail instead of returning a
#'   value from a nonconverged solve. If `TRUE`, the result remains explicitly
#'   labeled uncertified.
#' @return An `rfugw_result`. `wasserstein_p_cost` is the transport cost,
#'   `wasserstein_p_distance` is its `1 / p` root, and `rfugw_value()` returns
#'   the p-cost for this helper.
#' @examples
#' out <- ot_wasserstein_cost(c(0, 2), c(1, 3), p = 2)
#' rfugw_value(out)
#' out$wasserstein_p_distance
#' @export
ot_wasserstein_cost <- function(
    source = NULL,
    target = NULL,
    p = 2,
    source_weights = NULL,
    target_weights = NULL,
    cost = NULL,
    cost_power = NULL,
    metric = c("euclidean", "manhattan"),
    solver = c("exact", "sinkhorn"),
    epsilon = NULL,
    sinkhorn_method = NULL,
    max_iter = NULL,
    tol = NULL,
    init_plan = NULL,
    init_duals = NULL,
    allow_uncertified = FALSE) {
  .ot_wasserstein(
    value = "cost", source = source, target = target, p = p,
    source_weights = source_weights, target_weights = target_weights,
    cost = cost, cost_power = cost_power, metric = metric, solver = solver,
    epsilon = epsilon, sinkhorn_method = sinkhorn_method,
    max_iter = max_iter, tol = tol, init_plan = init_plan,
    init_duals = init_duals, allow_uncertified = allow_uncertified
  )
}

#' Wasserstein p-distance with explicit root semantics
#'
#' Calls the same certified transport path as [ot_wasserstein_cost()] but makes
#' the returned primary value `(minimum d^p cost)^(1/p)`. For `p = 1` this is
#' exactly the transport cost; no square root is applied.
#'
#' @inheritParams ot_wasserstein_cost
#' @return An `rfugw_result` whose `rfugw_value()` is
#'   `wasserstein_p_distance`. The unrooted value remains available as
#'   `wasserstein_p_cost`.
#' @examples
#' out <- ot_wasserstein_distance(c(0, 2), c(1, 3), p = 1)
#' rfugw_value(out)
#' @export
ot_wasserstein_distance <- function(
    source = NULL,
    target = NULL,
    p = 2,
    source_weights = NULL,
    target_weights = NULL,
    cost = NULL,
    cost_power = NULL,
    metric = c("euclidean", "manhattan"),
    solver = c("exact", "sinkhorn"),
    epsilon = NULL,
    sinkhorn_method = NULL,
    max_iter = NULL,
    tol = NULL,
    init_plan = NULL,
    init_duals = NULL,
    allow_uncertified = FALSE) {
  .ot_wasserstein(
    value = "distance", source = source, target = target, p = p,
    source_weights = source_weights, target_weights = target_weights,
    cost = cost, cost_power = cost_power, metric = metric, solver = solver,
    epsilon = epsilon, sinkhorn_method = sinkhorn_method,
    max_iter = max_iter, tol = tol, init_plan = init_plan,
    init_duals = init_duals, allow_uncertified = allow_uncertified
  )
}

.prepare_sinkhorn_divergence_problem <- function(
    source,
    target,
    p,
    source_weights,
    target_weights,
    cost,
    source_cost,
    target_cost,
    cost_power,
    metric) {
  if (is.null(cost)) {
    if (!is.null(source_cost) || !is.null(target_cost)) {
      stop(
        "`source_cost` and `target_cost` require a supplied cross `cost`.",
        call. = FALSE
      )
    }
    source <- .as_wasserstein_support(source, "source")
    target <- .as_wasserstein_support(target, "target")
    metric <- match.arg(metric, c("euclidean", "manhattan"))
    dat <- .prepare_wasserstein_problem(
      source, target, NULL, p, NULL, metric,
      source_weights, target_weights
    )
    dat$source_cost <- .wasserstein_support_cost(source, source, metric, dat$p)
    dat$target_cost <- .wasserstein_support_cost(target, target, metric, dat$p)
    return(dat)
  }

  if (is.null(source_cost) || is.null(target_cost)) {
    stop(
      "A supplied cross `cost` requires both `source_cost` and `target_cost`.",
      call. = FALSE
    )
  }
  dat <- .prepare_wasserstein_problem(
    NULL, NULL, cost, p, cost_power, metric,
    source_weights, target_weights
  )
  source_cost <- .validate_finite_matrix(source_cost, "source_cost", square = TRUE)
  target_cost <- .validate_finite_matrix(target_cost, "target_cost", square = TRUE)
  if (nrow(source_cost) != nrow(dat$cost)) {
    stop("`source_cost` must match the source dimension of `cost`.", call. = FALSE)
  }
  if (nrow(target_cost) != ncol(dat$cost)) {
    stop("`target_cost` must match the target dimension of `cost`.", call. = FALSE)
  }
  if (any(source_cost < 0) || any(target_cost < 0)) {
    stop("Self-cost matrices must be nonnegative.", call. = FALSE)
  }
  if (!.is_symmetric_cost(source_cost) || !.is_symmetric_cost(target_cost)) {
    stop("Self-cost matrices must satisfy the symmetric-cost contract.", call. = FALSE)
  }
  diagonal_tolerance <- 1e-12 * (1 + max(source_cost, target_cost))
  if (max(abs(diag(source_cost)), abs(diag(target_cost))) > diagonal_tolerance) {
    stop("Self-cost matrices must have a zero diagonal.", call. = FALSE)
  }
  dat$source_cost <- source_cost^(dat$p / dat$input_cost_power)
  dat$target_cost <- target_cost^(dat$p / dat$input_cost_power)
  dat
}

.timed_balanced_sinkhorn <- function(...) {
  started <- proc.time()[["elapsed"]]
  result <- ot_sinkhorn(...)
  list(result = result, seconds = proc.time()[["elapsed"]] - started)
}

#' Debiased Sinkhorn divergence with certified component solves
#'
#' Computes the debiased quantity
#' `OTe(source, target) - 0.5 * OTe(source, source) -
#' 0.5 * OTe(target, target)`, where every `OTe` uses the same convention
#' `<cost, plan> + epsilon * KL(plan || source_weights %o% target_weights)`.
#' All three component solves must earn fresh feasibility, objective, and
#' regularized primal-dual certificates.
#'
#' @inheritParams ot_wasserstein_cost
#' @param source_cost Optional source self-cost matrix. Required with `cost`.
#' @param target_cost Optional target self-cost matrix. Required with `cost`.
#' @param method Balanced Sinkhorn backend: `"auto"`, `"log"`, or `"scaling"`.
#' @param init_state Optional list with reusable `cross`, `source_self`, and/or
#'   `target_self` dual states or prior results.
#' @return An `rfugw_result` whose primary value is `sinkhorn_divergence`.
#'   `component_solves`, `component_status`, `component_values`,
#'   `component_residuals`, and `component_runtime_seconds` retain the complete
#'   three-solve provenance.
#' @examples
#' x <- c(0, 1, 3)
#' y <- c(0.5, 2, 4)
#' out <- ot_sinkhorn_divergence(x, y, p = 2, epsilon = 1)
#' rfugw_value(out)
#' @export
ot_sinkhorn_divergence <- function(
    source = NULL,
    target = NULL,
    p = 2,
    source_weights = NULL,
    target_weights = NULL,
    cost = NULL,
    source_cost = NULL,
    target_cost = NULL,
    cost_power = NULL,
    metric = c("euclidean", "manhattan"),
    epsilon = 0.05,
    method = c("auto", "log", "scaling"),
    max_iter = 1000L,
    tol = 1e-9,
    init_state = NULL,
    allow_uncertified = FALSE) {
  epsilon <- .validate_positive_scalar(epsilon, "epsilon")
  method <- match.arg(method)
  max_iter <- .validate_count(max_iter, "max_iter")
  tol <- .validate_positive_scalar(tol, "tol")
  if (length(allow_uncertified) != 1L || is.na(allow_uncertified) ||
      !is.logical(allow_uncertified)) {
    stop("`allow_uncertified` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(init_state) && !is.list(init_state)) {
    stop("`init_state` must be a list of component dual states.", call. = FALSE)
  }
  dat <- .prepare_sinkhorn_divergence_problem(
    source = source, target = target, p = p,
    source_weights = source_weights, target_weights = target_weights,
    cost = cost, source_cost = source_cost, target_cost = target_cost,
    cost_power = cost_power, metric = metric
  )

  cross <- .timed_balanced_sinkhorn(
    dat$cost, dat$source_weights, dat$target_weights,
    epsilon = epsilon, method = method, max_iter = max_iter, tol = tol,
    init_duals = init_state$cross %||% NULL
  )
  source_self <- .timed_balanced_sinkhorn(
    dat$source_cost, dat$source_weights, dat$source_weights,
    epsilon = epsilon, method = method, max_iter = max_iter, tol = tol,
    init_duals = init_state$source_self %||% NULL
  )
  target_self <- .timed_balanced_sinkhorn(
    dat$target_cost, dat$target_weights, dat$target_weights,
    epsilon = epsilon, method = method, max_iter = max_iter, tol = tol,
    init_duals = init_state$target_self %||% NULL
  )
  components <- list(
    cross = cross$result,
    source_self = source_self$result,
    target_self = target_self$result
  )
  component_status <- vapply(components, function(x) x$status, character(1))
  component_converged <- vapply(
    components,
    function(x) isTRUE(x$converged) && isTRUE(x$regularized_dual_consistent),
    logical(1)
  )
  component_values <- vapply(
    components, function(x) x$regularized_objective, numeric(1)
  )
  component_residuals <- vapply(
    components, function(x) x$residual %||% x$error, numeric(1)
  )
  component_runtime <- c(
    cross = cross$seconds,
    source_self = source_self$seconds,
    target_self = target_self$seconds
  )
  divergence <- component_values[["cross"]] -
    0.5 * component_values[["source_self"]] -
    0.5 * component_values[["target_self"]]
  nonnegative_tolerance <- max(
    1e-10,
    20 * tol * (1 + sum(abs(component_values)))
  )
  nonnegative <- is.finite(divergence) && divergence >= -nonnegative_tolerance
  certified <- all(component_converged) && nonnegative

  if (!all(component_converged) && !allow_uncertified) {
    failed <- names(component_converged)[!component_converged]
    stop(
      sprintf(
        "Sinkhorn divergence is uncertified because component solve(s) failed: %s.",
        paste(sprintf("%s=%s", failed, component_status[failed]), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (all(component_converged) && !nonnegative && !allow_uncertified) {
    stop(
      sprintf(
        "Sinkhorn divergence %s is below the nonnegativity tolerance %s.",
        format(divergence, digits = 6), format(nonnegative_tolerance, digits = 6)
      ),
      call. = FALSE
    )
  }

  out <- components$cross
  out$formulation <- "sinkhorn_divergence"
  out$base_formulation <- "ot_sinkhorn"
  out$sinkhorn_divergence <- divergence
  out$wasserstein_power <- dat$p
  out$value_kind <- "debiased_regularized_ot"
  out$value_certified <- certified
  out$value_certification <- if (certified) {
    "three_certified_regularized_ot_solves"
  } else {
    "uncertified_component_or_nonnegativity_failure"
  }
  out$entropy_convention <- "epsilon_kl_plan_to_product_measure"
  out$component_solves <- components
  out$component_status <- component_status
  out$component_converged <- component_converged
  out$component_values <- component_values
  out$component_residuals <- component_residuals
  out$component_iterations <- vapply(
    components, function(x) x$iterations, integer(1)
  )
  out$component_runtime_seconds <- component_runtime
  out$runtime_seconds <- sum(component_runtime)
  out$iterations <- sum(out$component_iterations)
  out$max_iter <- 3L * max_iter
  out$residual <- max(component_residuals)
  out$nonnegative_tolerance <- nonnegative_tolerance
  out$nonnegative_within_tolerance <- nonnegative
  out$input_mode <- dat$input_mode
  out$metric <- dat$metric
  out$metric_certified <- dat$metric_certified
  out$input_cost_power <- dat$input_cost_power
  out$effective_cost_power <- dat$p
  out$input_cost_scale <- dat$input_cost_scale
  out$effective_cost_scale <- dat$effective_cost_scale
  out$measure_normalization <- "separate_probability"
  out$source_measure_original <- dat$source_weights_original
  out$target_measure_original <- dat$target_weights_original
  out$source_measure_effective <- dat$source_weights
  out$target_measure_effective <- dat$target_weights
  out$original_source_mass <- dat$source_mass_original
  out$original_target_mass <- dat$target_mass_original
  out$effective_source_mass <- 1
  out$effective_target_mass <- 1
  out$converged <- certified
  out$status <- if (certified) {
    "converged"
  } else if (!all(component_converged)) {
    "component_failure"
  } else {
    "objective_mismatch"
  }
  out$termination_reason <- if (certified) {
    "all_components_certified"
  } else if (!all(component_converged)) {
    "uncertified_component_solve"
  } else {
    "negative_beyond_tolerance"
  }
  out$warning_payload <- if (certified) NULL else list(
    code = out$termination_reason,
    message = "Sinkhorn divergence was returned only because `allow_uncertified = TRUE`."
  )
  out
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

.partial_sinkhorn_active_problem <- function(M, p, q, mass, epsilon) {
  active_source <- which(p > 0)
  active_target <- which(q > 0)
  active_cost <- M[active_source, active_target, drop = FALSE]
  shifted_cost <- active_cost - min(active_cost)
  scaled_cost <- shifted_cost / epsilon
  if (any(!is.finite(scaled_cost))) {
    stop(
      "The cost range divided by `epsilon` exceeds finite double precision.",
      call. = FALSE
    )
  }
  log_reference <- -scaled_cost
  log_reference <- log_reference + log(mass) -
    .log_sum_exp(as.numeric(log_reference))
  list(
    active_source = active_source,
    active_target = active_target,
    cost = active_cost,
    p = p[active_source],
    q = q[active_target],
    mass = mass,
    epsilon = epsilon,
    log_reference = log_reference
  )
}

.linear_partial_sinkhorn_dispatch <- function(method, M, epsilon) {
  method <- match.arg(method, c("auto", "scaling", "log"))
  dynamic_range <- diff(range(M)) / epsilon
  threshold <- 100
  unsafe <- !is.finite(dynamic_range) || dynamic_range > threshold
  if (identical(method, "scaling") && unsafe) {
    stop(
      sprintf(
        paste0(
          "Entropic partial OT scaling is outside its certified regime: ",
          "shift-invariant dynamic range %.6g exceeds %.6g. Use ",
          "`method = \"log\"` or increase `epsilon`."
        ),
        dynamic_range, threshold
      ),
      call. = FALSE
    )
  }
  effective <- if (identical(method, "auto")) {
    if (unsafe) "log" else "scaling"
  } else {
    method
  }
  list(
    requested = method,
    effective = effective,
    metric = dynamic_range,
    threshold = threshold,
    reason = if (!identical(method, "auto")) {
      "explicit_request"
    } else if (unsafe) {
      "dynamic_range_exceeds_scaling_threshold"
    } else {
      "dynamic_range_within_scaling_threshold"
    },
    transition = if (identical(method, effective)) {
      "none"
    } else {
      paste0("auto_to_", effective)
    }
  )
}

.partial_sinkhorn_state_matches <- function(x, y) {
  isTRUE(all.equal(x, y, tolerance = 0, check.attributes = FALSE))
}

.validate_partial_sinkhorn_state <- function(init_state, problem, method) {
  if (is.null(init_state)) return(NULL)
  if (inherits(init_state, "rfugw_result")) {
    if (!isTRUE(init_state$converged)) {
      stop("`init_state` result must be certified and converged.", call. = FALSE)
    }
    init_state <- init_state$warm_state
  }
  if (!is.list(init_state) || !isTRUE(init_state$certified)) {
    stop("`init_state` must be a certified partial-Sinkhorn warm state.", call. = FALSE)
  }
  required <- c(
    "method", "active_source", "active_target", "cost", "p", "q",
    "mass", "epsilon", "plan_state", "corrections", "log_reference"
  )
  if (!all(required %in% names(init_state))) {
    stop("`init_state` is missing required partial-Sinkhorn fields.", call. = FALSE)
  }
  comparisons <- list(
    method = identical(init_state$method, method),
    active_source = identical(init_state$active_source, problem$active_source),
    active_target = identical(init_state$active_target, problem$active_target),
    cost = .partial_sinkhorn_state_matches(init_state$cost, problem$cost),
    p = .partial_sinkhorn_state_matches(init_state$p, problem$p),
    q = .partial_sinkhorn_state_matches(init_state$q, problem$q),
    mass = .partial_sinkhorn_state_matches(init_state$mass, problem$mass),
    epsilon = .partial_sinkhorn_state_matches(
      init_state$epsilon, problem$epsilon
    ),
    log_reference = .partial_sinkhorn_state_matches(
      init_state$log_reference, problem$log_reference
    )
  )
  failed <- names(comparisons)[!vapply(comparisons, isTRUE, logical(1))]
  if (length(failed)) {
    stop(
      sprintf(
        "`init_state` does not match the current problem: %s.",
        paste(failed, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  dims <- dim(problem$cost)
  if (!is.list(init_state$corrections) ||
      length(init_state$corrections) != 3L ||
      !identical(dim(init_state$plan_state), dims) ||
      any(!vapply(
        init_state$corrections,
        function(x) identical(dim(x), dims), logical(1)
      ))) {
    stop("`init_state` has incompatible state dimensions.", call. = FALSE)
  }
  values <- c(init_state$plan_state, unlist(init_state$corrections))
  if (any(!is.finite(values))) {
    stop("`init_state` must contain only finite active-support state.", call. = FALSE)
  }
  if (identical(method, "scaling") && any(values <= 0)) {
    stop("Scaling-domain `init_state` entries must be strictly positive.", call. = FALSE)
  }
  log_plan <- if (identical(method, "log")) {
    init_state$plan_state
  } else {
    log(init_state$plan_state)
  }
  log_corrections <- if (identical(method, "log")) {
    init_state$corrections
  } else {
    lapply(init_state$corrections, log)
  }
  invariant <- log_plan + log_corrections[[1L]] +
    log_corrections[[2L]] + log_corrections[[3L]] -
    problem$log_reference
  if (max(abs(invariant)) > 1e-7) {
    stop("`init_state` violates the Dykstra reference invariant.", call. = FALSE)
  }
  active_plan <- exp(log_plan)
  feasibility <- max(
    max(c(rowSums(active_plan) - problem$p, 0)),
    max(c(colSums(active_plan) - problem$q, 0)),
    abs(sum(active_plan) - problem$mass)
  )
  if (!is.finite(feasibility) || feasibility > 1e-8) {
    stop("`init_state` plan is infeasible for the current problem.", call. = FALSE)
  }
  init_state
}

.partial_sinkhorn_scaling_core <- function(
    problem, max_iter, tol, check_every, init_state, verbose, certify) {
  if (is.null(init_state)) {
    plan <- exp(problem$log_reference)
    corrections <- replicate(3L, matrix(1, nrow(plan), ncol(plan)), FALSE)
  } else {
    plan <- init_state$plan_state
    corrections <- init_state$corrections
  }
  trace <- data.frame(
    iteration = integer(), update_residual = numeric(),
    row_violation = numeric(), col_violation = numeric(),
    mass_residual = numeric()
  )
  update_residual <- Inf
  row_violation <- Inf
  col_violation <- Inf
  mass_residual <- Inf
  numerical_failure <- FALSE
  iterations <- 0L

  for (k in seq_len(max_iter)) {
    previous <- plan

    working <- plan * corrections[[1L]]
    row_scale <- pmin(problem$p / rowSums(working), 1)
    projected <- sweep(working, 1L, row_scale, `*`)
    corrections[[1L]] <- corrections[[1L]] * plan / projected
    plan <- projected

    previous_projection <- plan
    working <- plan * corrections[[2L]]
    col_scale <- pmin(problem$q / colSums(working), 1)
    projected <- sweep(working, 2L, col_scale, `*`)
    corrections[[2L]] <- corrections[[2L]] *
      previous_projection / projected
    plan <- projected

    previous_projection <- plan
    working <- plan * corrections[[3L]]
    projected <- working * (problem$mass / sum(working))
    corrections[[3L]] <- corrections[[3L]] *
      previous_projection / projected
    plan <- projected
    iterations <- k

    state_values <- c(plan, unlist(corrections))
    if (any(!is.finite(state_values)) || any(state_values <= 0)) {
      numerical_failure <- TRUE
      break
    }
    if (k == 1L || k %% check_every == 0L || k == max_iter) {
      update_residual <- max(abs(plan - previous))
      row_violation <- max(c(rowSums(plan) - problem$p, 0))
      col_violation <- max(c(colSums(plan) - problem$q, 0))
      mass_residual <- abs(sum(plan) - problem$mass)
      trace <- rbind(trace, data.frame(
        iteration = k,
        update_residual = update_residual,
        row_violation = row_violation,
        col_violation = col_violation,
        mass_residual = mass_residual
      ))
      if (isTRUE(verbose)) {
        cat(sprintf(
          "it=%d update=%.3e row=%.3e col=%.3e mass=%.3e\n",
          k, update_residual, row_violation, col_violation, mass_residual
        ))
      }
      if (max(
        update_residual, row_violation, col_violation, mass_residual
      ) <= tol && certify(list(
        plan = plan, log_plan = log(plan),
        log_corrections = lapply(corrections, log),
        row_violation = row_violation, col_violation = col_violation,
        mass_residual = mass_residual
      ))) {
        break
      }
    }
  }
  list(
    plan = plan,
    log_plan = log(plan),
    corrections = corrections,
    log_corrections = lapply(corrections, log),
    iterations = iterations,
    update_residual = update_residual,
    row_violation = row_violation,
    col_violation = col_violation,
    mass_residual = mass_residual,
    trace = trace,
    numerical_failure = numerical_failure
  )
}

.partial_sinkhorn_log_core <- function(
    problem, max_iter, tol, check_every, init_state, verbose, certify) {
  if (is.null(init_state)) {
    log_plan <- problem$log_reference
    log_corrections <- replicate(
      3L, matrix(0, nrow(log_plan), ncol(log_plan)), FALSE
    )
  } else {
    log_plan <- init_state$plan_state
    log_corrections <- init_state$corrections
  }
  trace <- data.frame(
    iteration = integer(), update_residual = numeric(),
    row_violation = numeric(), col_violation = numeric(),
    mass_residual = numeric()
  )
  update_residual <- Inf
  row_violation <- Inf
  col_violation <- Inf
  mass_residual <- Inf
  numerical_failure <- FALSE
  iterations <- 0L

  for (k in seq_len(max_iter)) {
    previous <- log_plan

    working <- log_plan + log_corrections[[1L]]
    row_totals <- apply(working, 1L, .log_sum_exp)
    row_scale <- pmin(log(problem$p) - row_totals, 0)
    projected <- sweep(working, 1L, row_scale, `+`)
    log_corrections[[1L]] <- log_corrections[[1L]] +
      log_plan - projected
    log_plan <- projected

    previous_projection <- log_plan
    working <- log_plan + log_corrections[[2L]]
    col_totals <- apply(working, 2L, .log_sum_exp)
    col_scale <- pmin(log(problem$q) - col_totals, 0)
    projected <- sweep(working, 2L, col_scale, `+`)
    log_corrections[[2L]] <- log_corrections[[2L]] +
      previous_projection - projected
    log_plan <- projected

    previous_projection <- log_plan
    working <- log_plan + log_corrections[[3L]]
    mass_scale <- log(problem$mass) - .log_sum_exp(as.numeric(working))
    projected <- working + mass_scale
    log_corrections[[3L]] <- log_corrections[[3L]] +
      previous_projection - projected
    log_plan <- projected
    iterations <- k

    state_values <- c(log_plan, unlist(log_corrections))
    if (any(!is.finite(state_values))) {
      numerical_failure <- TRUE
      break
    }
    if (k == 1L || k %% check_every == 0L || k == max_iter) {
      plan <- exp(log_plan)
      previous_plan <- exp(previous)
      update_residual <- max(abs(plan - previous_plan))
      row_violation <- max(c(rowSums(plan) - problem$p, 0))
      col_violation <- max(c(colSums(plan) - problem$q, 0))
      mass_residual <- abs(sum(plan) - problem$mass)
      trace <- rbind(trace, data.frame(
        iteration = k,
        update_residual = update_residual,
        row_violation = row_violation,
        col_violation = col_violation,
        mass_residual = mass_residual
      ))
      if (isTRUE(verbose)) {
        cat(sprintf(
          "it=%d update=%.3e row=%.3e col=%.3e mass=%.3e\n",
          k, update_residual, row_violation, col_violation, mass_residual
        ))
      }
      if (max(
        update_residual, row_violation, col_violation, mass_residual
      ) <= tol && certify(list(
        plan = plan, log_plan = log_plan,
        log_corrections = log_corrections,
        row_violation = row_violation, col_violation = col_violation,
        mass_residual = mass_residual
      ))) {
        break
      }
    }
  }
  list(
    plan = exp(log_plan),
    log_plan = log_plan,
    corrections = lapply(log_corrections, exp),
    log_corrections = log_corrections,
    iterations = iterations,
    update_residual = update_residual,
    row_violation = row_violation,
    col_violation = col_violation,
    mass_residual = mass_residual,
    trace = trace,
    numerical_failure = numerical_failure
  )
}

.partial_sinkhorn_certificate <- function(
    core, problem, M, p, q, epsilon, mass, tol) {
  active_plan <- core$plan
  plan <- matrix(0, nrow(M), ncol(M))
  plan[problem$active_source, problem$active_target] <- active_plan
  transport <- sum(M * plan)
  entropy <- sum(active_plan * core$log_plan)
  entropy_minus_one <- entropy - mass
  objective <- transport + epsilon * entropy_minus_one

  source_dual_active <- epsilon * rowMeans(core$log_corrections[[1L]])
  target_dual_active <- epsilon * colMeans(core$log_corrections[[2L]])
  stationarity_base <- problem$cost + epsilon * core$log_plan +
    outer(source_dual_active, target_dual_active, "+")
  mass_dual <- -mean(stationarity_base)
  stationarity_residual <- max(abs(stationarity_base + mass_dual))
  source_dual <- numeric(length(p))
  target_dual <- numeric(length(q))
  source_dual[problem$active_source] <- source_dual_active
  target_dual[problem$active_target] <- target_dual_active
  dual_objective <- -epsilon * mass - sum(source_dual * p) -
    sum(target_dual * q) - mass_dual * mass
  duality_gap <- objective - dual_objective
  source_slack <- p - rowSums(plan)
  target_slack <- q - colSums(plan)
  complementarity_residual <- max(
    abs(source_dual * source_slack),
    abs(target_dual * target_slack)
  )
  dual_feasibility_residual <- max(c(-source_dual, -target_dual, 0))
  feasibility_residual <- max(
    core$row_violation, core$col_violation, core$mass_residual,
    max(c(-plan, 0))
  )
  feasibility_tolerance <- max(1e-10, 10 * tol)
  certificate_tolerance <- max(
    1e-8,
    50 * tol * (1 + abs(objective) + abs(dual_objective) +
      max(abs(M)))
  )
  list(
    plan = plan,
    transport = transport,
    entropy = entropy,
    entropy_minus_one = entropy_minus_one,
    objective = objective,
    source_dual = source_dual,
    target_dual = target_dual,
    mass_dual = mass_dual,
    dual_objective = dual_objective,
    duality_gap = duality_gap,
    stationarity_residual = stationarity_residual,
    complementarity_residual = complementarity_residual,
    dual_feasibility_residual = dual_feasibility_residual,
    feasibility_residual = feasibility_residual,
    feasibility_tolerance = feasibility_tolerance,
    certificate_tolerance = certificate_tolerance,
    feasible = is.finite(feasibility_residual) &&
      feasibility_residual <= feasibility_tolerance,
    certified = all(is.finite(c(
      objective, dual_objective, duality_gap, stationarity_residual,
      complementarity_residual, dual_feasibility_residual
    ))) && duality_gap >= -certificate_tolerance &&
      abs(duality_gap) <= certificate_tolerance &&
      stationarity_residual <= certificate_tolerance &&
      complementarity_residual <= certificate_tolerance &&
      dual_feasibility_residual <= certificate_tolerance
  )
}

#' Entropy-Regularized Fixed-Mass Partial Optimal Transport
#'
#' Solves
#' `sum(M * G) + epsilon * sum(G * (log(G) - 1))` over nonnegative
#' subcouplings with `rowSums(G) <= p`, `colSums(G) <= q`, and
#' `sum(G) = mass`. Zero plan entries contribute zero to the entropy. This
#' counting-measure, entropy-minus-one convention is explicit and differs from
#' product-reference KL conventions.
#'
#' A genuine log-domain Dykstra backend is available for large dynamic range.
#' The scaling backend is accepted only inside its declared range; `"auto"`
#' dispatches safely. Returned convergence additionally requires primal
#' feasibility, independent objective recomputation, dual feasibility,
#' complementarity, stationarity, and a primal-dual gap certificate.
#'
#' @param M Finite source-by-target transport cost.
#' @param p,q Nonnegative source and target weights. Each positive-mass vector
#'   is normalized to a probability measure, consistently with
#'   [ot_partial_emd()].
#' @param mass Requested transported probability mass in `[0, 1]`.
#' @param epsilon Positive entropy coefficient.
#' @param method `"auto"`, `"scaling"`, or genuine `"log"` Dykstra.
#' @param max_iter Maximum Dykstra cycles.
#' @param tol Requested update and feasibility tolerance.
#' @param check_every Iteration interval for convergence diagnostics.
#' @param init_state Optional certified `warm_state` or prior converged result
#'   from the identical problem and effective backend. Raw plan starts are not
#'   accepted because they do not determine valid Dykstra corrections.
#' @param verbose If `TRUE`, print checked residuals.
#' @return An `rfugw_result` with the plan, exact entropy convention and
#'   objective decomposition, fixed-mass feasibility, KKT/duality certificate,
#'   dynamic-range dispatch, reusable warm state, and diagnostic trace.
#' @examples
#' M <- matrix(c(0, 2, 1, 0.3), 2, 2)
#' out <- ot_partial_sinkhorn(M, mass = 0.7, epsilon = 0.2)
#' out$status
#' sum(rfugw_plan(out))
#' @export
ot_partial_sinkhorn <- function(
    M,
    p = NULL,
    q = NULL,
    mass = 1,
    epsilon = 0.1,
    method = c("auto", "scaling", "log"),
    max_iter = 10000L,
    tol = 1e-9,
    check_every = 10L,
    init_state = NULL,
    verbose = FALSE) {
  dat <- .prepare_linear_ot(M, p, q)
  if (!is.numeric(mass) || length(mass) != 1L ||
      !is.finite(mass) || mass < 0 || mass > 1) {
    stop("`mass` must be one finite number in [0, 1].", call. = FALSE)
  }
  epsilon <- .validate_positive_scalar(epsilon, "epsilon")
  max_iter <- .validate_count(max_iter, "max_iter")
  tol <- .validate_positive_scalar(tol, "tol")
  check_every <- .validate_count(check_every, "check_every")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  method <- match.arg(method)

  if (mass == 0) {
    if (!is.null(init_state)) {
      stop("`init_state` is not used for the analytic zero-mass problem.", call. = FALSE)
    }
    plan <- matrix(0, nrow(dat$M), ncol(dat$M))
    out <- list(
      plan = plan,
      partial_sinkhorn_objective = 0,
      regularized_objective = 0,
      ot_dist = 0,
      transport_objective = 0,
      entropy = 0,
      entropy_minus_one = 0,
      weighted_entropy_minus_one = 0,
      entropy_convention = "counting_measure_entropy_minus_one",
      entropy_reference = "counting_measure",
      transported_mass = 0,
      mass = 0,
      mass_target = 0,
      mass_residual = 0,
      row_residual = 0,
      col_residual = 0,
      feasibility = "partial_fixed_mass",
      feasibility_residual = 0,
      feasibility_tolerance = tol,
      feasible = TRUE,
      objective_recomputed = 0,
      objective_residual = 0,
      objective_tolerance = tol,
      objective_consistent = TRUE,
      objective_components_consistent = TRUE,
      primal_objective = 0,
      dual_objective = 0,
      duality_gap = 0,
      duality_gap_tolerance = tol,
      stationarity_residual = 0,
      complementarity_residual = 0,
      dual_feasibility_residual = 0,
      iterations = 0L,
      max_iter = max_iter,
      residual = 0,
      error = 0,
      requested_sinkhorn_method = method,
      effective_sinkhorn_method = "analytic",
      sinkhorn_backend_transition = "analytic_zero_mass",
      sinkhorn_dispatch_reason = "zero_mass",
      sinkhorn_dynamic_range = 0,
      sinkhorn_scaling_threshold = 100,
      warm_started = FALSE,
      warm_state = NULL,
      initialization = "analytic_zero_mass",
      formulation = "ot_partial_sinkhorn",
      backend = "analytic_zero_mass",
      regularization = epsilon,
      status = "converged",
      converged = TRUE,
      termination_reason = "zero_mass",
      warning_payload = NULL
    )
    out$runtime_provenance <- .solver_runtime_provenance(out)
    class(out) <- "rfugw_result"
    return(out)
  }

  problem <- .partial_sinkhorn_active_problem(
    dat$M, dat$p, dat$q, mass, epsilon
  )
  dispatch <- .linear_partial_sinkhorn_dispatch(method, problem$cost, epsilon)
  state <- .validate_partial_sinkhorn_state(
    init_state, problem, dispatch$effective
  )
  # Small changes in tiny plan entries can precede KKT convergence. Do not
  # discard the remaining budget merely because primal updates look settled.
  certify <- function(core) {
    certificate <- .partial_sinkhorn_certificate(
      core, problem, dat$M, dat$p, dat$q, epsilon, mass, tol
    )
    isTRUE(certificate$feasible) && isTRUE(certificate$certified)
  }
  core <- if (identical(dispatch$effective, "log")) {
    .partial_sinkhorn_log_core(
      problem, max_iter, tol, check_every, state, verbose, certify
    )
  } else {
    .partial_sinkhorn_scaling_core(
      problem, max_iter, tol, check_every, state, verbose, certify
    )
  }
  certificate <- .partial_sinkhorn_certificate(
    core, problem, dat$M, dat$p, dat$q, epsilon, mass, tol
  )
  algorithm_residual <- max(
    core$update_residual, core$row_violation,
    core$col_violation, core$mass_residual
  )
  algorithm_converged <- !core$numerical_failure &&
    is.finite(algorithm_residual) && algorithm_residual <= tol
  objective_recomputed <- ot_linear_cost(dat$M, certificate$plan) +
    epsilon * (ot_entropy(certificate$plan) - sum(certificate$plan))
  objective_residual <- abs(certificate$objective - objective_recomputed)
  objective_tolerance <- max(
    1e-10,
    20 * tol * (1 + abs(certificate$objective) +
      abs(objective_recomputed))
  )
  objective_consistent <- is.finite(objective_residual) &&
    objective_residual <= objective_tolerance
  certified <- algorithm_converged && certificate$feasible &&
    certificate$certified && objective_consistent
  status <- if (core$numerical_failure) {
    "numerical_failure"
  } else if (!algorithm_converged) {
    "max_iter"
  } else if (!certificate$feasible) {
    "infeasible"
  } else if (!certificate$certified || !objective_consistent) {
    "objective_mismatch"
  } else {
    "converged"
  }
  termination_reason <- switch(
    status,
    converged = "dykstra_kkt_and_duality_certified",
    max_iter = "maximum_iterations",
    numerical_failure = "nonfinite_dykstra_state",
    infeasible = "partial_capacity_or_mass_failure",
    objective_mismatch = "kkt_or_objective_certificate_failure"
  )
  warm_state <- list(
    certified = certified,
    method = dispatch$effective,
    active_source = problem$active_source,
    active_target = problem$active_target,
    cost = problem$cost,
    p = problem$p,
    q = problem$q,
    mass = mass,
    epsilon = epsilon,
    plan_state = if (identical(dispatch$effective, "log")) {
      core$log_plan
    } else {
      core$plan
    },
    corrections = if (identical(dispatch$effective, "log")) {
      core$log_corrections
    } else {
      core$corrections
    },
    log_reference = problem$log_reference
  )
  out <- list(
    plan = certificate$plan,
    partial_sinkhorn_objective = certificate$objective,
    regularized_objective = certificate$objective,
    ot_dist = certificate$transport,
    transport_objective = certificate$transport,
    entropy = certificate$entropy,
    entropy_minus_one = certificate$entropy_minus_one,
    weighted_entropy_minus_one = epsilon * certificate$entropy_minus_one,
    entropy_convention = "counting_measure_entropy_minus_one",
    entropy_reference = "counting_measure",
    transported_mass = sum(certificate$plan),
    mass = sum(certificate$plan),
    mass_target = mass,
    mass_residual = abs(sum(certificate$plan) - mass),
    transported_mass_target = mass,
    transported_mass_defaulted = FALSE,
    row_residual = core$row_violation,
    col_residual = core$col_violation,
    feasibility = "partial_fixed_mass",
    feasibility_residual = certificate$feasibility_residual,
    feasibility_tolerance = certificate$feasibility_tolerance,
    feasible = certificate$feasible,
    objective_recomputed = objective_recomputed,
    objective_residual = objective_residual,
    objective_tolerance = objective_tolerance,
    objective_consistent = objective_consistent,
    objective_components_consistent = all(is.finite(c(
      certificate$transport, certificate$entropy,
      certificate$entropy_minus_one, certificate$objective
    ))) && abs(
      certificate$transport + epsilon * certificate$entropy_minus_one -
        certificate$objective
    ) <= objective_tolerance,
    primal_objective = certificate$objective,
    dual_objective = certificate$dual_objective,
    duality_gap = certificate$duality_gap,
    duality_gap_tolerance = certificate$certificate_tolerance,
    source_capacity_potential = certificate$source_dual,
    target_capacity_potential = certificate$target_dual,
    mass_potential = certificate$mass_dual,
    stationarity_residual = certificate$stationarity_residual,
    complementarity_residual = certificate$complementarity_residual,
    dual_feasibility_residual = certificate$dual_feasibility_residual,
    kkt_tolerance = certificate$certificate_tolerance,
    iterations = as.integer(core$iterations),
    max_iter = max_iter,
    residual = algorithm_residual,
    error = algorithm_residual,
    update_residual = core$update_residual,
    convergence_trace = core$trace,
    requested_sinkhorn_method = dispatch$requested,
    effective_sinkhorn_method = dispatch$effective,
    sinkhorn_backend_transition = dispatch$transition,
    sinkhorn_dispatch_reason = dispatch$reason,
    sinkhorn_dynamic_range = dispatch$metric,
    sinkhorn_scaling_threshold = dispatch$threshold,
    warm_started = !is.null(state),
    warm_state = warm_state,
    initialization = if (is.null(state)) "cold_reference_kernel" else "certified_state",
    source_measure_original = p %||% dat$p,
    target_measure_original = q %||% dat$q,
    source_measure_effective = dat$p,
    target_measure_effective = dat$q,
    original_source_mass = if (is.null(p)) 1 else sum(p),
    original_target_mass = if (is.null(q)) 1 else sum(q),
    effective_source_mass = 1,
    effective_target_mass = 1,
    measure_normalization = "separate_probability",
    formulation = "ot_partial_sinkhorn",
    backend = paste0("dykstra_", dispatch$effective),
    regularization = epsilon,
    status = status,
    converged = identical(status, "converged"),
    termination_reason = termination_reason,
    warning_payload = if (identical(status, "converged")) NULL else list(
      code = status,
      message = sprintf("Entropic partial OT is not certified: %s.", status)
    )
  )
  out$runtime_provenance <- .solver_runtime_provenance(out)
  class(out) <- unique(c("rfugw_result", class(out)))
  out
}

#' Exact penalized variable-mass partial optimal transport
#'
#' Minimizes
#' `<M,G> + discard_penalty * (sum(p) - sum(G)) +
#' discard_penalty * (sum(q) - sum(G))`
#' over nonnegative plans with `rowSums(G) <= p` and `colSums(G) <= q`.
#' Thus a larger discard penalty weakly favors transporting more mass. Unlike
#' [ot_partial_emd()], total transported mass is chosen by the optimization;
#' unlike KL-unbalanced OT, the marginal penalty is linear in discarded mass.
#'
#' The implementation adds one dummy source and target. Real-to-dummy and
#' dummy-to-real edges cost exactly `discard_penalty`, while the dummy-to-dummy
#' edge costs zero. No hidden big-M or cost-derived penalty is used.
#'
#' @param M Finite nonnegative source-by-target transport cost.
#' @param p Source finite nonnegative measure (default uniform probability).
#' @param q Target finite nonnegative measure (default uniform probability).
#' @param discard_penalty Finite nonnegative cost per discarded unit on each
#'   marginal. Transporting one unit avoids two discard penalties.
#' @param max_iter Maximum augmented simplex iterations.
#' @param tol Augmented exact-solver tolerance.
#' @return An `rfugw_result` with the optimized `plan`, transported and
#'   discarded masses, transport/penalty objective decomposition, and the
#'   scaled primal-dual/reduced-cost certificate for the equivalent augmented
#'   problem.
#' @examples
#' M <- matrix(c(0.2, 3, 2, 0.1), 2, 2)
#' out <- ot_partial_penalized(M, discard_penalty = 0.5)
#' out$transported_mass
#' rfugw_value(out)
#' @export
ot_partial_penalized <- function(
    M,
    p = NULL,
    q = NULL,
    discard_penalty,
    max_iter = 20000L,
    tol = 1e-12) {
  M <- .validate_finite_matrix(M, "M")
  if (any(M < 0)) {
    stop("`M` must be nonnegative for penalized partial transport.", call. = FALSE)
  }
  if (!is.numeric(discard_penalty) || length(discard_penalty) != 1L ||
      !is.finite(discard_penalty) || discard_penalty < 0) {
    stop("`discard_penalty` must be one finite nonnegative number.", call. = FALSE)
  }
  max_iter <- .validate_count(max_iter, "max_iter")
  tol <- .validate_positive_scalar(tol, "tol")
  ns <- nrow(M)
  nt <- ncol(M)
  p <- .validate_finite_measure(p, ns, "p", rep(1 / ns, ns))
  q <- .validate_finite_measure(q, nt, "q", rep(1 / nt, nt))
  source_mass <- sum(p)
  target_mass <- sum(q)
  augmented_total <- source_mass + target_mass

  if (augmented_total == 0) {
    out <- list(
      plan = matrix(0, ns, nt),
      penalized_partial_objective = 0,
      ot_dist = 0,
      transport_term = 0,
      discard_penalty_term = 0,
      transported_mass = 0,
      discarded_source_mass = 0,
      discarded_target_mass = 0,
      source_measure_original = p,
      target_measure_original = q,
      original_source_mass = 0,
      original_target_mass = 0,
      discard_penalty = discard_penalty,
      formulation = "ot_partial_penalized",
      backend = "analytic_zero_measures",
      regularization = 0,
      status = "converged",
      converged = TRUE,
      termination_reason = "zero_measures",
      iterations = 0L,
      max_iter = max_iter,
      residual = 0,
      error = 0,
      row_residual = 0,
      col_residual = 0,
      feasibility = "partial_variable_mass",
      feasibility_residual = 0,
      feasibility_tolerance = tol,
      feasible = TRUE,
      objective_recomputed = 0,
      objective_residual = 0,
      objective_tolerance = tol,
      objective_consistent = TRUE,
      objective_components_consistent = TRUE,
      primal_objective = 0,
      dual_objective = 0,
      duality_gap = 0,
      duality_gap_tolerance = tol,
      lp_ok = TRUE,
      warning_payload = NULL,
      mass = 0,
      mass_target = NA_real_,
      mass_residual = NA_real_,
      mass_certified = TRUE,
      mass_certification = "optimized_zero_measure"
    )
    out$runtime_provenance <- .solver_runtime_provenance(out)
    class(out) <- "rfugw_result"
    return(out)
  }

  augmented_cost <- rbind(
    cbind(M, rep(discard_penalty, ns)),
    c(rep(discard_penalty, nt), 0)
  )
  augmented_source <- c(p, target_mass)
  augmented_target <- c(q, source_mass)
  native <- cpp_ot_emd(
    M = augmented_cost,
    p = augmented_source / augmented_total,
    q = augmented_target / augmented_total,
    max_iter = max_iter,
    tol = tol
  )
  augmented_plan <- native$plan * augmented_total
  plan <- augmented_plan[seq_len(ns), seq_len(nt), drop = FALSE]
  transported_mass <- sum(plan)
  discarded_source <- source_mass - transported_mass
  discarded_target <- target_mass - transported_mass
  transport_term <- ot_linear_cost(M, plan)
  penalty_term <- discard_penalty * (discarded_source + discarded_target)
  objective <- transport_term + penalty_term
  augmented_primal <- ot_linear_cost(augmented_cost, augmented_plan)
  objective_tolerance <- max(
    tol * augmented_total,
    1e-10 * (1 + abs(objective) + abs(augmented_primal))
  )
  row_residual <- max(c(rowSums(plan) - p, 0))
  col_residual <- max(c(colSums(plan) - q, 0))
  mass_bounds_residual <- max(
    c(-transported_mass, transported_mass - min(source_mass, target_mass), 0)
  )
  feasibility_tolerance <- max(tol * augmented_total, 1e-10 * augmented_total)
  feasibility_residual <- max(row_residual, col_residual, mass_bounds_residual)
  objective_residual <- abs(objective - augmented_primal)
  scaled_dual <- native$dual_objective * augmented_total
  scaled_gap <- augmented_primal - scaled_dual
  scaled_gap_tolerance <- native$duality_gap_tolerance * augmented_total
  certificate_ok <- isTRUE(native$lp_ok) &&
    is.finite(feasibility_residual) &&
    feasibility_residual <= feasibility_tolerance &&
    objective_residual <= objective_tolerance &&
    is.finite(scaled_gap) &&
    scaled_gap >= -scaled_gap_tolerance &&
    abs(scaled_gap) <= scaled_gap_tolerance

  native$augmented_plan <- augmented_plan
  native$plan <- plan
  native$augmented_cost <- augmented_cost
  native$augmented_source_measure <- augmented_source
  native$augmented_target_measure <- augmented_target
  native$augmented_normalization_scale <- augmented_total
  native$normalized_augmented_primal_objective <- native$primal_objective
  native$normalized_augmented_dual_objective <- native$dual_objective
  native$normalized_augmented_duality_gap <- native$duality_gap
  native$primal_objective <- augmented_primal
  native$dual_objective <- scaled_dual
  native$duality_gap <- scaled_gap
  native$duality_gap_tolerance <- scaled_gap_tolerance
  native$ot_dist <- objective
  native$penalized_partial_objective <- objective
  native$transport_term <- transport_term
  native$discard_penalty_term <- penalty_term
  native$discard_penalty <- discard_penalty
  native$transported_mass <- transported_mass
  native$discarded_source_mass <- discarded_source
  native$discarded_target_mass <- discarded_target
  native$source_measure_original <- p
  native$target_measure_original <- q
  native$original_source_mass <- source_mass
  native$original_target_mass <- target_mass
  native$effective_source_mass <- source_mass
  native$effective_target_mass <- target_mass
  native$measure_normalization <- "finite_measure_preserved"
  native$penalty_direction <- "larger_discard_penalty_weakly_favors_more_transport"
  native$formulation <- "ot_partial_penalized"
  native$backend <- "cpp_transport_symmetric_dummy_reduction"
  native$regularization <- 0
  native$row_residual <- row_residual
  native$col_residual <- col_residual
  native$mass <- transported_mass
  native$mass_target <- NA_real_
  native$mass_residual <- NA_real_
  native$mass_certified <- certificate_ok
  native$mass_certification <- "optimized_and_bounded"
  native$feasibility <- "partial_variable_mass"
  native$feasibility_residual <- feasibility_residual
  native$feasibility_tolerance <- feasibility_tolerance
  native$feasible <- feasibility_residual <= feasibility_tolerance
  native$objective_recomputed <- objective
  native$objective_residual <- objective_residual
  native$objective_tolerance <- objective_tolerance
  native$objective_consistent <- objective_residual <= objective_tolerance
  native$objective_components_consistent <- all(is.finite(c(
    transport_term, penalty_term, objective
  ))) && abs(transport_term + penalty_term - objective) <= objective_tolerance
  native$residual <- max(feasibility_residual, objective_residual, abs(scaled_gap))
  native$error <- native$residual
  native$iterations <- as.integer(native$iterations)
  native$max_iter <- max_iter
  native$converged <- certificate_ok
  native$status <- if (certificate_ok) "converged" else if (!isTRUE(native$lp_ok)) {
    native$termination_reason
  } else if (!native$feasible) {
    "infeasible"
  } else {
    "objective_mismatch"
  }
  native$termination_reason <- if (certificate_ok) {
    "optimal_augmented_transport"
  } else {
    native$termination_reason
  }
  native$warning_payload <- if (certificate_ok) NULL else list(
    code = native$status,
    message = sprintf("Penalized partial OT is not certified: %s.", native$status)
  )
  native$runtime_provenance <- .solver_runtime_provenance(native)
  class(native) <- unique(c("rfugw_result", class(native)))
  native
}

#' KL-unbalanced entropic optimal transport
#'
#' @inheritParams ot_sinkhorn
#' @param p Source finite nonnegative measure (default uniform).
#' @param q Target finite nonnegative measure (default uniform).
#' @param rho Marginal generalized-KL penalties, length 1 or 2.
#' @param normalization Measure normalization policy. `"none"` preserves the
#'   supplied finite measures. `"joint"` applies one common factor so the mean
#'   source/target mass is one, preserving their mass ratio. `"separate"`
#'   normalizes each measure to mass one and is the backward-compatible default
#'   during the 0.1 transition.
#' @param init_plan Optional warm start.
#' @param method `"scaling"`, genuine `"log"`, or `"auto"`. Auto uses the
#'   documented dynamic-range criterion from `ot_sinkhorn()`.
#' @return An `rfugw_result`. `ot_dist` remains the unregularized transport
#'   term for compatibility. `regularized_objective`, all three generalized-KL
#'   terms, original/effective measures and masses, transported mass,
#'   normalization policy, and a fixed-point certificate are also returned.
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
    method = c("scaling", "log", "auto"),
    normalization = c("separate", "none", "joint")) {
  requested_method <- match.arg(method)
  normalization_defaulted <- missing(normalization)
  normalization <- match.arg(normalization)
  dat <- .prepare_unbalanced_linear_ot(M, p, q, normalization)
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
  zero_measure <- sum(dat$p) == 0 || sum(dat$q) == 0
  out <- if (zero_measure) {
    list(
      plan = matrix(0, nrow(dat$M), ncol(dat$M)),
      ot_dist = 0,
      iterations = 0L,
      error = 0,
      inner_residual = 0,
      max_inner_residual = 0,
      inner_iterations = 0L,
      inner_converged = TRUE,
      inner_status = "zero_measure_closed_form"
    )
  } else if (identical(dispatch$effective, "log")) {
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
  out$backend <- if (zero_measure) {
    "r_zero_measure_closed_form"
  } else if (identical(dispatch$effective, "log")) {
    "r_log"
  } else {
    "cpp_scaling"
  }
  out$normalization_defaulted <- normalization_defaulted
  out$normalization_source <- if (normalization_defaulted) {
    "backward_compatible_0.1_default"
  } else {
    "explicit_argument"
  }
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
  if (zero_measure) ans$termination_reason <- "zero_measure_closed_form"
  .attach_unbalanced_measure_contract(ans, dat, rho, epsilon, tol)
}
