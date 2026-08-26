.solver_protocol_version <- "1.0"

.require_mass_policy <- function(mass_policy, choices) {
  if (is.null(mass_policy) || length(mass_policy) != 1L ||
      !is.character(mass_policy) || is.na(mass_policy)) {
    stop(
      sprintf(
        "`mass_policy` must be explicit: choose one of %s.",
        paste(sprintf("\"%s\"", choices), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  match.arg(mass_policy, choices)
}

.protocol_probability_measures <- function(M, p, q, mass_policy) {
  M <- .validate_finite_matrix(M, "M")
  ns <- nrow(M)
  nt <- ncol(M)
  source <- .validate_finite_measure(p, ns, "p", rep(1 / ns, ns))
  target <- .validate_finite_measure(q, nt, "q", rep(1 / nt, nt))
  totals <- c(source = sum(source), target = sum(target))
  if (any(totals <= 0)) {
    stop(
      "Probability-measure problems require positive source and target mass.",
      call. = FALSE
    )
  }
  mass_policy <- .require_mass_policy(
    mass_policy, c("probability", "normalize")
  )
  if (identical(mass_policy, "probability") &&
      any(abs(totals - 1) > 1e-12)) {
    stop(
      paste0(
        "`mass_policy = \"probability\"` requires each supplied measure to ",
        "sum to one; use `mass_policy = \"normalize\"` to request separate ",
        "normalization explicitly."
      ),
      call. = FALSE
    )
  }
  effective_source <- source / totals[["source"]]
  effective_target <- target / totals[["target"]]
  list(
    cost = unname(M),
    source_original = source,
    target_original = target,
    source = effective_source,
    target = effective_target,
    original_mass = totals,
    effective_mass = c(source = 1, target = 1),
    mass_policy = mass_policy,
    normalization = if (identical(mass_policy, "probability")) {
      "none_probability_verified"
    } else {
      "separate_probability"
    }
  )
}

.new_transport_problem <- function(
    estimand, solver, measures, controls, objective, supports_state = FALSE) {
  structure(
    list(
      protocol_version = .solver_protocol_version,
      estimand = estimand,
      solver = solver,
      orientation = "rows_source_columns_target",
      cost = measures$cost,
      source_measure_original = measures$source_original,
      target_measure_original = measures$target_original,
      source_measure = measures$source,
      target_measure = measures$target,
      original_mass = unname(measures$original_mass),
      effective_mass = unname(measures$effective_mass),
      mass_policy = measures$mass_policy,
      normalization = measures$normalization,
      objective = objective,
      controls = controls,
      supports_state = isTRUE(supports_state)
    ),
    class = c("rfugw_transport_problem", "list")
  )
}

#' Construct an explicit balanced entropic transport problem
#'
#' These problem constructors are the versioned client boundary for
#' [transport_solve()]. Rows of `M` are source support and columns are target
#' support. `mass_policy` is mandatory so a downstream caller cannot
#' accidentally confuse already-normalized probabilities with weights that it
#' intends rfugw to normalize.
#'
#' @param M Finite source-by-target cost matrix.
#' @param p,q Nonnegative source and target weights.
#' @param epsilon Positive coefficient of
#'   `KL(plan || p %o% q)`, in the same units as `M`.
#' @param method Requested Sinkhorn method: `"scaling"`, `"log"`, or
#'   `"auto"`.
#' @param max_iter Maximum solver iterations.
#' @param tol Requested certificate tolerance.
#' @param mass_policy For balanced and partial problems, either
#'   `"probability"` (verify that both sums are one) or `"normalize"`
#'   (separately normalize and record the original masses). For KL-unbalanced
#'   problems, see [transport_problem_unbalanced()].
#' @return An `rfugw_transport_problem` for [transport_solve()].
#' @export
transport_problem_sinkhorn <- function(
    M, p = NULL, q = NULL, epsilon = 0.05,
    method = c("auto", "scaling", "log"), max_iter = 1000L, tol = 1e-9,
    mass_policy = NULL) {
  measures <- .protocol_probability_measures(M, p, q, mass_policy)
  controls <- list(
    epsilon = .validate_positive_scalar(epsilon, "epsilon"),
    method = match.arg(method),
    max_iter = .validate_count(max_iter, "max_iter"),
    tol = .validate_positive_scalar(tol, "tol")
  )
  .new_transport_problem(
    estimand = "balanced_entropic_linear_ot",
    solver = "ot_sinkhorn",
    measures = measures,
    controls = controls,
    objective = list(
      advertised = "transport_cost",
      certified = "transport_cost_plus_epsilon_product_measure_kl",
      units = "cost_times_probability_mass",
      entropy_reference = "source_measure_tensor_target_measure"
    ),
    supports_state = TRUE
  )
}

#' Construct an explicit exact balanced transport problem
#'
#' @inheritParams transport_problem_sinkhorn
#' @return An `rfugw_transport_problem` for [transport_solve()].
#' @export
transport_problem_emd <- function(
    M, p = NULL, q = NULL, max_iter = 20000L, tol = 1e-12,
    mass_policy = NULL) {
  measures <- .protocol_probability_measures(M, p, q, mass_policy)
  controls <- list(
    max_iter = .validate_count(max_iter, "max_iter"),
    tol = .validate_positive_scalar(tol, "tol")
  )
  .new_transport_problem(
    estimand = "exact_balanced_linear_ot",
    solver = "ot_emd",
    measures = measures,
    controls = controls,
    objective = list(
      advertised = "transport_cost",
      certified = "transport_cost_with_primal_dual_lp_certificate",
      units = "cost_times_probability_mass",
      entropy_reference = NA_character_
    )
  )
}

#' Construct an explicit exact fixed-mass partial transport problem
#'
#' @inheritParams transport_problem_sinkhorn
#' @param mass Transported probability mass in `[0, 1]`.
#' @return An `rfugw_transport_problem` for [transport_solve()].
#' @export
transport_problem_partial_emd <- function(
    M, p = NULL, q = NULL, mass = 1, max_iter = 20000L, tol = 1e-12,
    mass_policy = NULL) {
  measures <- .protocol_probability_measures(M, p, q, mass_policy)
  if (any(measures$cost < 0)) {
    stop("`M` must be nonnegative for exact partial transport.", call. = FALSE)
  }
  if (!is.numeric(mass) || length(mass) != 1L || !is.finite(mass) ||
      mass < 0 || mass > 1) {
    stop("`mass` must be one finite number in [0, 1].", call. = FALSE)
  }
  controls <- list(
    mass = unname(mass),
    max_iter = .validate_count(max_iter, "max_iter"),
    tol = .validate_positive_scalar(tol, "tol")
  )
  .new_transport_problem(
    estimand = "exact_fixed_mass_partial_linear_ot",
    solver = "ot_partial_emd",
    measures = measures,
    controls = controls,
    objective = list(
      advertised = "partial_transport_cost",
      certified = "partial_transport_cost_with_augmented_lp_certificate",
      units = "cost_times_probability_mass",
      entropy_reference = NA_character_
    )
  )
}

#' Construct an explicit entropic fixed-mass partial transport problem
#'
#' @inheritParams transport_problem_sinkhorn
#' @param mass Transported probability mass in `[0, 1]`.
#' @param check_every Iteration interval for diagnostics.
#' @return An `rfugw_transport_problem` for [transport_solve()].
#' @export
transport_problem_partial_sinkhorn <- function(
    M, p = NULL, q = NULL, mass = 1, epsilon = 0.1,
    method = c("auto", "scaling", "log"), max_iter = 10000L, tol = 1e-9,
    check_every = 10L, mass_policy = NULL) {
  measures <- .protocol_probability_measures(M, p, q, mass_policy)
  if (!is.numeric(mass) || length(mass) != 1L || !is.finite(mass) ||
      mass < 0 || mass > 1) {
    stop("`mass` must be one finite number in [0, 1].", call. = FALSE)
  }
  controls <- list(
    mass = unname(mass),
    epsilon = .validate_positive_scalar(epsilon, "epsilon"),
    method = match.arg(method),
    max_iter = .validate_count(max_iter, "max_iter"),
    tol = .validate_positive_scalar(tol, "tol"),
    check_every = .validate_count(check_every, "check_every")
  )
  .new_transport_problem(
    estimand = "entropic_fixed_mass_partial_linear_ot",
    solver = "ot_partial_sinkhorn",
    measures = measures,
    controls = controls,
    objective = list(
      advertised = "transport_plus_counting_measure_entropy_minus_one",
      certified = "partial_entropic_primal_dual_objective",
      units = "cost_times_probability_mass",
      entropy_reference = "counting_measure"
    ),
    supports_state = TRUE
  )
}

#' Construct an explicit KL-unbalanced transport problem
#'
#' Unlike the probability-only problem constructors, this constructor can
#' preserve finite measures. The required `mass_policy` maps to the direct
#' solver without a hidden default: `"finite_measure"` preserves both
#' measures, `"joint"` applies one common scale, and
#' `"separate_probability"` normalizes them separately.
#'
#' @inheritParams transport_problem_sinkhorn
#' @param rho One or two positive marginal-KL penalties, in cost units.
#' @param mass_policy One of `"finite_measure"`, `"joint"`, or
#'   `"separate_probability"`.
#' @return An `rfugw_transport_problem` for [transport_solve()].
#' @export
transport_problem_unbalanced <- function(
    M, p = NULL, q = NULL, epsilon = 0.05, rho = 10,
    method = c("auto", "scaling", "log"), max_iter = 500L, tol = 1e-7,
    mass_policy = NULL) {
  mass_policy <- .require_mass_policy(
    mass_policy, c("finite_measure", "joint", "separate_probability")
  )
  normalization <- switch(
    mass_policy,
    finite_measure = "none",
    joint = "joint",
    separate_probability = "separate"
  )
  prepared <- .prepare_unbalanced_linear_ot(M, p, q, normalization)
  if (length(rho) == 1L) rho <- rep(rho, 2L)
  if (!is.numeric(rho) || length(rho) != 2L || any(!is.finite(rho)) ||
      any(rho <= 0)) {
    stop("`rho` must be one or two finite positive numbers.", call. = FALSE)
  }
  measures <- list(
    cost = prepared$M,
    source_original = prepared$original_p,
    target_original = prepared$original_q,
    source = prepared$p,
    target = prepared$q,
    original_mass = prepared$original_mass,
    effective_mass = prepared$effective_mass,
    mass_policy = mass_policy,
    normalization = normalization
  )
  controls <- list(
    epsilon = .validate_positive_scalar(epsilon, "epsilon"),
    rho = unname(as.numeric(rho)),
    method = match.arg(method),
    max_iter = .validate_count(max_iter, "max_iter"),
    tol = .validate_positive_scalar(tol, "tol")
  )
  .new_transport_problem(
    estimand = "kl_unbalanced_entropic_linear_ot",
    solver = "ot_sinkhorn_unbalanced",
    measures = measures,
    controls = controls,
    objective = list(
      advertised = "finite_measure_transport_plus_three_generalized_kl_terms",
      certified = "finite_measure_regularized_objective_and_fixed_point",
      units = "cost_times_measure_mass",
      entropy_reference = "source_measure_tensor_target_measure"
    )
  )
}

.revalidate_transport_problem <- function(problem) {
  if (!inherits(problem, "rfugw_transport_problem") || !is.list(problem)) {
    stop(
      "`problem` must come from an exported `transport_problem_*()` constructor.",
      call. = FALSE
    )
  }
  if (!identical(problem$protocol_version, .solver_protocol_version)) {
    stop(
      sprintf(
        "Unsupported transport-problem protocol version `%s`; this rfugw build supports `%s`.",
        problem$protocol_version %||% "missing", .solver_protocol_version
      ),
      call. = FALSE
    )
  }
  c <- problem$controls
  switch(
    problem$solver,
    ot_sinkhorn = transport_problem_sinkhorn(
      problem$cost, problem$source_measure_original,
      problem$target_measure_original, c$epsilon, c$method, c$max_iter,
      c$tol, problem$mass_policy
    ),
    ot_emd = transport_problem_emd(
      problem$cost, problem$source_measure_original,
      problem$target_measure_original, c$max_iter, c$tol,
      problem$mass_policy
    ),
    ot_partial_emd = transport_problem_partial_emd(
      problem$cost, problem$source_measure_original,
      problem$target_measure_original, c$mass, c$max_iter, c$tol,
      problem$mass_policy
    ),
    ot_partial_sinkhorn = transport_problem_partial_sinkhorn(
      problem$cost, problem$source_measure_original,
      problem$target_measure_original, c$mass, c$epsilon, c$method,
      c$max_iter, c$tol, c$check_every, problem$mass_policy
    ),
    ot_sinkhorn_unbalanced = transport_problem_unbalanced(
      problem$cost, problem$source_measure_original,
      problem$target_measure_original, c$epsilon, c$rho, c$method,
      c$max_iter, c$tol, problem$mass_policy
    ),
    stop(
      sprintf("Unsupported protocol solver `%s`.", problem$solver %||% "missing"),
      call. = FALSE
    )
  )
}

.state_rejection <- function(code, message) {
  structure(list(code = code, message = message), class = "rfugw_state_rejection")
}

.new_solver_state <- function(type, formulation, payload, dimensions,
                              source_support, target_support, origin,
                              certified) {
  structure(
    list(
      protocol_version = .solver_protocol_version,
      state_type = type,
      formulation = formulation,
      certified = isTRUE(certified),
      dimensions = as.integer(dimensions),
      source_support = as.logical(source_support),
      target_support = as.logical(target_support),
      origin = origin,
      payload = payload
    ),
    class = c("rfugw_solver_state", "list")
  )
}

#' Extract reusable solver state
#'
#' Returned state is opaque, versioned, serializable with `saveRDS()`, and
#' certified for reuse only when the originating solve converged. Balanced
#' dual-potential state supports scaling/log interoperability and epsilon
#' continuation. Partial-Sinkhorn Dykstra state is intentionally bound to the
#' identical problem and effective method.
#'
#' @param x An `rfugw_result` or `rfugw_solver_state`.
#' @return An `rfugw_solver_state`.
#' @export
rfugw_state <- function(x) {
  if (inherits(x, "rfugw_solver_state")) return(x)
  if (!inherits(x, "rfugw_result")) {
    stop("`x` must be an rfugw result or solver state.", call. = FALSE)
  }
  if (!isTRUE(x$converged)) {
    stop("Reusable state requires a certified converged result.", call. = FALSE)
  }
  if (identical(x$formulation, "ot_sinkhorn")) {
    payload <- x$dual_state
    if (!is.list(payload) || is.null(payload$source) || is.null(payload$target)) {
      stop("The balanced result does not contain reusable dual state.", call. = FALSE)
    }
    return(.new_solver_state(
      type = "balanced_sinkhorn_dual_potentials_v1",
      formulation = "balanced_entropic_linear_ot",
      payload = payload,
      dimensions = c(length(payload$source), length(payload$target)),
      source_support = payload$source_support %||% rep(TRUE, length(payload$source)),
      target_support = payload$target_support %||% rep(TRUE, length(payload$target)),
      origin = list(
        solver = "ot_sinkhorn",
        backend = x$backend,
        requested_method = x$requested_sinkhorn_method,
        effective_method = x$effective_sinkhorn_method,
        epsilon = x$regularization,
        status = x$status
      ),
      certified = TRUE
    ))
  }
  if (identical(x$formulation, "ot_partial_sinkhorn")) {
    payload <- x$warm_state
    if (!is.list(payload) || !isTRUE(payload$certified)) {
      stop("The partial result does not contain certified Dykstra state.", call. = FALSE)
    }
    return(.new_solver_state(
      type = "partial_sinkhorn_dykstra_v1",
      formulation = "entropic_fixed_mass_partial_linear_ot",
      payload = payload,
      dimensions = dim(x$plan),
      source_support = seq_len(nrow(x$plan)) %in% payload$active_source,
      target_support = seq_len(ncol(x$plan)) %in% payload$active_target,
      origin = list(
        solver = "ot_partial_sinkhorn",
        backend = x$backend,
        requested_method = x$requested_sinkhorn_method,
        effective_method = x$effective_sinkhorn_method,
        epsilon = x$regularization,
        status = x$status
      ),
      certified = TRUE
    ))
  }
  stop(
    sprintf("Formulation `%s` does not expose reusable protocol state.",
            x$formulation %||% "unknown"),
    call. = FALSE
  )
}

.check_protocol_state <- function(problem, state) {
  if (inherits(state, "rfugw_result")) {
    state <- tryCatch(
      rfugw_state(state),
      error = function(e) .state_rejection("uncertified_result", conditionMessage(e))
    )
  }
  if (inherits(state, "rfugw_state_rejection")) {
    return(list(accepted = FALSE, state = NULL, reason = state))
  }
  reject <- function(code, message) {
    list(accepted = FALSE, state = NULL, reason = .state_rejection(code, message))
  }
  if (!inherits(state, "rfugw_solver_state") || !is.list(state)) {
    return(reject("invalid_state_class", "State must come from `rfugw_state()`."))
  }
  if (!identical(state$protocol_version, .solver_protocol_version)) {
    return(reject("unsupported_state_version", "State protocol version is incompatible."))
  }
  if (!isTRUE(state$certified)) {
    return(reject("uncertified_state", "State is not certified for reuse."))
  }
  expected_dim <- c(length(problem$source_measure), length(problem$target_measure))
  if (!identical(as.integer(state$dimensions), as.integer(expected_dim))) {
    return(reject(
      "dimension_mismatch",
      sprintf("State dimensions must be %d x %d.", expected_dim[[1]], expected_dim[[2]])
    ))
  }
  if (!identical(state$source_support, problem$source_measure > 0) ||
      !identical(state$target_support, problem$target_measure > 0)) {
    return(reject(
      "support_mismatch", "State support is incompatible with the requested measures."
    ))
  }

  if (identical(problem$solver, "ot_sinkhorn")) {
    if (!identical(state$state_type, "balanced_sinkhorn_dual_potentials_v1") ||
        !identical(state$formulation, problem$estimand)) {
      return(reject(
        "formulation_mismatch", "State is not balanced-Sinkhorn dual state."
      ))
    }
    payload <- state$payload
    if (!is.list(payload) || !is.numeric(payload$source) ||
        !is.numeric(payload$target) ||
        length(payload$source) != expected_dim[[1]] ||
        length(payload$target) != expected_dim[[2]]) {
      return(reject("malformed_payload", "Balanced dual payload has invalid dimensions."))
    }
    if (any(!is.finite(c(payload$source, payload$target)))) {
      return(reject("nonfinite_payload", "Balanced dual potentials must be finite."))
    }
    dispatch <- .select_sinkhorn_method(
      problem$controls$method, problem$cost, problem$controls$epsilon,
      context = "Balanced Sinkhorn"
    )
    if (identical(dispatch$effective, "scaling")) {
      exponent_max <- max(c(payload$source, payload$target)) /
        problem$controls$epsilon
      if (!is.finite(exponent_max) ||
          exponent_max > log(.Machine$double.xmax) - 2) {
        return(reject(
          "backend_overflow",
          "State potentials overflow the requested scaling backend; use `method = \"log\"` or `\"auto\"`."
        ))
      }
    }
    return(list(accepted = TRUE, state = state, reason = NULL))
  }

  if (identical(problem$solver, "ot_partial_sinkhorn")) {
    if (!identical(state$state_type, "partial_sinkhorn_dykstra_v1") ||
        !identical(state$formulation, problem$estimand)) {
      return(reject(
        "formulation_mismatch", "State is not partial-Sinkhorn Dykstra state."
      ))
    }
    active <- .partial_sinkhorn_active_problem(
      problem$cost, problem$source_measure, problem$target_measure,
      problem$controls$mass, problem$controls$epsilon
    )
    dispatch <- .linear_partial_sinkhorn_dispatch(
      problem$controls$method, active$cost, problem$controls$epsilon
    )
    check <- tryCatch(
      {
        .validate_partial_sinkhorn_state(state$payload, active, dispatch$effective)
        NULL
      },
      error = function(e) conditionMessage(e)
    )
    if (!is.null(check)) {
      code <- if (grepl("method", check, fixed = TRUE)) {
        "backend_mismatch"
      } else if (grepl("finite", check, fixed = TRUE)) {
        "nonfinite_payload"
      } else {
        "problem_mismatch"
      }
      return(reject(code, check))
    }
    return(list(accepted = TRUE, state = state, reason = NULL))
  }

  reject(
    "state_unsupported",
    sprintf("Formulation `%s` does not support iterative protocol state.",
            problem$estimand)
  )
}

.attach_protocol_result <- function(result, problem, state_check) {
  result$solver_protocol_version <- .solver_protocol_version
  result$solver_problem <- problem
  result$warm_started <- isTRUE(state_check$accepted)
  result$warm_start_accepted <- isTRUE(state_check$accepted)
  result$warm_start_rejection <- state_check$reason
  result$source_measure_original <- problem$source_measure_original
  result$target_measure_original <- problem$target_measure_original
  result$source_measure_effective <- problem$source_measure
  result$target_measure_effective <- problem$target_measure
  result$original_source_mass <- problem$original_mass[[1L]]
  result$original_target_mass <- problem$original_mass[[2L]]
  result$effective_source_mass <- problem$effective_mass[[1L]]
  result$effective_target_mass <- problem$effective_mass[[2L]]
  result$measure_normalization <- problem$normalization
  result$mass_policy <- problem$mass_policy
  result$requested_controls <- problem$controls
  result$effective_controls <- problem$controls
  result$effective_controls$method <-
    result$effective_sinkhorn_method %||% NA_character_
  result$effective_controls$tol <-
    result$feasibility_tolerance %||% problem$controls$tol
  result$effective_controls$max_iter <- result$max_iter
  result$effective_controls$backend <- result$backend
  result$plan_representation <- transport_plan_representation(rfugw_plan(result))
  class(result) <- unique(c("rfugw_result", class(result)))
  result
}

#' Solve a versioned transport problem
#'
#' `transport_solve()` validates the complete problem and any reusable state
#' before entering a numerical backend. When `state_policy = "cold"`, an
#' incompatible state is rejected with a structured reason and the problem is
#' solved cold; `"error"` fails closed instead. A state is never inherited as
#' evidence that the new result converged.
#'
#' @param problem An `rfugw_transport_problem`.
#' @param init_state Optional state returned by [rfugw_state()] or a certified
#'   prior `rfugw_result`.
#' @param state_policy What to do with invalid or unsupported state: error or
#'   explicitly record rejection and solve cold.
#' @return An `rfugw_result` with protocol, problem, initialization, mass,
#'   control, plan-representation, status, and runtime provenance.
#' @export
transport_solve <- function(
    problem, init_state = NULL, state_policy = c("error", "cold")) {
  state_policy <- match.arg(state_policy)
  problem <- .revalidate_transport_problem(problem)
  state_check <- if (is.null(init_state)) {
    list(
      accepted = FALSE,
      state = NULL,
      reason = .state_rejection("not_requested", "No warm state was requested.")
    )
  } else {
    .check_protocol_state(problem, init_state)
  }
  if (!is.null(init_state) && !isTRUE(state_check$accepted) &&
      identical(state_policy, "error")) {
    stop(
      sprintf(
        "Warm state rejected [%s]: %s",
        state_check$reason$code, state_check$reason$message
      ),
      call. = FALSE
    )
  }
  payload <- if (isTRUE(state_check$accepted)) state_check$state$payload else NULL
  c <- problem$controls
  result <- switch(
    problem$solver,
    ot_sinkhorn = ot_sinkhorn(
      problem$cost, problem$source_measure, problem$target_measure,
      epsilon = c$epsilon, method = c$method, max_iter = c$max_iter,
      tol = c$tol, init_duals = payload
    ),
    ot_emd = ot_emd(
      problem$cost, problem$source_measure, problem$target_measure,
      max_iter = c$max_iter, tol = c$tol
    ),
    ot_partial_emd = ot_partial_emd(
      problem$cost, problem$source_measure, problem$target_measure,
      mass = c$mass, max_iter = c$max_iter, tol = c$tol
    ),
    ot_partial_sinkhorn = ot_partial_sinkhorn(
      problem$cost, problem$source_measure, problem$target_measure,
      mass = c$mass, epsilon = c$epsilon, method = c$method,
      max_iter = c$max_iter, tol = c$tol, check_every = c$check_every,
      init_state = payload
    ),
    ot_sinkhorn_unbalanced = ot_sinkhorn_unbalanced(
      problem$cost, problem$source_measure_original,
      problem$target_measure_original, epsilon = c$epsilon, rho = c$rho,
      method = c$method, max_iter = c$max_iter, tol = c$tol,
      normalization = problem$normalization
    )
  )
  .attach_protocol_result(result, problem, state_check)
}

#' Extract the explicit transport problem contract
#'
#' @param x An `rfugw_transport_problem` or a result from [transport_solve()].
#' @return An `rfugw_transport_problem`.
#' @export
rfugw_problem <- function(x) {
  if (inherits(x, "rfugw_transport_problem")) return(x)
  if (inherits(x, "rfugw_result") &&
      inherits(x$solver_problem, "rfugw_transport_problem")) {
    return(x$solver_problem)
  }
  stop("`x` does not contain a solver-protocol problem contract.", call. = FALSE)
}

#' Extract source, target, and transported-mass provenance
#'
#' @param x An `rfugw_transport_problem` or `rfugw_result`.
#' @return A list separating input and effective measures and masses.
#' @export
rfugw_masses <- function(x) {
  problem <- if (inherits(x, "rfugw_transport_problem")) x else {
    tryCatch(rfugw_problem(x), error = function(e) NULL)
  }
  if (!is.null(problem)) {
    transported <- if (inherits(x, "rfugw_result")) {
      transport_plan_mass(rfugw_plan(x), "total")
    } else {
      NA_real_
    }
    return(list(
      source = list(
        input = problem$source_measure_original,
        effective = problem$source_measure,
        input_mass = problem$original_mass[[1L]],
        effective_mass = problem$effective_mass[[1L]]
      ),
      target = list(
        input = problem$target_measure_original,
        effective = problem$target_measure,
        input_mass = problem$original_mass[[2L]],
        effective_mass = problem$effective_mass[[2L]]
      ),
      transported_mass = transported,
      mass_policy = problem$mass_policy,
      normalization = problem$normalization
    ))
  }
  if (!inherits(x, "rfugw_result")) {
    stop("`x` must be a transport problem or rfugw result.", call. = FALSE)
  }
  plan <- rfugw_plan(x)
  source_mass <- transport_plan_mass(plan, "source")
  target_mass <- transport_plan_mass(plan, "target")
  total_mass <- transport_plan_mass(plan, "total")
  list(
    source = list(
      input = x$source_measure_original,
      effective = x$source_measure_effective %||% source_mass,
      input_mass = x$original_source_mass %||% sum(source_mass),
      effective_mass = x$effective_source_mass %||% sum(source_mass)
    ),
    target = list(
      input = x$target_measure_original,
      effective = x$target_measure_effective %||% target_mass,
      input_mass = x$original_target_mass %||% sum(target_mass),
      effective_mass = x$effective_target_mass %||% sum(target_mass)
    ),
    transported_mass = total_mass,
    mass_policy = x$mass_policy %||% NA_character_,
    normalization = x$measure_normalization %||% NA_character_
  )
}

#' Extract stable solver and runtime provenance
#'
#' @param x An `rfugw_transport_problem` or `rfugw_result`.
#' @return A list identifying the protocol, estimand, orientation, solver,
#'   backend, controls, objective convention, plan representation, and runtime.
#' @export
rfugw_provenance <- function(x) {
  problem <- rfugw_problem(x)
  is_result <- inherits(x, "rfugw_result")
  list(
    protocol_version = problem$protocol_version,
    estimand = problem$estimand,
    orientation = problem$orientation,
    solver = problem$solver,
    backend = if (is_result) x$backend else NA_character_,
    requested_controls = problem$controls,
    effective_controls = if (is_result) x$effective_controls else NULL,
    objective = problem$objective,
    mass_policy = problem$mass_policy,
    normalization = problem$normalization,
    plan_representation = if (is_result) {
      transport_plan_representation(rfugw_plan(x))
    } else {
      NULL
    },
    status = if (is_result) x$status else NA_character_,
    converged = if (is_result) x$converged else NA,
    runtime = if (is_result) x$runtime_provenance else NULL
  )
}

#' Query public transport-solver capabilities
#'
#' @param problem Optional `rfugw_transport_problem`; when supplied, return
#'   the row for its estimand.
#' @return A data frame with one row per distinct estimand family. Comma-
#'   separated fields enumerate methods and plan representations; they are
#'   descriptive, not accepted as solver-selection strings. `maturity` and
#'   `coverage_family` link every row to the executable numerical-path matrix
#'   used by release evidence.
#' @export
transport_capabilities <- function(problem = NULL) {
  capability <- function(
      estimand, public_solver, protocol_constructor, mass_policy, warm_state,
      state_compatibility, plan_representation, exact_certificate,
      relational_costs, maturity, coverage_family) {
    data.frame(
      estimand = estimand,
      public_solver = public_solver,
      protocol_constructor = protocol_constructor,
      mass_policy = mass_policy,
      warm_state = warm_state,
      state_compatibility = state_compatibility,
      plan_representation = plan_representation,
      exact_certificate = exact_certificate,
      relational_costs = relational_costs,
      maturity = maturity,
      coverage_family = coverage_family,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, list(
    capability(
      "balanced_entropic_linear_ot", "ot_sinkhorn", TRUE,
      "probability|normalize", TRUE,
      "dual potentials: scaling/log/auto and epsilon continuation", "dense",
      FALSE, FALSE, "flagship", "balanced_linear_ot"
    ),
    capability(
      "sinkhorn_divergence", "ot_sinkhorn_divergence", FALSE,
      "probability", FALSE, "component solves only", "dense", FALSE, FALSE,
      "supported", "sinkhorn_divergence"
    ),
    capability(
      "exact_balanced_linear_ot", "ot_emd", TRUE,
      "probability|normalize", FALSE, "none", "dense", TRUE, FALSE,
      "flagship", "balanced_linear_ot"
    ),
    capability(
      "exact_fixed_mass_partial_linear_ot", "ot_partial_emd", TRUE,
      "probability|normalize", FALSE, "none", "dense", TRUE, FALSE,
      "supported", "partial_linear_ot"
    ),
    capability(
      "entropic_fixed_mass_partial_linear_ot", "ot_partial_sinkhorn", TRUE,
      "probability|normalize", TRUE,
      "identical problem and effective Dykstra method", "dense", FALSE,
      FALSE, "supported", "partial_linear_ot"
    ),
    capability(
      "penalized_variable_mass_partial_linear_ot", "ot_partial_penalized",
      FALSE, "finite_measure", FALSE, "none", "dense", TRUE, FALSE,
      "supported", "penalized_partial_ot"
    ),
    capability(
      "kl_unbalanced_entropic_linear_ot", "ot_sinkhorn_unbalanced", TRUE,
      "finite_measure|joint|separate_probability", FALSE, "none", "dense",
      FALSE, FALSE, "flagship", "kl_uot"
    ),
    capability(
      "translation_invariant_kl_unbalanced_linear_ot",
      "ot_sinkhorn_unbalanced_ti", FALSE, "finite_measure", FALSE, "none",
      "implicit|sparse_edges|dense", FALSE, FALSE, "supported", "ti_kl_uot"
    ),
    capability(
      "classical_wasserstein_summary",
      "ot_wasserstein_cost|ot_wasserstein_distance", FALSE, "probability",
      FALSE, "underlying solver only", "dense", FALSE, FALSE, "supported",
      "wasserstein_summary"
    ),
    capability(
      "fixed_support_wasserstein_barycenter_weights", "ot_barycenter_weights",
      FALSE, "probability", TRUE, "component Sinkhorn states", "dense", TRUE,
      FALSE, "supported", "fixed_support_barycenter"
    ),
    capability(
      "gromov_wasserstein", "gromov_wasserstein|fgw_exact_cg", FALSE,
      "probability", FALSE, "plan starts only", "dense", FALSE, TRUE,
      "flagship", "balanced_fgw"
    ),
    capability(
      "fused_gromov_wasserstein", "fgw_entropic|fgw_exact_cg", FALSE,
      "probability", FALSE, "plan starts only", "dense", FALSE, TRUE,
      "flagship", "balanced_fgw"
    ),
    capability(
      "fixed_mass_partial_gromov_wasserstein",
      "partial_gromov_wasserstein|entropic_partial_gromov_wasserstein",
      FALSE, "probability", FALSE, "plan starts only", "dense", FALSE, TRUE,
      "supported", "partial_fgw"
    ),
    capability(
      "fixed_mass_partial_fused_gromov_wasserstein",
      "partial_fused_gromov_wasserstein|entropic_partial_fused_gromov_wasserstein",
      FALSE, "probability", FALSE, "plan starts only", "dense", FALSE, TRUE,
      "supported", "partial_fgw"
    ),
    capability(
      "penalized_variable_mass_partial_fused_gromov_wasserstein",
      "penalized_partial_fused_gromov_wasserstein", FALSE, "finite_measure",
      FALSE, "plan starts only", "dense", FALSE, TRUE, "supported",
      "penalized_partial_fgw"
    ),
    capability(
      "semirelaxed_gromov_wasserstein",
      "semirelaxed_gromov_wasserstein|entropic_semirelaxed_gromov_wasserstein",
      FALSE, "fixed_source_probability", FALSE, "plan starts only", "dense",
      FALSE, TRUE, "supported", "semirelaxed_fgw"
    ),
    capability(
      "semirelaxed_fused_gromov_wasserstein",
      "semirelaxed_fused_gromov_wasserstein|entropic_semirelaxed_fused_gromov_wasserstein",
      FALSE, "fixed_source_probability", FALSE, "plan starts only", "dense",
      FALSE, TRUE, "supported", "semirelaxed_fgw"
    ),
    capability(
      "fused_unbalanced_gromov_wasserstein", "fugw_kl", FALSE,
      "finite_measure", FALSE, "plan starts only", "dense", FALSE, TRUE,
      "flagship", "fugw"
    ),
    capability(
      "across_spaces_unbalanced_ot",
      "fused_unbalanced_across_spaces_divergence", FALSE, "finite_measure",
      FALSE, "none", "dense", FALSE, TRUE, "supported", "ucoot"
    ),
    capability(
      "unbalanced_co_optimal_transport", "unbalanced_co_optimal_transport",
      FALSE, "finite_measure", FALSE, "none", "dense", FALSE, TRUE,
      "supported", "ucoot"
    ),
    capability(
      "fixed_support_gromov_wasserstein_barycenter", "fgw_barycenters", FALSE,
      "probability", TRUE, "constituent plan starts", "dense", FALSE, TRUE,
      "supported", "barycenter"
    ),
    capability(
      "multi_collection_alignment", "multialign_fit", FALSE, "probability",
      TRUE, "constituent plan starts", "dense", FALSE, TRUE, "supported",
      "multialign"
    ),
    capability(
      "sampled_gromov_wasserstein",
      "sampled_gromov_wasserstein|sampled_gromov_wasserstein_coords|sampled_gw_from_graphs",
      FALSE, "probability", FALSE, "none", "dense", FALSE, TRUE,
      "experimental", "sampled_gw"
    ),
    capability(
      "dense_plan_gromov_wasserstein_approximation",
      "dense_gromov_wasserstein_plan_svd|lowrank_gromov_wasserstein_samples",
      FALSE, "probability", FALSE, "none", "dense_input_low_rank_output",
      FALSE, TRUE, "experimental", "lowrank_gw"
    )
  ))
  row.names(out) <- NULL
  if (is.null(problem)) return(out)
  problem <- .revalidate_transport_problem(problem)
  out[out$estimand == problem$estimand, , drop = FALSE]
}

.validate_capability_path_matrix <- function(
    capabilities = transport_capabilities(), path_matrix) {
  capability_fields <- c("estimand", "maturity", "coverage_family")
  path_fields <- c(
    "family", "dimension", "path", "comparison", "scope", "maturity"
  )
  missing_capability <- setdiff(capability_fields, names(capabilities))
  missing_path <- setdiff(path_fields, names(path_matrix))
  if (length(missing_capability)) {
    stop(
      "Capability inventory is missing: ",
      paste(missing_capability, collapse = ", "), call. = FALSE
    )
  }
  if (length(missing_path)) {
    stop(
      "Numerical-path matrix is missing: ",
      paste(missing_path, collapse = ", "), call. = FALSE
    )
  }
  if (anyDuplicated(capabilities$estimand)) {
    stop("Capability estimands must be unique.", call. = FALSE)
  }
  if (anyDuplicated(path_matrix[c("family", "dimension", "path")])) {
    stop("Numerical-path entries must be unique.", call. = FALSE)
  }
  maturity_levels <- c("flagship", "supported", "experimental")
  if (any(!capabilities$maturity %in% maturity_levels) ||
      any(!path_matrix$maturity %in% maturity_levels)) {
    stop("Unknown capability or path maturity.", call. = FALSE)
  }
  missing_families <- setdiff(
    unique(capabilities$coverage_family), unique(path_matrix$family)
  )
  if (length(missing_families)) {
    stop(
      "Capabilities lack numerical-path evidence: ",
      paste(missing_families, collapse = ", "), call. = FALSE
    )
  }
  unowned_families <- setdiff(
    unique(path_matrix$family), unique(capabilities$coverage_family)
  )
  if (length(unowned_families)) {
    stop(
      "Numerical-path families lack capability ownership: ",
      paste(unowned_families, collapse = ", "), call. = FALSE
    )
  }
  families <- unique(capabilities$coverage_family)
  for (family in families) {
    declared <- unique(capabilities$maturity[
      capabilities$coverage_family == family
    ])
    evidenced <- unique(path_matrix$maturity[path_matrix$family == family])
    if (length(declared) != 1L || length(evidenced) != 1L ||
        !identical(declared, evidenced)) {
      stop(
        sprintf(
          "Maturity mismatch for `%s`: capability=%s, path=%s.", family,
          paste(declared, collapse = "|"), paste(evidenced, collapse = "|")
        ),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

.release_support_evidence <- function(
    capabilities = transport_capabilities(), path_matrix) {
  .validate_capability_path_matrix(capabilities, path_matrix)
  describe <- function(rows) {
    sprintf(
      "%s via `%s` (%s; coverage family `%s`)",
      gsub("_", " ", rows$estimand, fixed = TRUE), rows$public_solver,
      rows$maturity, rows$coverage_family
    )
  }
  verified <- capabilities[
    capabilities$maturity %in% c("flagship", "supported"), , drop = FALSE
  ]
  experimental <- capabilities[
    capabilities$maturity == "experimental", , drop = FALSE
  ]
  list(
    verified_support = describe(verified),
    experimental_boundaries = c(
      describe(experimental),
      paste0(
        "Directed-KL structural GW/FGW loss is deferred because an implicit ",
        "logarithm floor changes the estimand"
      ),
      paste0(
        "No end-to-end scalable relational-OT path is promoted; sampled and ",
        "dense-plan compression paths retain experimental maturity"
      )
    )
  )
}

#' Print a transport-problem contract
#'
#' @param x An `rfugw_transport_problem`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.rfugw_transport_problem <- function(x, ...) {
  cat("<rfugw_transport_problem v", x$protocol_version, ">\n", sep = "")
  cat("  estimand:   ", x$estimand, "\n", sep = "")
  cat("  solver:     ", x$solver, "\n", sep = "")
  cat("  orientation:", x$orientation, "\n", sep = "")
  cat("  shape:      ", nrow(x$cost), " x ", ncol(x$cost), "\n", sep = "")
  cat("  mass policy:", x$mass_policy, "\n", sep = "")
  invisible(x)
}

#' Print reusable solver state
#'
#' @param x An `rfugw_solver_state`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.rfugw_solver_state <- function(x, ...) {
  cat("<rfugw_solver_state v", x$protocol_version, ">\n", sep = "")
  cat("  type:       ", x$state_type, "\n", sep = "")
  cat("  formulation:", x$formulation, "\n", sep = "")
  cat("  shape:      ", x$dimensions[[1L]], " x ", x$dimensions[[2L]], "\n", sep = "")
  cat("  certified:  ", x$certified, "\n", sep = "")
  invisible(x)
}
