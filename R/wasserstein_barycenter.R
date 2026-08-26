.project_simplex_lower <- function(x, lower = 0) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x))) {
    stop("Simplex projection input must be finite and nonempty.", call. = FALSE)
  }
  if (!is.numeric(lower) || length(lower) != 1L || !is.finite(lower) ||
      lower < 0 || lower * length(x) >= 1) {
    stop("`min_weight` must be nonnegative and smaller than 1 / support size.",
         call. = FALSE)
  }
  mass <- 1 - lower * length(x)
  shifted <- x - lower
  sorted <- sort(shifted, decreasing = TRUE)
  cumulative <- cumsum(sorted)
  rho <- max(which(sorted - (cumulative - mass) / seq_along(sorted) > 0))
  theta <- (cumulative[[rho]] - mass) / rho
  pmax(shifted - theta, 0) + lower
}

.prepare_barycenter_problem <- function(costs, measures, coefficients) {
  if (!is.list(costs) || !length(costs)) {
    stop("`costs` must be a nonempty list of cost matrices.", call. = FALSE)
  }
  if (!is.list(measures) || length(measures) != length(costs)) {
    stop("`measures` must be a list with one entry per cost matrix.",
         call. = FALSE)
  }
  costs <- lapply(seq_along(costs), function(s) {
    cost <- .validate_finite_matrix(costs[[s]], sprintf("costs[[%d]]", s))
    if (any(cost < 0)) {
      stop("Wasserstein barycenter costs must be nonnegative.", call. = FALSE)
    }
    unname(cost)
  })
  support_size <- ncol(costs[[1L]])
  if (any(vapply(costs, ncol, integer(1)) != support_size)) {
    stop("Every cost matrix must have the same target-support column count.",
         call. = FALSE)
  }
  measures <- lapply(seq_along(measures), function(s) {
    .assert_prob(
      measures[[s]], nrow(costs[[s]]), sprintf("measures[[%d]]", s)
    )
  })
  if (is.null(coefficients)) {
    coefficients <- rep(1 / length(costs), length(costs))
  }
  if (!is.numeric(coefficients) || length(coefficients) != length(costs) ||
      any(!is.finite(coefficients)) || any(coefficients < 0) ||
      sum(coefficients) <= 0) {
    stop(
      "`coefficients` must be finite, nonnegative, and have positive total mass.",
      call. = FALSE
    )
  }
  coefficients <- unname(as.numeric(coefficients / sum(coefficients)))
  list(
    costs = costs,
    measures = measures,
    coefficients = coefficients,
    support_size = support_size,
    source_sizes = vapply(costs, nrow, integer(1))
  )
}

.build_exact_barycenter_lp <- function(problem) {
  component_sizes <- vapply(problem$costs, length, integer(1))
  offsets <- c(0L, cumsum(component_sizes))
  plan_variables <- sum(component_sizes)
  q_index <- plan_variables + seq_len(problem$support_size)
  variable_count <- plan_variables + problem$support_size
  constraint_count <- sum(problem$source_sizes) +
    length(problem$costs) * problem$support_size
  constraints <- matrix(0, constraint_count, variable_count)
  rhs <- numeric(constraint_count)
  directions <- rep("=", constraint_count)
  constraint <- 0L

  for (s in seq_along(problem$costs)) {
    ns <- problem$source_sizes[[s]]
    nt <- problem$support_size
    for (i in seq_len(ns)) {
      constraint <- constraint + 1L
      index <- offsets[[s]] + i + (seq_len(nt) - 1L) * ns
      constraints[constraint, index] <- 1
      rhs[[constraint]] <- problem$measures[[s]][[i]]
    }
    for (j in seq_len(nt)) {
      constraint <- constraint + 1L
      index <- offsets[[s]] + (j - 1L) * ns + seq_len(ns)
      constraints[constraint, index] <- 1
      constraints[constraint, q_index[[j]]] <- -1
    }
  }
  objective <- c(
    unlist(Map(
      function(cost, coefficient) coefficient * as.vector(cost),
      problem$costs, problem$coefficients
    ), use.names = FALSE),
    numeric(problem$support_size)
  )
  list(
    objective = objective,
    constraints = constraints,
    directions = directions,
    rhs = rhs,
    q_index = q_index,
    offsets = offsets,
    component_sizes = component_sizes,
    variable_count = variable_count,
    constraint_count = constraint_count
  )
}

.solve_exact_barycenter <- function(problem, max_iter, tol) {
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop(
      "Exact fixed-support Wasserstein barycenter weights require the suggested package `lpSolve`; install it or use `mode = \"regularized\"`.",
      call. = FALSE
    )
  }
  max_iter <- .validate_count(max_iter, "max_iter")
  tol <- .validate_positive_scalar(tol, "tol")
  lp_problem <- .build_exact_barycenter_lp(problem)
  allocation_bytes <- as.numeric(object.size(lp_problem$constraints)) +
    as.numeric(object.size(lp_problem$objective))
  lp <- lpSolve::lp(
    direction = "min",
    objective.in = lp_problem$objective,
    const.mat = lp_problem$constraints,
    const.dir = lp_problem$directions,
    const.rhs = lp_problem$rhs,
    compute.sens = 1
  )
  if (!identical(lp$status, 0L)) {
    out <- list(
      weights = rep(NA_real_, problem$support_size),
      barycenter_weights = rep(NA_real_, problem$support_size),
      component_results = vector("list", length(problem$costs)),
      component_values = rep(NA_real_, length(problem$costs)),
      coefficients = problem$coefficients,
      barycenter_objective = Inf,
      objective_recomputed = Inf,
      objective_residual = Inf,
      objective_tolerance = tol,
      objective_consistent = FALSE,
      simplex_residual = Inf,
      feasible = FALSE,
      kkt_residual = Inf,
      kkt_tolerance = tol,
      kkt_consistent = FALSE,
      status = "lp_failure",
      converged = FALSE,
      termination_reason = paste0("lpSolve_status_", lp$status),
      formulation = "fixed_support_wasserstein_barycenter_exact",
      mode = "exact",
      backend = "lpSolve_joint_lp",
      iterations = 0L,
      max_iter = max_iter,
      component_solves = 0L,
      total_component_iterations = 0L,
      allocation_evidence = list(
        joint_lp_problem_bytes = allocation_bytes,
        variable_count = lp_problem$variable_count,
        constraint_count = lp_problem$constraint_count
      )
    )
    class(out) <- c("rfugw_barycenter_result", "list")
    return(out)
  }

  raw_weights <- lp$solution[lp_problem$q_index]
  weights <- pmax(raw_weights, 0)
  weights <- weights / sum(weights)
  component_results <- Map(
    function(cost, measure) {
      ot_emd(cost, measure, weights, max_iter = max_iter, tol = tol)
    },
    problem$costs,
    problem$measures
  )
  component_values <- vapply(component_results, rfugw_value, numeric(1))
  objective <- sum(problem$coefficients * component_values)
  objective_residual <- abs(objective - lp$objval)
  objective_tolerance <- max(
    1e-9,
    50 * tol * (1 + abs(objective) + abs(lp$objval))
  )
  primal_residual <- max(abs(
    as.numeric(lp_problem$constraints %*% lp$solution) - lp_problem$rhs
  ))
  simplex_residual <- max(
    abs(sum(weights) - 1), max(c(-weights, 0)), primal_residual
  )
  constraint_duals <- head(lp$duals, lp_problem$constraint_count)
  variable_reduced_costs <- tail(lp$duals, lp_problem$variable_count)
  weight_reduced_costs <- variable_reduced_costs[lp_problem$q_index]
  active <- weights > max(1e-10, 10 * tol)
  kkt_residual <- max(
    if (any(active)) max(abs(weight_reduced_costs[active])) else 0,
    if (any(!active)) max(c(-weight_reduced_costs[!active], 0)) else 0
  )
  cost_scale <- max(vapply(problem$costs, max, numeric(1)))
  kkt_tolerance <- max(1e-8, 50 * tol * (1 + cost_scale))
  dual_objective <- sum(constraint_duals * lp_problem$rhs)
  duality_gap <- objective - dual_objective
  all_components <- all(vapply(
    component_results,
    function(x) isTRUE(x$converged) && isTRUE(x$feasible) &&
      isTRUE(x$objective_consistent),
    logical(1)
  ))
  feasible <- is.finite(simplex_residual) &&
    simplex_residual <= max(1e-9, 20 * tol)
  objective_consistent <- is.finite(objective_residual) &&
    objective_residual <= objective_tolerance
  kkt_consistent <- is.finite(kkt_residual) && kkt_residual <= kkt_tolerance &&
    is.finite(duality_gap) && abs(duality_gap) <= objective_tolerance
  converged <- all_components && feasible && objective_consistent && kkt_consistent
  status <- if (!all_components) {
    "inner_failure"
  } else if (!feasible) {
    "infeasible"
  } else if (!objective_consistent) {
    "objective_mismatch"
  } else if (!kkt_consistent) {
    "stationarity_failure"
  } else {
    "converged"
  }
  out <- list(
    weights = weights,
    barycenter_weights = weights,
    source_measures = problem$measures,
    coefficients = problem$coefficients,
    costs = problem$costs,
    component_results = component_results,
    component_values = component_values,
    component_status = vapply(component_results, rfugw_status, character(1)),
    component_residuals = vapply(
      component_results, function(x) x$residual, numeric(1)
    ),
    component_solves = length(component_results),
    total_component_iterations = sum(vapply(
      component_results, function(x) x$iterations, integer(1)
    )),
    barycenter_objective = objective,
    objective_components = problem$coefficients * component_values,
    objective_recomputed = objective,
    joint_lp_objective = lp$objval,
    objective_residual = objective_residual,
    objective_tolerance = objective_tolerance,
    objective_consistent = objective_consistent,
    primal_objective = objective,
    dual_objective = dual_objective,
    duality_gap = duality_gap,
    simplex_residual = simplex_residual,
    feasibility_residual = simplex_residual,
    feasibility_tolerance = max(1e-9, 20 * tol),
    feasible = feasible,
    weight_reduced_costs = weight_reduced_costs,
    kkt_residual = kkt_residual,
    kkt_tolerance = kkt_tolerance,
    kkt_consistent = kkt_consistent,
    objective_trace = objective,
    objective_monotone = TRUE,
    warm_started = FALSE,
    warm_start_accepted = FALSE,
    warm_state = NULL,
    status = status,
    converged = converged,
    termination_reason = if (converged) {
      "joint_lp_primal_dual_and_component_certified"
    } else {
      status
    },
    warning_payload = if (converged) NULL else list(
      code = status,
      message = sprintf("Exact barycenter result is not certified: %s.", status)
    ),
    formulation = "fixed_support_wasserstein_barycenter_exact",
    mode = "exact",
    objective_convention = "weighted_sum_exact_transport_costs",
    objective_units = "cost_times_probability_mass",
    regularization = 0,
    backend = "lpSolve_joint_lp_plus_cpp_transport_certificates",
    iterations = 1L,
    max_iter = max_iter,
    allocation_evidence = list(
      joint_lp_problem_bytes = allocation_bytes,
      variable_count = lp_problem$variable_count,
      constraint_count = lp_problem$constraint_count
    )
  )
  class(out) <- c("rfugw_barycenter_result", "list")
  out
}

.barycenter_state_matches <- function(x, y) {
  isTRUE(all.equal(x, y, tolerance = 0, check.attributes = FALSE))
}

.validate_barycenter_state <- function(
    init_state, problem, support_cost, epsilon, min_weight) {
  if (is.null(init_state)) return(NULL)
  if (inherits(init_state, "rfugw_barycenter_result")) {
    if (!isTRUE(init_state$converged)) {
      stop("`init_state` result must be certified and converged.", call. = FALSE)
    }
    init_state <- init_state$warm_state
  }
  if (!is.list(init_state) || !isTRUE(init_state$certified) ||
      !identical(init_state$state_type, "regularized_barycenter_duals_v1")) {
    stop("`init_state` must be certified regularized barycenter state.",
         call. = FALSE)
  }
  required <- c(
    "weights", "cross_states", "self_state", "costs", "measures",
    "coefficients", "support_cost", "epsilon", "min_weight"
  )
  if (!all(required %in% names(init_state))) {
    stop("`init_state` is missing required barycenter fields.", call. = FALSE)
  }
  comparisons <- list(
    costs = .barycenter_state_matches(init_state$costs, problem$costs),
    measures = .barycenter_state_matches(init_state$measures, problem$measures),
    coefficients = .barycenter_state_matches(
      init_state$coefficients, problem$coefficients
    ),
    support_cost = .barycenter_state_matches(
      init_state$support_cost, support_cost
    ),
    epsilon = .barycenter_state_matches(init_state$epsilon, epsilon),
    min_weight = .barycenter_state_matches(init_state$min_weight, min_weight)
  )
  failed <- names(comparisons)[!vapply(comparisons, isTRUE, logical(1))]
  if (length(failed)) {
    stop(
      sprintf("`init_state` does not match the barycenter problem: %s.",
              paste(failed, collapse = ", ")),
      call. = FALSE
    )
  }
  if (!is.numeric(init_state$weights) ||
      length(init_state$weights) != problem$support_size ||
      any(!is.finite(init_state$weights)) ||
      any(init_state$weights < min_weight) ||
      abs(sum(init_state$weights) - 1) > 1e-12) {
    stop("`init_state$weights` is not feasible for the barycenter simplex.",
         call. = FALSE)
  }
  if (!is.list(init_state$cross_states) ||
      length(init_state$cross_states) != length(problem$costs) ||
      !is.list(init_state$self_state)) {
    stop("`init_state` has incompatible component state.", call. = FALSE)
  }
  init_state
}

.evaluate_regularized_barycenter <- function(
    problem, support_cost, weights, epsilon, sinkhorn_method,
    sinkhorn_max_iter, sinkhorn_tol, states) {
  cross_results <- vector("list", length(problem$costs))
  for (s in seq_along(problem$costs)) {
    cross_results[[s]] <- ot_sinkhorn(
      problem$costs[[s]], problem$measures[[s]], weights,
      epsilon = epsilon,
      method = sinkhorn_method,
      max_iter = sinkhorn_max_iter,
      tol = sinkhorn_tol,
      init_duals = states$cross[[s]] %||% NULL
    )
  }
  self_result <- ot_sinkhorn(
    support_cost, weights, weights,
    epsilon = epsilon,
    method = sinkhorn_method,
    max_iter = sinkhorn_max_iter,
    tol = sinkhorn_tol,
    init_duals = states$self %||% NULL
  )
  all_results <- c(cross_results, list(self_result))
  certified <- all(vapply(
    all_results,
    function(x) isTRUE(x$converged) && isTRUE(x$feasible) &&
      isTRUE(x$objective_consistent) &&
      isTRUE(x$regularized_dual_consistent),
    logical(1)
  ))
  cross_values <- vapply(
    cross_results, function(x) x$regularized_objective, numeric(1)
  )
  objective <- sum(problem$coefficients * cross_values) -
    0.5 * self_result$regularized_objective
  gradient <- Reduce(`+`, Map(
    function(result, coefficient) {
      coefficient * result$regularized_target_potential
    },
    cross_results, problem$coefficients
  )) - 0.5 * (
    self_result$regularized_source_potential +
      self_result$regularized_target_potential
  )
  list(
    objective = objective,
    gradient = unname(as.numeric(gradient)),
    cross_results = cross_results,
    self_result = self_result,
    cross_values = cross_values,
    certified = certified,
    states = list(
      cross = lapply(cross_results, function(x) x$dual_state),
      self = self_result$dual_state
    ),
    component_iterations = sum(vapply(
      all_results, function(x) x$iterations, integer(1)
    )),
    warm_accepted = sum(vapply(
      all_results,
      function(x) !identical(x$initialization, "cold"),
      logical(1)
    )),
    solves = length(all_results)
  )
}

.solve_regularized_barycenter <- function(
    problem, support_cost, epsilon, init_weights, init_state, min_weight,
    max_iter, tol, step_size, backtracking, armijo, sinkhorn_method,
    sinkhorn_max_iter, sinkhorn_tol) {
  support_cost <- .validate_finite_matrix(
    support_cost, "support_cost", square = TRUE
  )
  if (nrow(support_cost) != problem$support_size) {
    stop("`support_cost` must match the common target-support size.",
         call. = FALSE)
  }
  if (any(support_cost < 0)) {
    stop("`support_cost` must be nonnegative.", call. = FALSE)
  }
  if (!isTRUE(all.equal(
    support_cost, t(support_cost), tolerance = 1e-12,
    check.attributes = FALSE
  ))) {
    stop("`support_cost` must be symmetric for the debiased objective.",
         call. = FALSE)
  }
  epsilon <- .validate_positive_scalar(epsilon, "epsilon")
  min_weight <- .validate_nonneg_scalar(min_weight, "min_weight")
  if (min_weight * problem$support_size >= 1) {
    stop("`min_weight` must be smaller than 1 / support size.", call. = FALSE)
  }
  max_iter <- .validate_count(max_iter, "max_iter")
  tol <- .validate_positive_scalar(tol, "tol")
  step_size <- .validate_positive_scalar(step_size, "step_size")
  backtracking <- .validate_count(backtracking, "backtracking")
  armijo <- .validate_positive_scalar(armijo, "armijo")
  if (armijo >= 1) stop("`armijo` must be smaller than one.", call. = FALSE)
  sinkhorn_method <- match.arg(sinkhorn_method, c("auto", "scaling", "log"))
  sinkhorn_max_iter <- .validate_count(sinkhorn_max_iter, "sinkhorn_max_iter")
  sinkhorn_tol <- .validate_positive_scalar(sinkhorn_tol, "sinkhorn_tol")
  if (!is.null(init_weights) && !is.null(init_state)) {
    stop("Supply only one of `init_weights` and `init_state`.", call. = FALSE)
  }
  state <- .validate_barycenter_state(
    init_state, problem, support_cost, epsilon, min_weight
  )
  if (!is.null(state)) {
    weights <- state$weights
    states <- list(cross = state$cross_states, self = state$self_state)
  } else {
    if (is.null(init_weights)) {
      init_weights <- rep(1 / problem$support_size, problem$support_size)
    }
    if (!is.numeric(init_weights) ||
        length(init_weights) != problem$support_size ||
        any(!is.finite(init_weights)) || any(init_weights < 0) ||
        sum(init_weights) <= 0) {
      stop(
        "`init_weights` must be finite, nonnegative, and match the support size.",
        call. = FALSE
      )
    }
    weights <- .project_simplex_lower(init_weights / sum(init_weights), min_weight)
    states <- list(
      cross = vector("list", length(problem$costs)),
      self = NULL
    )
  }

  total_solves <- 0L
  total_component_iterations <- 0L
  warm_accepted <- 0L
  evaluate <- function(candidate, candidate_states) {
    value <- .evaluate_regularized_barycenter(
      problem, support_cost, candidate, epsilon, sinkhorn_method,
      sinkhorn_max_iter, sinkhorn_tol, candidate_states
    )
    total_solves <<- total_solves + value$solves
    total_component_iterations <<-
      total_component_iterations + value$component_iterations
    warm_accepted <<- warm_accepted + value$warm_accepted
    value
  }

  current <- evaluate(weights, states)
  trace <- data.frame(
    iteration = 0L,
    objective = current$objective,
    kkt_residual = NA_real_,
    step_size = NA_real_,
    component_solves = total_solves,
    component_iterations = total_component_iterations
  )
  status <- if (current$certified) "max_iter" else "inner_failure"
  termination_reason <- if (current$certified) {
    "maximum_iterations"
  } else {
    "uncertified_initial_component"
  }
  iterations <- 0L
  current_step <- step_size

  if (current$certified) {
    for (iteration in seq_len(max_iter)) {
      projected <- .project_simplex_lower(
        weights - current$gradient, min_weight
      )
      kkt_residual <- max(abs(weights - projected))
      trace$kkt_residual[nrow(trace)] <- kkt_residual
      if (is.finite(kkt_residual) && kkt_residual <= tol) {
        status <- "converged"
        termination_reason <- "projected_gradient_kkt_certified"
        break
      }

      accepted <- FALSE
      uncertified_candidates <- 0L
      local_step <- current_step
      for (line_search in seq_len(backtracking)) {
        candidate_weights <- .project_simplex_lower(
          weights - local_step * current$gradient, min_weight
        )
        direction <- candidate_weights - weights
        candidate <- evaluate(candidate_weights, current$states)
        if (!candidate$certified) {
          uncertified_candidates <- uncertified_candidates + 1L
          local_step <- local_step / 2
          next
        }
        armijo_bound <- current$objective +
          armijo * sum(current$gradient * direction)
        if (candidate$objective <= armijo_bound +
            10 * .Machine$double.eps * (1 + abs(current$objective))) {
          accepted <- TRUE
          break
        }
        local_step <- local_step / 2
      }
      if (!accepted) {
        status <- if (uncertified_candidates == backtracking) {
          "inner_failure"
        } else {
          "line_search_failure"
        }
        termination_reason <- if (identical(status, "inner_failure")) {
          "all_line_search_components_uncertified"
        } else {
          "armijo_backtracking_exhausted"
        }
        break
      }

      weights <- candidate_weights
      current <- candidate
      current_step <- min(step_size, 2 * local_step)
      iterations <- iteration
      projected <- .project_simplex_lower(
        weights - current$gradient, min_weight
      )
      kkt_residual <- max(abs(weights - projected))
      trace <- rbind(trace, data.frame(
        iteration = iteration,
        objective = current$objective,
        kkt_residual = kkt_residual,
        step_size = local_step,
        component_solves = total_solves,
        component_iterations = total_component_iterations
      ))
      if (kkt_residual <= tol) {
        status <- "converged"
        termination_reason <- "projected_gradient_kkt_certified"
        break
      }
    }
  }

  final_projected <- .project_simplex_lower(
    weights - current$gradient, min_weight
  )
  kkt_residual <- max(abs(weights - final_projected))
  kkt_tolerance <- tol
  simplex_residual <- max(
    abs(sum(weights) - 1), max(c(min_weight - weights, 0))
  )
  objective_monotone <- all(diff(trace$objective) <=
    1e-10 * (1 + abs(head(trace$objective, -1L))))
  objective_recomputed <- sum(
    problem$coefficients * vapply(
      current$cross_results,
      function(x) x$regularized_objective,
      numeric(1)
    )
  ) - 0.5 * current$self_result$regularized_objective
  objective_residual <- abs(current$objective - objective_recomputed)
  objective_tolerance <- max(
    1e-10,
    20 * sinkhorn_tol * (1 + abs(current$objective))
  )
  objective_consistent <- is.finite(objective_residual) &&
    objective_residual <= objective_tolerance
  feasible <- simplex_residual <= max(1e-12, 10 * tol)
  kkt_consistent <- kkt_residual <= kkt_tolerance
  converged <- identical(status, "converged") && current$certified &&
    feasible && objective_consistent && objective_monotone && kkt_consistent
  if (identical(status, "converged") && !converged) {
    status <- if (!feasible) {
      "infeasible"
    } else if (!objective_consistent || !objective_monotone) {
      "objective_mismatch"
    } else {
      "stationarity_failure"
    }
    termination_reason <- status
  }
  final_states <- current$states
  warm_state <- list(
    certified = converged,
    state_type = "regularized_barycenter_duals_v1",
    weights = weights,
    cross_states = final_states$cross,
    self_state = final_states$self,
    costs = problem$costs,
    measures = problem$measures,
    coefficients = problem$coefficients,
    support_cost = support_cost,
    epsilon = epsilon,
    min_weight = min_weight
  )
  out <- list(
    weights = weights,
    barycenter_weights = weights,
    source_measures = problem$measures,
    coefficients = problem$coefficients,
    costs = problem$costs,
    support_cost = support_cost,
    component_results = current$cross_results,
    support_self_result = current$self_result,
    component_values = current$cross_values,
    component_status = vapply(
      current$cross_results, rfugw_status, character(1)
    ),
    component_residuals = vapply(
      current$cross_results, function(x) x$residual, numeric(1)
    ),
    component_solves = total_solves,
    total_component_iterations = total_component_iterations,
    warm_start_reuse_count = warm_accepted,
    barycenter_objective = current$objective,
    objective_components = problem$coefficients * current$cross_values,
    support_self_correction = -0.5 * current$self_result$regularized_objective,
    objective_recomputed = objective_recomputed,
    objective_residual = objective_residual,
    objective_tolerance = objective_tolerance,
    objective_consistent = objective_consistent,
    objective_trace = trace,
    objective_monotone = objective_monotone,
    gradient = current$gradient,
    projected_gradient = weights - final_projected,
    kkt_residual = kkt_residual,
    kkt_tolerance = kkt_tolerance,
    kkt_consistent = kkt_consistent,
    simplex_residual = simplex_residual,
    feasibility_residual = simplex_residual,
    feasibility_tolerance = max(1e-12, 10 * tol),
    feasible = feasible,
    min_weight = min_weight,
    warm_started = !is.null(state),
    warm_start_accepted = !is.null(state),
    warm_state = warm_state,
    status = status,
    converged = converged,
    termination_reason = termination_reason,
    warning_payload = if (converged) NULL else list(
      code = status,
      message = sprintf("Regularized barycenter is not certified: %s.", status)
    ),
    formulation = "fixed_support_sinkhorn_divergence_barycenter",
    mode = "regularized",
    objective_convention = paste0(
      "weighted_cross_product_kl_minus_half_support_self_product_kl; ",
      "source_self_constants_omitted"
    ),
    objective_units = "cost_times_probability_mass",
    regularization = epsilon,
    backend = paste0("projected_gradient_", current$self_result$backend),
    iterations = iterations,
    max_iter = max_iter,
    sinkhorn_method = sinkhorn_method,
    sinkhorn_max_iter = sinkhorn_max_iter,
    sinkhorn_tol = sinkhorn_tol
  )
  class(out) <- c("rfugw_barycenter_result", "list")
  out
}

#' Fixed-support Wasserstein barycenter weights
#'
#' Optimizes one probability vector on a user-supplied common target support.
#' Every `costs[[s]]` has source measure `measures[[s]]` on its rows and the
#' common barycenter support on its columns. This is a linear-Wasserstein
#' estimand, not a GW/FGW support-learning update.
#'
#' Exact mode solves one joint linear program and then independently certifies
#' every fixed-weight component with [ot_emd()]. Regularized mode minimizes a
#' semi-debiased Sinkhorn objective: the weighted cross product-reference-KL
#' objectives minus half the common-support self objective. Source-self terms
#' are constant in the barycenter weights and intentionally omitted. A
#' projected-gradient KKT map, safeguarded Armijo descent, and every component
#' certificate are required for convergence.
#'
#' @param costs Nonempty list of finite nonnegative source-by-common-support
#'   cost matrices. All matrices must have the same column count.
#' @param measures List of source probability weights, one per cost matrix.
#'   Finite nonnegative vectors with positive mass are normalized separately.
#' @param coefficients Optional nonnegative barycenter coefficients, normalized
#'   to sum one. Zero coefficients are allowed.
#' @param mode `"exact"` or `"regularized"`.
#' @param support_cost In regularized mode, the finite nonnegative symmetric
#'   cost matrix within the common support. It supplies the debiasing self term.
#' @param epsilon Positive regularization in cost units for regularized mode.
#' @param init_weights Optional deterministic starting weights. Zeros are
#'   accepted and projected onto the declared numerical interior.
#' @param init_state Optional certified `warm_state` or prior converged
#'   regularized barycenter result from the identical problem.
#' @param min_weight Numerical lower bound for regularized weights. Exact mode
#'   has no such bound and can return exact zeros.
#' @param max_iter Maximum joint-simplex iterations in regularized mode and
#'   component transport-simplex iterations in exact mode. `NULL` selects 500
#'   and 20,000, respectively.
#' @param tol KKT and simplex tolerance. `NULL` selects `1e-6` in regularized
#'   mode and `1e-10` in exact mode.
#' @param step_size Initial projected-gradient step size.
#' @param backtracking Maximum Armijo reductions per outer iteration.
#' @param armijo Armijo sufficient-decrease coefficient in `(0, 1)`.
#' @param sinkhorn_method `"auto"`, `"scaling"`, or `"log"`.
#' @param sinkhorn_max_iter Maximum iterations for each component Sinkhorn solve.
#' @param sinkhorn_tol Marginal tolerance for each component Sinkhorn solve.
#' @return An `rfugw_barycenter_result` with `weights`, component results,
#'   objective convention/components, simplex and KKT certificates, complete
#'   solve counts, warm state, status, and provenance.
#' @examples
#' support <- matrix(1:3, ncol = 1)
#' C <- outer(1:3, 1:3, function(x, y) (x - y)^2)
#' C <- C / max(C)
#' out <- ot_barycenter_weights(
#'   list(C, C),
#'   list(c(0.2, 0.5, 0.3), c(0.4, 0.2, 0.4)),
#'   mode = "regularized",
#'   support_cost = C,
#'   epsilon = 0.1
#' )
#' out$weights
#' @export
ot_barycenter_weights <- function(
    costs,
    measures,
    coefficients = NULL,
    mode = c("exact", "regularized"),
    support_cost = NULL,
    epsilon = 0.05,
    init_weights = NULL,
    init_state = NULL,
    min_weight = 1e-12,
    max_iter = NULL,
    tol = NULL,
    step_size = 1,
    backtracking = 30L,
    armijo = 1e-4,
    sinkhorn_method = c("auto", "scaling", "log"),
    sinkhorn_max_iter = 5000L,
    sinkhorn_tol = 1e-9) {
  mode <- match.arg(mode)
  problem <- .prepare_barycenter_problem(costs, measures, coefficients)
  if (identical(mode, "exact")) {
    if (is.null(max_iter)) max_iter <- 20000L
    if (is.null(tol)) tol <- 1e-10
    if (!is.null(support_cost)) {
      stop("`support_cost` is used only in `mode = \"regularized\"`.",
           call. = FALSE)
    }
    if (!is.null(init_state)) {
      stop("Exact joint-LP mode does not accept iterative `init_state`.",
           call. = FALSE)
    }
    if (!is.null(init_weights)) {
      stop("Exact joint-LP mode does not use `init_weights`.", call. = FALSE)
    }
    return(.solve_exact_barycenter(problem, max_iter, tol))
  }
  if (is.null(support_cost)) {
    stop("`support_cost` is required in `mode = \"regularized\"`.",
         call. = FALSE)
  }
  if (is.null(max_iter)) max_iter <- 500L
  if (is.null(tol)) tol <- 1e-6
  .solve_regularized_barycenter(
    problem = problem,
    support_cost = support_cost,
    epsilon = epsilon,
    init_weights = init_weights,
    init_state = init_state,
    min_weight = min_weight,
    max_iter = max_iter,
    tol = tol,
    step_size = step_size,
    backtracking = backtracking,
    armijo = armijo,
    sinkhorn_method = match.arg(sinkhorn_method),
    sinkhorn_max_iter = sinkhorn_max_iter,
    sinkhorn_tol = sinkhorn_tol
  )
}

#' Print a fixed-support Wasserstein barycenter result
#'
#' @param x An `rfugw_barycenter_result`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.rfugw_barycenter_result <- function(x, ...) {
  cat("<rfugw_barycenter_result>\n")
  cat(sprintf("  formulation: %s\n", x$formulation))
  cat(sprintf("  status:      %s\n", x$status))
  cat(sprintf("  support:     %d weights\n", length(x$weights)))
  cat(sprintf("  value:       %s\n", format(x$barycenter_objective, digits = 6)))
  cat(sprintf("  KKT residual:%s\n", format(x$kkt_residual, digits = 4)))
  cat(sprintf("  solves:      %s\n", x$component_solves))
  invisible(x)
}
