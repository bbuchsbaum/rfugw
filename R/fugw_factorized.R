.factorized_finalize_plan_stats <- function(stats, state) {
  stats$source_marginal <- as.numeric(stats$source_marginal)
  stats$target_marginal <- as.numeric(stats$target_marginal)
  source <- state$source_measure
  target <- state$target_measure
  if (any(stats$source_marginal > 0 & source == 0) ||
      any(stats$target_marginal > 0 & target == 0)) {
    stats$plan_product_kl <- Inf
  } else {
    source_active <- stats$source_marginal > 0
    target_active <- stats$target_marginal > 0
    stats$plan_product_kl <- stats$entropy_sum -
      sum(stats$source_marginal[source_active] * log(source[source_active])) -
      sum(stats$target_marginal[target_active] * log(target[target_active])) -
      stats$mass + sum(source) * sum(target)
  }
  stats
}

.factorized_plan_stats <- function(state) {
  terms <- state$cost
  stats <- cpp_factorized_plan_stats(
    terms$row,
    terms$column,
    terms$left,
    terms$right,
    state$source_measure,
    state$target_measure,
    state$source_bar,
    state$target_bar,
    state$epsilon,
    state$block_size
  )
  .factorized_finalize_plan_stats(stats, state)
}

.factorized_plan_stats_moments <- function(
    state, source_right, target_right) {
  terms <- state$cost
  native <- cpp_factorized_plan_stats_moments(
    terms$row,
    terms$column,
    terms$left,
    terms$right,
    state$source_measure,
    state$target_measure,
    state$source_bar,
    state$target_bar,
    state$epsilon,
    source_right,
    target_right,
    state$block_size
  )
  cross <- native$cross_moment
  workspace_elements <- native$workspace_elements
  native$cross_moment <- NULL
  native$workspace_elements <- NULL
  stats <- .factorized_finalize_plan_stats(native, state)
  source_gram <- if (!ncol(source_right)) {
    matrix(0, 0L, 0L)
  } else {
    crossprod(source_right, source_right * stats$source_marginal)
  }
  target_gram <- if (!ncol(target_right)) {
    matrix(0, 0L, 0L)
  } else {
    crossprod(target_right, target_right * stats$target_marginal)
  }
  list(
    stats = stats,
    moments = list(
      source_mass = stats$source_marginal,
      target_mass = stats$target_marginal,
      H = unname(cross),
      G_source = unname(source_gram),
      G_target = unname(target_gram),
      transported_mass = stats$mass
    ),
    workspace_elements = workspace_elements
  )
}

.factorized_rescale_moments <- function(moments, scale) {
  if (identical(scale, 1) || isTRUE(all.equal(scale, 1))) return(moments)
  moments$source_mass <- scale * moments$source_mass
  moments$target_mass <- scale * moments$target_mass
  moments$H <- scale * moments$H
  moments$G_source <- scale * moments$G_source
  moments$G_target <- scale * moments$G_target
  moments$transported_mass <- scale * moments$transported_mass
  moments
}

.factorized_plan_apply_state <- function(state, values, adjoint = FALSE) {
  terms <- state$cost
  cpp_factorized_plan_apply(
    terms$row,
    terms$column,
    terms$left,
    terms$right,
    state$source_measure,
    state$target_measure,
    state$source_bar,
    state$target_bar,
    state$epsilon,
    values,
    adjoint,
    state$block_size
  )
}

.factorized_plan_materialize_state <- function(state) {
  terms <- state$cost
  cpp_factorized_plan_materialize(
    terms$row,
    terms$column,
    terms$left,
    terms$right,
    state$source_measure,
    state$target_measure,
    state$source_bar,
    state$target_bar,
    state$epsilon,
    state$block_size
  )
}

.factorized_plan_from_state <- function(state) {
  stats <- state$stats
  plan <- transport_operator(
    length(state$source_measure),
    length(state$target_measure),
    apply = function(x) .factorized_plan_apply_state(state, x, FALSE),
    adjoint = function(x) .factorized_plan_apply_state(state, x, TRUE),
    source_mass = stats$source_marginal,
    target_mass = stats$target_marginal,
    materialize = function() .factorized_plan_materialize_state(state),
    metadata = list(
      formulation = state$formulation,
      cost_representation = "affine_bilinear",
      coupling_stored = FALSE,
      full_support = TRUE,
      factor_rank = ncol(state$cost$left),
      block_size = state$block_size,
      max_tile_elements = state$max_tile_elements,
      full_matrix_elements = state$full_matrix_elements,
      exact_operator_actions = TRUE
    )
  )
  plan$data$factorized_state <- state
  plan
}

.factorized_plan_state <- function(plan, name = "plan") {
  if (!.is_transport_plan(plan) ||
      !identical(plan$representation, "implicit_operator") ||
      is.null(plan$data$factorized_state)) {
    stop(sprintf("`%s` is not an implicit factorized rfugw plan.", name),
         call. = FALSE)
  }
  plan$data$factorized_state
}

.factorized_product_plan <- function(source, target, block_size) {
  n_source <- length(source)
  n_target <- length(target)
  source_active <- source > 0
  target_active <- target > 0
  entropy_sum <- sum(source[source_active] * log(source[source_active])) +
    sum(target[target_active] * log(target[target_active]))
  state <- list(
    formulation = "factorized_product_initialization",
    cost = list(
      row = numeric(n_source),
      column = numeric(n_target),
      left = matrix(numeric(), nrow = n_source, ncol = 0L),
      right = matrix(numeric(), nrow = n_target, ncol = 0L)
    ),
    source_measure = source,
    target_measure = target,
    source_bar = numeric(n_source),
    target_bar = numeric(n_target),
    epsilon = 1,
    block_size = block_size,
    max_tile_elements = min(n_source, block_size) * min(n_target, block_size),
    full_matrix_elements = as.double(n_source) * as.double(n_target),
    stats = list(
      source_marginal = source,
      target_marginal = target,
      mass = 1,
      transport_cost = 0,
      entropy_sum = entropy_sum,
      plan_product_kl = 0,
      minimum_log_weight = min(
        outer(log(source[source_active]), log(target[target_active]), "+")
      ),
      maximum_log_weight = max(
        outer(log(source[source_active]), log(target[target_active]), "+")
      ),
      finite_plan = TRUE
    )
  )
  .factorized_plan_from_state(state)
}

.factorized_rescale_plan <- function(plan, scale) {
  if (!is.numeric(scale) || length(scale) != 1L ||
      !is.finite(scale) || scale <= 0) {
    stop("Implicit plan rescaling requires a finite positive scale.",
         call. = FALSE)
  }
  state <- .factorized_plan_state(plan)
  log_scale <- log(scale)
  shift <- 0.5 * state$epsilon * log_scale
  state$source_bar <- state$source_bar + shift
  state$target_bar <- state$target_bar + shift
  old_mass <- state$stats$mass
  state$stats$source_marginal <- scale * state$stats$source_marginal
  state$stats$target_marginal <- scale * state$stats$target_marginal
  state$stats$mass <- scale * old_mass
  state$stats$transport_cost <- scale * state$stats$transport_cost
  state$stats$entropy_sum <- scale * state$stats$entropy_sum +
    scale * old_mass * log_scale
  state$stats$minimum_log_weight <-
    state$stats$minimum_log_weight + log_scale
  state$stats$maximum_log_weight <-
    state$stats$maximum_log_weight + log_scale
  source <- state$source_measure
  target <- state$target_measure
  source_active <- state$stats$source_marginal > 0
  target_active <- state$stats$target_marginal > 0
  state$stats$plan_product_kl <- state$stats$entropy_sum -
    sum(state$stats$source_marginal[source_active] *
          log(source[source_active])) -
    sum(state$stats$target_marginal[target_active] *
          log(target[target_active])) -
    state$stats$mass + sum(source) * sum(target)
  .factorized_plan_from_state(state)
}

.factorized_plan_difference <- function(left_plan, right_plan, block_size) {
  left_state <- .factorized_plan_state(left_plan, "left_plan")
  right_state <- .factorized_plan_state(right_plan, "right_plan")
  left_cost <- left_state$cost
  right_cost <- right_state$cost
  cpp_factorized_plan_difference(
    left_cost$row,
    left_cost$column,
    left_cost$left,
    left_cost$right,
    left_state$source_measure,
    left_state$target_measure,
    left_state$source_bar,
    left_state$target_bar,
    left_state$epsilon,
    right_cost$row,
    right_cost$column,
    right_cost$left,
    right_cost$right,
    right_state$source_measure,
    right_state$target_measure,
    right_state$source_bar,
    right_state$target_bar,
    right_state$epsilon,
    block_size
  )
}

.factorized_relative_change <- function(current, previous) {
  if (!identical(dim(current), dim(previous)) ||
      length(current) != length(previous)) return(Inf)
  if (!length(current)) return(0)
  max(abs(current - previous)) /
    max(1, max(abs(current)), max(abs(previous)))
}

.factorized_plan_screening_residual <- function(current_plan, previous_plan) {
  if (is.null(current_plan) || is.null(previous_plan)) return(Inf)
  current <- .factorized_plan_state(current_plan)
  previous <- .factorized_plan_state(previous_plan)
  if (length(current$source_bar) != length(previous$source_bar) ||
      length(current$target_bar) != length(previous$target_bar)) return(Inf)
  source_delta <- current$source_bar - previous$source_bar
  target_delta <- current$target_bar - previous$target_bar
  gauge_shift <- 0.5 * (mean(source_delta) - mean(target_delta))
  potential_change <- max(
    abs(source_delta - gauge_shift),
    abs(target_delta + gauge_shift)
  ) / max(
    1,
    abs(current$source_bar), abs(previous$source_bar),
    abs(current$target_bar), abs(previous$target_bar)
  )
  cost_change <- max(
    .factorized_relative_change(current$cost$row, previous$cost$row),
    .factorized_relative_change(
      current$cost$column, previous$cost$column
    ),
    .factorized_relative_change(current$cost$left, previous$cost$left),
    .factorized_relative_change(current$cost$right, previous$cost$right)
  )
  marginal_change <- max(
    sum(abs(
      current$stats$source_marginal - previous$stats$source_marginal
    )),
    sum(abs(
      current$stats$target_marginal - previous$stats$target_marginal
    )),
    abs(current$stats$mass - previous$stats$mass)
  ) / max(1, current$stats$mass, previous$stats$mass)
  max(
    potential_change,
    cost_change,
    marginal_change,
    abs(current$epsilon - previous$epsilon) /
      max(1, current$epsilon, previous$epsilon)
  )
}

.factorized_ti_kkt_residual <- function(marginal, reference, potential, rho) {
  active <- reference > 0
  if (any(marginal[active] <= 0) || any(marginal[!active] > 0)) return(Inf)
  max(abs(
    log(marginal[active]) -
      (log(reference[active]) - potential[active] / rho)
  ))
}

.factorized_ti_certificate <- function(
    state, native, rho, tolerance) {
  stats <- state$stats
  source <- state$source_measure
  target <- state$target_measure
  source_kl <- .generalized_kl_vector(stats$source_marginal, source)
  target_kl <- .generalized_kl_vector(stats$target_marginal, target)
  source_kkt <- .factorized_ti_kkt_residual(
    stats$source_marginal, source, native$source_potential, rho[[1L]]
  )
  target_kkt <- .factorized_ti_kkt_residual(
    stats$target_marginal, target, native$target_potential, rho[[2L]]
  )
  kkt_residual <- max(source_kkt, target_kkt)
  primal <- stats$transport_cost + rho[[1L]] * source_kl +
    rho[[2L]] * target_kl + state$epsilon * stats$plan_product_kl
  dual_source_mass <- sum(source * exp(-native$source_potential / rho[[1L]]))
  dual_target_mass <- sum(target * exp(-native$target_potential / rho[[2L]]))
  dual <- rho[[1L]] * (sum(source) - dual_source_mass) +
    rho[[2L]] * (sum(target) - dual_target_mass) +
    state$epsilon * (sum(source) * sum(target) - stats$mass)
  gap <- primal - dual
  objective_tolerance <- max(1e-8, 100 * tolerance) * max(1, abs(primal))
  kkt_tolerance <- max(1e-7, 100 * tolerance)
  fixed_point_tolerance <- max(1e-8, 20 * tolerance)
  list(
    source_marginal = stats$source_marginal,
    target_marginal = stats$target_marginal,
    finite_plan = isTRUE(stats$finite_plan),
    fixed_point_residual = native$fixed_point_residual,
    fixed_point_tolerance = fixed_point_tolerance,
    fixed_point_consistent = is.finite(native$fixed_point_residual) &&
      native$fixed_point_residual <= fixed_point_tolerance,
    source_kkt_residual = source_kkt,
    target_kkt_residual = target_kkt,
    kkt_residual = kkt_residual,
    kkt_tolerance = kkt_tolerance,
    kkt_consistent = is.finite(kkt_residual) &&
      kkt_residual <= kkt_tolerance,
    transport = stats$transport_cost,
    source_kl = source_kl,
    target_kl = target_kl,
    product_kl = stats$plan_product_kl,
    primal = primal,
    dual = dual,
    gap = gap,
    objective_tolerance = objective_tolerance,
    objective_consistent = is.finite(gap) &&
      gap >= -objective_tolerance && abs(gap) <= objective_tolerance,
    transported_mass = stats$mass,
    dual_source_mass = dual_source_mass,
    dual_target_mass = dual_target_mass,
    gauge_residual = 0,
    gauge_invariant = TRUE,
    minimum_log_weight = stats$minimum_log_weight,
    maximum_log_weight = stats$maximum_log_weight
  )
}

.factorized_terms_block <- function(terms, rows, columns) {
  value <- outer(terms$row[rows], terms$column[columns], "+")
  if (ncol(terms$left)) {
    value <- value + tcrossprod(
      terms$left[rows, , drop = FALSE],
      terms$right[columns, , drop = FALSE]
    )
  }
  unname(value)
}

.factorized_block_kkt_audit <- function(
    plan,
    final_cost,
    epsilon,
    rho,
    tolerance,
    block_size,
    post_solve_scale = 1) {
  state <- .factorized_plan_state(plan)
  final_terms <- .cost_native_terms(final_cost, "final_cost")
  stored_terms <- state$cost
  source <- state$source_measure
  target <- state$target_measure
  stats <- state$stats
  block_size <- .validate_count(block_size, "block_size")
  shape <- cost_shape(final_cost)
  shape_consistent <- identical(
    unname(shape), as.integer(c(length(source), length(target)))
  )
  parameter_consistent <- is.numeric(epsilon) && length(epsilon) == 1L &&
    is.finite(epsilon) && epsilon > 0 &&
    is.numeric(rho) && length(rho) == 2L &&
    all(is.finite(rho)) && all(rho > 0)

  source_active <- source > 0
  target_active <- target > 0
  support_consistent <-
    length(stats$source_marginal) == length(source) &&
    length(stats$target_marginal) == length(target) &&
    all(is.finite(stats$source_marginal)) &&
    all(is.finite(stats$target_marginal)) &&
    all(stats$source_marginal[source_active] > 0) &&
    all(stats$target_marginal[target_active] > 0) &&
    all(stats$source_marginal[!source_active] == 0) &&
    all(stats$target_marginal[!target_active] == 0)
  state_finite <- isTRUE(stats$finite_plan) &&
    is.finite(stats$mass) && stats$mass > 0 &&
    is.finite(stats$plan_product_kl) &&
    all(is.finite(state$source_bar)) && all(is.finite(state$target_bar)) &&
    is.finite(stats$minimum_log_weight) &&
    is.finite(stats$maximum_log_weight)

  kkt_residual <- Inf
  maximum_cost_change <- Inf
  maximum_final_cost <- Inf
  maximum_stored_cost <- Inf
  tile_count <- 0L
  tiles_finite <- FALSE
  transport <- Inf
  workspace_elements <- NA_real_
  if (shape_consistent && parameter_consistent && support_consistent &&
      state_finite) {
    native_audit <- cpp_factorized_plan_kkt_audit(
      final_terms$row,
      final_terms$column,
      final_terms$left,
      final_terms$right,
      stored_terms$row,
      stored_terms$column,
      stored_terms$left,
      stored_terms$right,
      source,
      target,
      state$source_bar,
      state$target_bar,
      stats$source_marginal,
      stats$target_marginal,
      epsilon,
      rho[[1L]],
      rho[[2L]],
      block_size
    )
    tiles_finite <- isTRUE(native_audit$finite)
    kkt_residual <- native_audit$kkt_residual
    maximum_cost_change <- native_audit$maximum_cost_change
    maximum_final_cost <- native_audit$maximum_final_cost
    maximum_stored_cost <- native_audit$maximum_stored_cost
    transport <- native_audit$transport
    tile_count <- as.integer(native_audit$tile_count)
    workspace_elements <- native_audit$workspace_elements
  }

  source_kl <- tryCatch(
    .generalized_kl_vector(stats$source_marginal, source),
    error = function(...) Inf
  )
  target_kl <- tryCatch(
    .generalized_kl_vector(stats$target_marginal, target),
    error = function(...) Inf
  )
  block_objective <- transport + rho[[1L]] * source_kl +
    rho[[2L]] * target_kl + epsilon * stats$plan_product_kl
  objective_finite <- all(is.finite(c(
    transport, source_kl, target_kl, stats$plan_product_kl, block_objective
  )))
  scale <- max(1, epsilon, rho)
  normalized_kkt_residual <- kkt_residual / scale
  normalized_tolerance <- max(1e-7, 100 * tolerance)
  roundoff_tolerance <- 1000 * .Machine$double.eps * max(
    1,
    maximum_final_cost,
    maximum_stored_cost,
    abs(state$source_bar),
    abs(state$target_bar),
    rho
  )
  effective_tolerance <- max(
    normalized_tolerance * scale,
    roundoff_tolerance
  )
  numerical_finite <- shape_consistent && parameter_consistent &&
    support_consistent && state_finite && tiles_finite && objective_finite &&
    is.finite(kkt_residual) && is.finite(effective_tolerance)
  certified <- numerical_finite && kkt_residual <= effective_tolerance
  underflow_log_threshold <- log(.Machine$double.xmin) +
    log(.Machine$double.eps)

  list(
    certified = certified,
    method = "post_rescale_full_support_primal_kkt",
    kkt_residual = kkt_residual,
    normalized_kkt_residual = normalized_kkt_residual,
    effective_tolerance = effective_tolerance,
    normalized_tolerance = normalized_tolerance,
    transport = transport,
    source_kl = source_kl,
    target_kl = target_kl,
    product_kl = stats$plan_product_kl,
    block_objective = block_objective,
    objective_finite = objective_finite,
    numerical_finite = numerical_finite,
    support_consistent = support_consistent,
    maximum_cost_change = maximum_cost_change,
    maximum_final_cost = maximum_final_cost,
    maximum_stored_cost = maximum_stored_cost,
    post_solve_scale = post_solve_scale,
    transported_mass = stats$mass,
    minimum_log_weight = stats$minimum_log_weight,
    maximum_log_weight = stats$maximum_log_weight,
    active_entry_underflow_possible =
      stats$minimum_log_weight < underflow_log_threshold,
    active_marginals_representable = support_consistent && stats$mass > 0,
    tile_count = tile_count,
    workspace_elements = workspace_elements,
    workspace_bound = paste0(
      "2 * min(n_source, block_size) * min(n_target, block_size) + ",
      "2 * (n_source + n_target) doubles"
    ),
    workspace_scope = paste0(
      "Peak live native doubles for two reusable cost tiles, log-measure ",
      "vectors, and marginal KKT base vectors."
    ),
    max_tile_elements = min(sum(source_active), block_size) *
      min(sum(target_active), block_size),
    full_matrix_elements = as.double(length(source)) * as.double(length(target)),
    interpretation = paste0(
      "Entrywise first-order residual for the returned implicit coupling ",
      "under the cost induced by the other returned coupling. Individual ",
      "negligible entries may underflow, but every active marginal and total ",
      "mass must remain representable."
    )
  )
}

.fugw_final_block_audits <- function(
    sample_plan,
    feature_plan,
    sample_moments,
    feature_moments,
    source_left,
    target_left,
    feature_terms,
    source,
    target,
    rho,
    epsilon,
    feature_weight,
    structure_weight,
    tolerance,
    block_size,
    sample_scale = 1,
    feature_scale = 1,
    exact = TRUE) {
  sample_state <- .factorized_plan_state(sample_plan)
  feature_state <- .factorized_plan_state(feature_plan)
  feature_cost <- .fugw_dynamic_cost(
    sample_moments,
    sample_state$stats,
    source_left,
    target_left,
    feature_terms,
    source,
    target,
    rho,
    epsilon,
    feature_weight,
    structure_weight,
    list(exact = exact, role = "final_feature_coupling_audit")
  )
  sample_cost <- .fugw_dynamic_cost(
    feature_moments,
    feature_state$stats,
    source_left,
    target_left,
    feature_terms,
    source,
    target,
    rho,
    epsilon,
    feature_weight,
    structure_weight,
    list(exact = exact, role = "final_sample_coupling_audit")
  )
  feature_audit <- .factorized_block_kkt_audit(
    feature_plan,
    feature_cost,
    epsilon * sample_state$stats$mass,
    rho * sample_state$stats$mass,
    tolerance,
    block_size,
    post_solve_scale = feature_scale
  )
  sample_audit <- .factorized_block_kkt_audit(
    sample_plan,
    sample_cost,
    epsilon * feature_state$stats$mass,
    rho * feature_state$stats$mass,
    tolerance,
    block_size,
    post_solve_scale = sample_scale
  )
  list(
    certified = isTRUE(sample_audit$certified) &&
      isTRUE(feature_audit$certified),
    sample = sample_audit,
    feature = feature_audit,
    residual = max(
      sample_audit$normalized_kkt_residual,
      feature_audit$normalized_kkt_residual
    ),
    tolerance = max(
      sample_audit$normalized_tolerance,
      feature_audit$normalized_tolerance
    ),
    interpretation = paste0(
      "Both returned couplings are audited against the dynamic costs induced ",
      "by the other returned coupling after plan rescaling."
    )
  )
}

.validate_factorized_potential_init <- function(init, n_source, n_target) {
  if (is.null(init)) {
    return(list(source_bar = numeric(), target_bar = numeric()))
  }
  if (!is.list(init) ||
      !all(c("source_bar", "target_bar") %in% names(init))) {
    stop("`init_potentials` must contain `source_bar` and `target_bar`.",
         call. = FALSE)
  }
  source_bar <- init$source_bar
  target_bar <- init$target_bar
  if (!is.numeric(source_bar) || length(source_bar) != n_source ||
      any(!is.finite(source_bar)) ||
      !is.numeric(target_bar) || length(target_bar) != n_target ||
      any(!is.finite(target_bar))) {
    stop("Initial potentials must be finite and match the cost shape.",
         call. = FALSE)
  }
  list(source_bar = as.numeric(source_bar), target_bar = as.numeric(target_bar))
}

.ot_sinkhorn_unbalanced_ti_factorized <- function(
    cost,
    p,
    q,
    epsilon,
    rho,
    max_iter,
    tol,
    plan,
    init_potentials = NULL,
    block_size = 256L,
    moment_source = NULL,
    moment_target = NULL) {
  shape <- cost_shape(cost)
  terms <- .cost_native_terms(cost)
  p <- .validate_finite_measure(
    p, shape[["source"]], "p", rep(1 / shape[["source"]], shape[["source"]])
  )
  q <- .validate_finite_measure(
    q, shape[["target"]], "q", rep(1 / shape[["target"]], shape[["target"]])
  )
  if (sum(p) <= 0 || sum(q) <= 0) {
    stop("Factorized TI-UOT requires positive total mass on both measures.",
         call. = FALSE)
  }
  if (identical(plan, "sparse")) {
    stop(
      "A full-support factorized cost cannot return a sparse plan without an explicit support-refinement contract; use `plan = \"operator\"` or `\"dense\"`.",
      call. = FALSE
    )
  }
  block_size <- .validate_count(block_size, "block_size")
  init <- .validate_factorized_potential_init(
    init_potentials, shape[["source"]], shape[["target"]]
  )
  combined_moments <- !is.null(moment_source) || !is.null(moment_target)
  if (combined_moments) {
    if (!is.matrix(moment_source) || !is.numeric(moment_source) ||
        nrow(moment_source) != shape[["source"]] ||
        any(!is.finite(moment_source)) ||
        !is.matrix(moment_target) || !is.numeric(moment_target) ||
        nrow(moment_target) != shape[["target"]] ||
        any(!is.finite(moment_target))) {
      stop(
        "Combined moment factors must be finite matrices matching the cost shape.",
        call. = FALSE
      )
    }
    moment_source <- unname(moment_source)
    moment_target <- unname(moment_target)
  }
  native <- cpp_ot_sinkhorn_unbalanced_ti_factorized(
    terms$row,
    terms$column,
    terms$left,
    terms$right,
    p,
    q,
    epsilon,
    rho[[1L]],
    rho[[2L]],
    max_iter,
    tol,
    block_size,
    init$source_bar,
    init$target_bar
  )
  state <- list(
    formulation = "translation_invariant_kl_uot_factorized",
    cost = terms,
    cost_provenance = cost_provenance(cost),
    source_measure = p,
    target_measure = q,
    source_bar = as.numeric(native$source_bar),
    target_bar = as.numeric(native$target_bar),
    source_potential = as.numeric(native$source_potential),
    target_potential = as.numeric(native$target_potential),
    epsilon = epsilon,
    rho = rho,
    block_size = block_size,
    max_tile_elements = native$max_tile_elements,
    full_matrix_elements = native$full_matrix_elements
  )
  combined <- if (combined_moments) {
    .factorized_plan_stats_moments(
      state, moment_source, moment_target
    )
  } else NULL
  state$stats <- if (is.null(combined)) {
    .factorized_plan_stats(state)
  } else combined$stats
  certificate <- .factorized_ti_certificate(state, native, rho, tol)
  implicit_plan <- .factorized_plan_from_state(state)
  result_plan <- if (identical(plan, "dense")) {
    transport_plan_materialize(implicit_plan)
  } else {
    implicit_plan
  }
  inner_residual <- max(
    certificate$fixed_point_residual,
    certificate$kkt_residual
  )
  inner_tolerance <- max(
    certificate$fixed_point_tolerance,
    certificate$kkt_tolerance
  )
  inner_converged <- isTRUE(native$converged) &&
    isTRUE(certificate$fixed_point_consistent) &&
    isTRUE(certificate$kkt_consistent) &&
    isTRUE(certificate$objective_consistent) &&
    isTRUE(certificate$gauge_invariant)
  support <- list(
    uot_primal_feasible = TRUE,
    finite_potential_support = TRUE,
    uncovered_source = integer(),
    uncovered_target = integer(),
    active_support_size = as.double(shape[[1L]]) * as.double(shape[[2L]]),
    balanced_marginal_feasible = abs(sum(p) - sum(q)) <=
      max(1e-12, tol * max(1, sum(p), sum(q))),
    balanced_equal_total_mass = abs(sum(p) - sum(q)) <=
      max(1e-12, tol * max(1, sum(p), sum(q))),
    balanced_check_method = "complete_support_equal_mass",
    interpretation = paste0(
      "KL-UOT is solved over complete implicit support. No support edges or ",
      "coupling entries are stored."
    )
  )
  out <- list(
    plan = result_plan,
    implicit_plan = implicit_plan,
    implicit_state = state,
    formulation = "ot_sinkhorn_unbalanced_ti_factorized",
    backend = native$backend,
    objective = certificate$primal,
    regularized_objective = certificate$primal,
    dual_objective = certificate$dual,
    primal_dual_gap = certificate$gap,
    transport_cost = certificate$transport,
    source_marginal_kl = certificate$source_kl,
    target_marginal_kl = certificate$target_kl,
    plan_product_kl = certificate$product_kl,
    transported_mass = certificate$transported_mass,
    source_marginal = certificate$source_marginal,
    target_marginal = certificate$target_marginal,
    source_measure_effective = p,
    target_measure_effective = q,
    regularization = epsilon,
    rho = rho,
    source_bar = state$source_bar,
    target_bar = state$target_bar,
    fbar = state$source_bar,
    gbar = state$target_bar,
    source_potential = state$source_potential,
    target_potential = state$target_potential,
    f = state$source_potential,
    g = state$target_potential,
    fixed_point_residual = certificate$fixed_point_residual,
    fixed_point_tolerance = certificate$fixed_point_tolerance,
    fixed_point_consistent = certificate$fixed_point_consistent,
    source_kkt_residual = certificate$source_kkt_residual,
    target_kkt_residual = certificate$target_kkt_residual,
    kkt_residual = certificate$kkt_residual,
    kkt_tolerance = certificate$kkt_tolerance,
    kkt_consistent = certificate$kkt_consistent,
    gauge_residual = certificate$gauge_residual,
    gauge_invariant = certificate$gauge_invariant,
    support_certificate = support,
    cost_representation = "affine_bilinear_full_support",
    cost_provenance = cost_provenance(cost),
    support_size = support$active_support_size,
    support_density = 1,
    plan_representation = transport_plan_representation(result_plan),
    iterations = native$iterations,
    error = max(native$residual, inner_residual),
    native_residual = native$residual,
    inner_residual = inner_residual,
    max_inner_residual = inner_residual,
    inner_iterations = native$iterations,
    inner_converged = inner_converged,
    inner_status = if (inner_converged) "converged" else "certificate_failure",
    warm_started = !is.null(init_potentials),
    block_size = block_size,
    max_tile_elements = native$max_tile_elements,
    full_matrix_elements = native$full_matrix_elements
  )
  if (!is.null(combined)) {
    out$moments <- combined$moments
    out$moment_workspace_elements <- combined$workspace_elements
  }
  ans <- .attach_solver_diagnostics(
    out,
    residual = max(native$residual, inner_residual),
    converged = isTRUE(native$converged),
    iterations = native$iterations,
    max_iter = max_iter,
    plan = result_plan,
    inner_residual = inner_residual,
    max_inner_residual = inner_residual,
    inner_iterations = native$iterations,
    inner_converged = inner_converged,
    inner_status = out$inner_status,
    feasibility = "unbalanced",
    feasibility_tol = inner_tolerance,
    objective_recomputed = certificate$primal,
    objective_tolerance = certificate$objective_tolerance,
    objective_components = list(
      transport_cost = certificate$transport,
      source_marginal_kl = certificate$source_kl,
      target_marginal_kl = certificate$target_kl,
      plan_product_kl = certificate$product_kl
    ),
    numerical_ok = isTRUE(native$numerical_ok) &&
      isTRUE(certificate$finite_plan)
  )
  ans$termination_reason <- .termination_reason_from_result(ans, max_iter)
  ans$uot_certificate <- certificate
  ans$runtime_provenance$timing <- list(
    solve_seconds = unname(native$solve_seconds)
  )
  ans$runtime_provenance$memory_contract <- list(
    full_matrix_allocated = identical(plan, "dense"),
    max_tile_elements = native$max_tile_elements,
    full_matrix_elements = native$full_matrix_elements,
    factor_rank = ncol(terms$left)
  )
  if (!is.null(combined)) {
    ans$runtime_provenance$memory_contract$combined_stats_moments <- TRUE
    ans$runtime_provenance$memory_contract$moment_workspace_elements <-
      combined$workspace_elements
  }
  ans
}

.factorized_structure_moments <- function(plan, source_right, target_right) {
  source_mass <- transport_plan_mass(plan, "source")
  target_mass <- transport_plan_mass(plan, "target")
  cross <- if (!ncol(source_right) || !ncol(target_right)) {
    matrix(0, nrow = ncol(source_right), ncol = ncol(target_right))
  } else {
    crossprod(source_right, transport_plan_apply(plan, target_right))
  }
  source_gram <- if (!ncol(source_right)) {
    matrix(0, 0L, 0L)
  } else {
    crossprod(source_right, source_right * source_mass)
  }
  target_gram <- if (!ncol(target_right)) {
    matrix(0, 0L, 0L)
  } else {
    crossprod(target_right, target_right * target_mass)
  }
  list(
    source_mass = source_mass,
    target_mass = target_mass,
    H = unname(cross),
    G_source = unname(source_gram),
    G_target = unname(target_gram),
    transported_mass = sum(source_mass)
  )
}

.factorized_product_moments <- function(
    source, target, source_right, target_right) {
  source_mean <- if (ncol(source_right)) {
    as.numeric(crossprod(source_right, source))
  } else {
    numeric()
  }
  target_mean <- if (ncol(target_right)) {
    as.numeric(crossprod(target_right, target))
  } else {
    numeric()
  }
  list(
    source_mass = source,
    target_mass = target,
    H = outer(source_mean, target_mean),
    G_source = if (ncol(source_right)) {
      unname(crossprod(source_right, source_right * source))
    } else matrix(0, 0L, 0L),
    G_target = if (ncol(target_right)) {
      unname(crossprod(target_right, target_right * target))
    } else matrix(0, 0L, 0L),
    transported_mass = 1
  )
}

.fugw_scalar_shift <- function(stats, source, target, rho, epsilon) {
  source_marginal <- stats$source_marginal
  target_marginal <- stats$target_marginal
  rho[[1L]] * sum(source_marginal *
    (log(source_marginal + 1e-300) - log(source + 1e-300))) +
    rho[[2L]] * sum(target_marginal *
      (log(target_marginal + 1e-300) - log(target + 1e-300))) +
    epsilon * (
      stats$entropy_sum -
        sum(source_marginal * log(source + 1e-300)) -
        sum(target_marginal * log(target + 1e-300))
    )
}

.quadratic_rows <- function(left, gram) {
  if (!ncol(left)) return(numeric(nrow(left)))
  rowSums((left %*% gram) * left)
}

.fugw_dynamic_cost <- function(
    moments,
    stats,
    source_left,
    target_left,
    feature_terms,
    source,
    target,
    rho,
    epsilon,
    feature_weight,
    structure_weight,
    provenance) {
  row <- structure_weight *
    .quadratic_rows(source_left, moments$G_source) +
    (feature_weight / 2) * feature_terms$row
  column <- structure_weight *
    .quadratic_rows(target_left, moments$G_target) +
    (feature_weight / 2) * feature_terms$column
  row <- row + .fugw_scalar_shift(stats, source, target, rho, epsilon)
  left <- matrix(numeric(), nrow = length(source), ncol = 0L)
  right <- matrix(numeric(), nrow = length(target), ncol = 0L)
  if (structure_weight != 0 && ncol(source_left) && ncol(target_left)) {
    left <- cbind(left, -2 * structure_weight * source_left %*% moments$H)
    right <- cbind(right, target_left)
  }
  if (feature_weight != 0 && ncol(feature_terms$left)) {
    left <- cbind(left, (feature_weight / 2) * feature_terms$left)
    right <- cbind(right, feature_terms$right)
  }
  factorized_cost(
    left,
    right,
    row = row,
    column = column,
    exact = isTRUE(provenance$exact),
    provenance = provenance
  )
}

.cost_plan_inner_product <- function(cost, plan) {
  terms <- .cost_native_terms(cost)
  source_mass <- transport_plan_mass(plan, "source")
  target_mass <- transport_plan_mass(plan, "target")
  value <- sum(terms$row * source_mass) +
    sum(terms$column * target_mass)
  if (ncol(terms$left)) {
    value <- value + sum(
      terms$left * transport_plan_apply(plan, terms$right)
    )
  }
  as.numeric(value)
}

.divergence_between_products <- function(
    left, right, reference, left_kl = NULL, right_kl = NULL) {
  left_mass <- sum(left)
  right_mass <- sum(right)
  reference_mass <- sum(reference)
  if (is.null(left_kl)) left_kl <- .generalized_kl_vector(left, reference)
  if (is.null(right_kl)) right_kl <- .generalized_kl_vector(right, reference)
  right_mass * left_kl + left_mass * right_kl +
    (left_mass - reference_mass) * (right_mass - reference_mass)
}

.fugw_factorized_objective <- function(
    sample_plan,
    feature_plan,
    feature_cost,
    feature_moments,
    source_left,
    target_left,
    source,
    target,
    rho,
    epsilon,
    feature_weight,
    structure_weight) {
  sample_state <- .factorized_plan_state(sample_plan)
  feature_state <- .factorized_plan_state(feature_plan)
  structure_cost <- .fugw_dynamic_cost(
    feature_moments,
    feature_state$stats,
    source_left,
    target_left,
    .cost_native_terms(.zero_factorized_cost(length(source), length(target))),
    source,
    target,
    c(0, 0),
    0,
    0,
    1,
    list(exact = TRUE, component = "structure_objective")
  )
  structure_unweighted <- .cost_plan_inner_product(
    structure_cost, sample_plan
  )
  feature_unweighted <- 0.5 * (
    .cost_plan_inner_product(feature_cost, sample_plan) +
      .cost_plan_inner_product(feature_cost, feature_plan)
  )
  source_divergence <- .divergence_between_products(
    sample_state$stats$source_marginal,
    feature_state$stats$source_marginal,
    source
  )
  target_divergence <- .divergence_between_products(
    sample_state$stats$target_marginal,
    feature_state$stats$target_marginal,
    target
  )
  plan_divergence <- .divergence_between_products(
    sample_state$stats$mass,
    feature_state$stats$mass,
    1,
    sample_state$stats$plan_product_kl,
    feature_state$stats$plan_product_kl
  )
  regularization <- rho[[1L]] * source_divergence +
    rho[[2L]] * target_divergence + epsilon * plan_divergence
  list(
    structure_unweighted = structure_unweighted,
    feature_unweighted = feature_unweighted,
    regularization = regularization,
    source_marginal_divergence = source_divergence,
    target_marginal_divergence = target_divergence,
    plan_divergence = plan_divergence,
    fugw_cost = structure_weight * structure_unweighted +
      feature_weight * feature_unweighted + regularization
  )
}

.relative_matrix_change <- function(current, previous) {
  if (!length(current) && !length(previous)) return(0)
  max(abs(current - previous)) / max(1, max(abs(current)), max(abs(previous)))
}

.factorized_geometry_certificate <- function(
    source_cost,
    target_cost,
    source_weights = NULL,
    target_weights = NULL,
    block_size = 32L) {
  source_audit <- geometry_audit(
    source_cost, weights = source_weights, block_size = block_size
  )
  target_audit <- geometry_audit(
    target_cost, weights = target_weights, block_size = block_size
  )
  errors <- c(
    source_audit$relative_error,
    target_audit$relative_error
  )
  exact <- isTRUE(source_audit$exact) && isTRUE(target_audit$exact)
  certified <- isTRUE(source_audit$certified) &&
    isTRUE(target_audit$certified)
  heldout_errors <- c(
    source_audit$heldout_error,
    target_audit$heldout_error
  )
  list(
    exact = exact,
    certified = certified,
    status = if (!certified) {
      "geometry_uncertified"
    } else if (exact) {
      "exact_verified"
    } else {
      "approximate_verified"
    },
    relative_error = if (any(is.finite(errors))) max(errors, na.rm = TRUE) else {
      if (exact) 0 else NA_real_
    },
    heldout_error = if (any(is.finite(heldout_errors))) {
      max(heldout_errors[is.finite(heldout_errors)])
    } else NA_real_,
    source = source_audit,
    target = target_audit,
    interpretation = paste0(
      "Source and target geometry audits are independent of solver, support, ",
      "and hierarchy-transfer certificates."
    )
  )
}

.validate_fugw_init_entry <- function(
    entry, name, n_source, n_target, source_rank, target_rank) {
  required <- c("stats", "moments", "potentials")
  if (!is.list(entry) || !all(required %in% names(entry))) {
    stop(
      sprintf("`init$%s` must contain stats, moments, and potentials.", name),
      call. = FALSE
    )
  }
  stats <- entry$stats
  moments <- entry$moments
  potentials <- entry$potentials
  required_stats <- c(
    "source_marginal", "target_marginal", "mass", "entropy_sum",
    "plan_product_kl"
  )
  if (!is.list(stats) || !all(required_stats %in% names(stats)) ||
      !is.numeric(stats$source_marginal) ||
      length(stats$source_marginal) != n_source ||
      !is.numeric(stats$target_marginal) ||
      length(stats$target_marginal) != n_target ||
      any(!is.finite(stats$source_marginal)) ||
      any(!is.finite(stats$target_marginal)) ||
      any(stats$source_marginal < 0) || any(stats$target_marginal < 0) ||
      !is.numeric(stats$mass) || length(stats$mass) != 1L ||
      !is.finite(stats$mass) || stats$mass <= 0 ||
      !is.numeric(stats$entropy_sum) || length(stats$entropy_sum) != 1L ||
      !is.finite(stats$entropy_sum) ||
      !is.numeric(stats$plan_product_kl) ||
      length(stats$plan_product_kl) != 1L ||
      !is.finite(stats$plan_product_kl)) {
    stop(sprintf("`init$%s$stats` is invalid.", name), call. = FALSE)
  }
  mass_tolerance <- 1e-9 * max(1, stats$mass)
  if (abs(sum(stats$source_marginal) - stats$mass) > mass_tolerance ||
      abs(sum(stats$target_marginal) - stats$mass) > mass_tolerance) {
    stop(sprintf("`init$%s$stats` has inconsistent transported mass.", name),
         call. = FALSE)
  }
  if (!is.list(moments) ||
      !is.matrix(moments$H) ||
      !identical(dim(moments$H), c(source_rank, target_rank)) ||
      !is.matrix(moments$G_source) ||
      !identical(dim(moments$G_source), c(source_rank, source_rank)) ||
      !is.matrix(moments$G_target) ||
      !identical(dim(moments$G_target), c(target_rank, target_rank)) ||
      any(!is.finite(moments$H)) || any(!is.finite(moments$G_source)) ||
      any(!is.finite(moments$G_target))) {
    stop(sprintf("`init$%s$moments` is incompatible with structure factors.",
                 name), call. = FALSE)
  }
  potentials <- .validate_factorized_potential_init(
    potentials, n_source, n_target
  )
  stats$source_marginal <- as.numeric(stats$source_marginal)
  stats$target_marginal <- as.numeric(stats$target_marginal)
  moments$source_mass <- stats$source_marginal
  moments$target_mass <- stats$target_marginal
  moments$transported_mass <- stats$mass
  list(stats = stats, moments = moments, potentials = potentials)
}

.validate_fugw_factorized_init <- function(
    init, n_source, n_target, source_rank, target_rank) {
  if (is.null(init)) return(NULL)
  if (!is.list(init) || !all(c("sample", "feature") %in% names(init))) {
    stop("`init` must contain `sample` and `feature` summaries.",
         call. = FALSE)
  }
  list(
    sample = .validate_fugw_init_entry(
      init$sample, "sample", n_source, n_target, source_rank, target_rank
    ),
    feature = .validate_fugw_init_entry(
      init$feature, "feature", n_source, n_target, source_rank, target_rank
    ),
    provenance = init$provenance %||% list(origin = "user_summary")
  )
}

#' Factorized Matrix-Free Fused Unbalanced Gromov-Wasserstein
#'
#' Solves the same joint-KL two-coupling objective as [fugw_kl()] while
#' representing structure and feature costs by affine-bilinear factors and both
#' couplings by Sinkhorn potentials. The solver never allocates a dense source
#' structure, target structure, feature cost, sample coupling, or feature
#' coupling unless `plan = "dense"` is explicitly requested.
#'
#' For structure costs `Cx = Ax %*% t(Bx)` and
#' `Cy = Ay %*% t(By)`, each block cost is closed by the marginal moments
#' `Bx' diag(qx) Bx`, `By' diag(qy) By`, and `Bx' Q By`. Inner blocks are
#' translation-invariant KL-UOT solves over complete implicit support. The
#' returned certificate separates inner numerical convergence, outer
#' stationarity, geometry representation, and support completeness. Global
#' optimality is never claimed.
#'
#' @section Computational contract:
#' The default operator path has no `n_source * n_target` stored state and
#' prohibits dense `Cx`, `Cy`, `M`, `P`, and `Q` allocations. It is not a
#' subquadratic algorithm: each full-support sweep visits every source-target
#' pair. For effective affine-bilinear rank `r`, bilinear cost evaluation is
#' `O(n_source * n_target * r)` arithmetic, with factor storage linear in the
#' domain sizes and small `H`/Gram workspaces quadratic only in factor ranks.
#' `block_size` controls the reusable pair tile. `plan = "dense"` explicitly
#' leaves this memory contract at the returned-plan boundary.
#'
#' @section Experimental status:
#' This is an experimental candidate. Exact-factor numerical certificates are
#' per-fit claims. Geometry audits describe the represented objective, while
#' planted image/cortical validation and performance measurements are external
#' benchmark evidence; neither is a global-optimality certificate.
#'
#' @param Cx,Cy Square affine-bilinear [factorized_cost()] objects for
#'   source and target structure.
#' @param wx,wy Source and target probability weights. Defaults are uniform.
#' @param reg_marginals One or two positive marginal KL penalties.
#' @param epsilon Positive joint KL regularization.
#' @param alpha Legacy feature coefficient. The structure coefficient is one.
#' @param M Affine-bilinear cross-domain feature cost. `NULL` is the exact zero
#'   cost.
#' @param max_iter Maximum outer block-coordinate iterations.
#' @param tol Positive outer tolerance on both implicit coupling updates and
#'   normalized moment changes.
#' @param max_iter_ot Maximum translation-invariant UOT iterations per block.
#' @param tol_ot Positive inner UOT tolerance.
#' @param rescale_plan Rescale alternating couplings to equal mass, matching
#'   [fugw_kl()].
#' @param check_every Evaluate outer certificates every this many iterations.
#' @param plan Return implicit operators (default) or explicitly materialized
#'   dense couplings.
#' @param block_size Positive native tile size. Workspace is bounded by a
#'   constant number of tiles, not by the full source-by-target shape.
#' @param certify Whether to retain the complete layered certificate. Core
#'   numerical gates are always evaluated.
#' @param init Optional certified coupling summaries for potential, moment, and
#'   mass transfer. This is the internal contract used by
#'   [fugw_multiscale()]; arbitrary dense initial plans are deliberately not
#'   accepted by the matrix-free path.
#' @param feature_weight,structure_weight Explicit nonnegative FUGW
#'   coefficients. Supply both instead of `alpha`.
#' @return An `rfugw_result`. `plans$sample` and `plans$feature` retain the two
#'   couplings; `plan` is the conventional primary sample coupling.
#' @examples
#' x <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
#' y <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), ncol = 2, byrow = TRUE)
#' fit <- fugw_factorized(
#'   sqeuclidean_cost(x), sqeuclidean_cost(y), epsilon = 0.1,
#'   max_iter = 5, max_iter_ot = 200
#' )
#' transport_plan_apply(fit$plans$sample, rep(1, nrow(y)))
#' certificate_layers <- list(
#'   inner_uot = fit$certificate$inner_uot$certified,
#'   final_two_sided_stationarity =
#'     fit$certificate$outer_stationarity$certified &&
#'       fit$certificate$outer_stationarity$final_blocks$certified,
#'   geometry = fit$certificate$geometry$certified,
#'   hierarchy = "not_applicable_single_level",
#'   support = fit$certificate$support,
#'   scientific_validation = "external_frozen_benchmark_not_per_fit"
#' )
#' certificate_layers
#' @export
fugw_factorized <- function(
    Cx,
    Cy,
    wx = NULL,
    wy = NULL,
    reg_marginals = c(10, 10),
    epsilon = 1e-2,
    alpha = 0.5,
    M = NULL,
    max_iter = 100L,
    tol = 1e-7,
    max_iter_ot = 2000L,
    tol_ot = 1e-8,
    rescale_plan = TRUE,
    check_every = 1L,
    plan = c("operator", "dense"),
    block_size = 256L,
    certify = TRUE,
    init = NULL,
    feature_weight = NULL,
    structure_weight = NULL) {
  alpha_was_missing <- missing(alpha)
  plan <- match.arg(plan)
  if (!is.logical(certify) || length(certify) != 1L || is.na(certify)) {
    stop("`certify` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!inherits(Cx, "rfugw_cost_operator") ||
      !inherits(Cy, "rfugw_cost_operator")) {
    stop("`Cx` and `Cy` must be rfugw cost operators.", call. = FALSE)
  }
  source_shape <- cost_shape(Cx)
  target_shape <- cost_shape(Cy)
  if (source_shape[[1L]] != source_shape[[2L]] ||
      target_shape[[1L]] != target_shape[[2L]]) {
    stop("`Cx` and `Cy` must be square structure-cost operators.",
         call. = FALSE)
  }
  source_factors <- .cost_factor_pair(Cx, "Cx")
  target_factors <- .cost_factor_pair(Cy, "Cy")
  n_source <- source_shape[[1L]]
  n_target <- target_shape[[1L]]
  if (is.null(M)) {
    M <- .zero_factorized_cost(
      n_source, n_target, list(role = "feature_cost")
    )
  }
  feature_shape <- cost_shape(M)
  if (!identical(
    unname(feature_shape), as.integer(c(n_source, n_target))
  )) {
    stop("`M` must have shape nrow(Cx) by nrow(Cy).", call. = FALSE)
  }
  feature_terms <- .cost_native_terms(M, "M")
  weights <- .resolve_objective_weights(
    alpha = alpha,
    alpha_was_missing = alpha_was_missing,
    feature_weight = feature_weight,
    structure_weight = structure_weight,
    convention = "fugw_coefficients"
  )
  epsilon <- .validate_positive_scalar(epsilon, "epsilon")
  tol <- .validate_positive_scalar(tol, "tol")
  tol_ot <- .validate_positive_scalar(tol_ot, "tol_ot")
  max_iter <- .validate_count(max_iter, "max_iter")
  max_iter_ot <- .validate_count(max_iter_ot, "max_iter_ot")
  check_every <- .validate_count(check_every, "check_every")
  block_size <- .validate_count(block_size, "block_size")
  if (length(reg_marginals) == 1L) reg_marginals <- rep(reg_marginals, 2L)
  if (!is.numeric(reg_marginals) || length(reg_marginals) != 2L ||
      any(!is.finite(reg_marginals)) || any(reg_marginals <= 0)) {
    stop("`reg_marginals` must be one or two finite positive numbers.",
         call. = FALSE)
  }
  if (is.null(wx)) wx <- rep(1 / n_source, n_source)
  if (is.null(wy)) wy <- rep(1 / n_target, n_target)
  wx <- .assert_prob(wx, n_source, "wx")
  wy <- .assert_prob(wy, n_target, "wy")
  init <- .validate_fugw_factorized_init(
    init,
    n_source,
    n_target,
    ncol(source_factors$right),
    ncol(target_factors$right)
  )

  if (is.null(init)) {
    sample_plan <- .factorized_product_plan(wx, wy, block_size)
    feature_plan <- .factorized_product_plan(wx, wy, block_size)
    sample_moments <- .factorized_product_moments(
      wx, wy, source_factors$right, target_factors$right
    )
    feature_moments <- sample_moments
    sample_summary <- list(
      stats = .factorized_plan_state(sample_plan)$stats,
      potentials = list(
        source_bar = .factorized_plan_state(sample_plan)$source_bar,
        target_bar = .factorized_plan_state(sample_plan)$target_bar
      )
    )
    feature_summary <- list(
      stats = .factorized_plan_state(feature_plan)$stats,
      potentials = list(
        source_bar = .factorized_plan_state(feature_plan)$source_bar,
        target_bar = .factorized_plan_state(feature_plan)$target_bar
      )
    )
  } else {
    sample_plan <- NULL
    feature_plan <- NULL
    sample_moments <- init$sample$moments
    feature_moments <- init$feature$moments
    sample_summary <- init$sample
    feature_summary <- init$feature
  }
  objective_trace <- numeric()
  sample_update_trace <- numeric()
  feature_update_trace <- numeric()
  moment_trace <- numeric()
  plan_screening_trace <- numeric()
  plan_update_audited_trace <- logical()
  objective_check_certified <- logical()
  inner_trace <- list()
  historical_all_inner_certified <- TRUE
  current_inner_certified <- FALSE
  total_inner_iterations <- 0L
  historical_max_inner_residual <- 0
  outer_residual <- Inf
  moment_residual <- Inf
  converged <- FALSE
  final_block_audits <- NULL
  last_sample_fit <- NULL
  last_feature_fit <- NULL
  sample_scale <- 1
  feature_scale <- 1
  iteration_done <- 0L

  for (iteration in seq_len(max_iter)) {
    old_sample <- sample_plan
    old_feature <- feature_plan
    old_sample_moments <- sample_moments
    old_feature_moments <- feature_moments

    sample_state <- if (is.null(sample_plan)) {
      list(stats = sample_summary$stats)
    } else {
      .factorized_plan_state(sample_plan)
    }
    feature_state <- if (is.null(feature_plan)) {
      list(stats = feature_summary$stats)
    } else {
      .factorized_plan_state(feature_plan)
    }
    feature_block_cost <- .fugw_dynamic_cost(
      sample_moments,
      sample_state$stats,
      source_factors$left,
      target_factors$left,
      feature_terms,
      wx,
      wy,
      reg_marginals,
      epsilon,
      weights$feature_weight,
      weights$structure_weight,
      list(
        exact = isTRUE(cost_provenance(Cx)$exact) &&
          isTRUE(cost_provenance(Cy)$exact) && isTRUE(cost_provenance(M)$exact),
        role = "feature_coupling_block",
        outer_iteration = iteration
      )
    )
    feature_mass_reference <- sample_state$stats$mass
    last_feature_fit <- .ot_sinkhorn_unbalanced_ti_factorized(
      feature_block_cost,
      wx,
      wy,
      epsilon * feature_mass_reference,
      reg_marginals * feature_mass_reference,
      max_iter_ot,
      tol_ot,
      "operator",
      init_potentials = feature_summary$potentials,
      block_size = block_size,
      moment_source = source_factors$right,
      moment_target = target_factors$right
    )
    feature_plan <- last_feature_fit$implicit_plan
    feature_scale <- 1
    if (isTRUE(rescale_plan)) {
      feature_mass <- transport_plan_mass(feature_plan)
      feature_scale <- sqrt(feature_mass_reference / feature_mass)
      feature_plan <- .factorized_rescale_plan(feature_plan, feature_scale)
    }
    feature_moments <- .factorized_rescale_moments(
      last_feature_fit$moments, feature_scale
    )
    feature_summary <- list(
      stats = .factorized_plan_state(feature_plan)$stats,
      moments = feature_moments,
      potentials = list(
        source_bar = .factorized_plan_state(feature_plan)$source_bar,
        target_bar = .factorized_plan_state(feature_plan)$target_bar
      )
    )

    feature_state <- .factorized_plan_state(feature_plan)
    sample_block_cost <- .fugw_dynamic_cost(
      feature_moments,
      feature_state$stats,
      source_factors$left,
      target_factors$left,
      feature_terms,
      wx,
      wy,
      reg_marginals,
      epsilon,
      weights$feature_weight,
      weights$structure_weight,
      list(
        exact = isTRUE(cost_provenance(Cx)$exact) &&
          isTRUE(cost_provenance(Cy)$exact) && isTRUE(cost_provenance(M)$exact),
        role = "sample_coupling_block",
        outer_iteration = iteration
      )
    )
    sample_mass_reference <- feature_state$stats$mass
    last_sample_fit <- .ot_sinkhorn_unbalanced_ti_factorized(
      sample_block_cost,
      wx,
      wy,
      epsilon * sample_mass_reference,
      reg_marginals * sample_mass_reference,
      max_iter_ot,
      tol_ot,
      "operator",
      init_potentials = sample_summary$potentials,
      block_size = block_size,
      moment_source = source_factors$right,
      moment_target = target_factors$right
    )
    sample_plan <- last_sample_fit$implicit_plan
    sample_scale <- 1
    if (isTRUE(rescale_plan)) {
      sample_mass <- transport_plan_mass(sample_plan)
      sample_scale <- sqrt(sample_mass_reference / sample_mass)
      sample_plan <- .factorized_rescale_plan(sample_plan, sample_scale)
    }
    sample_moments <- .factorized_rescale_moments(
      last_sample_fit$moments, sample_scale
    )
    sample_summary <- list(
      stats = .factorized_plan_state(sample_plan)$stats,
      moments = sample_moments,
      potentials = list(
        source_bar = .factorized_plan_state(sample_plan)$source_bar,
        target_bar = .factorized_plan_state(sample_plan)$target_bar
      )
    )

    current_inner_residual <- max(
      last_feature_fit$inner_residual,
      last_sample_fit$inner_residual
    )
    historical_max_inner_residual <- max(
      historical_max_inner_residual, current_inner_residual
    )
    current_inner_certified <- isTRUE(last_feature_fit$inner_converged) &&
      isTRUE(last_sample_fit$inner_converged)
    historical_all_inner_certified <- historical_all_inner_certified &&
      current_inner_certified
    rescale_stationarity_tolerance <- max(1e-7, 100 * tol_ot)
    post_rescale_candidate <- current_inner_certified &&
      all(is.finite(c(feature_scale, sample_scale))) &&
      max(abs(log(c(feature_scale, sample_scale)))) <=
        rescale_stationarity_tolerance
    total_inner_iterations <- total_inner_iterations +
      last_feature_fit$iterations + last_sample_fit$iterations
    inner_trace[[length(inner_trace) + 1L]] <- data.frame(
      iteration = rep(iteration, 2L),
      block = c("feature", "sample"),
      iterations = c(
        last_feature_fit$iterations, last_sample_fit$iterations
      ),
      residual = c(
        last_feature_fit$inner_residual, last_sample_fit$inner_residual
      ),
      certified = c(
        last_feature_fit$inner_converged, last_sample_fit$inner_converged
      ),
      post_solve_scale = c(feature_scale, sample_scale),
      post_rescale_candidate = rep(post_rescale_candidate, 2L),
      stringsAsFactors = FALSE
    )

    iteration_done <- iteration
    if (iteration %% check_every == 0L || iteration == max_iter) {
      moment_residual <- max(
        .relative_matrix_change(sample_moments$H, old_sample_moments$H),
        .relative_matrix_change(feature_moments$H, old_feature_moments$H),
        .relative_matrix_change(
          sample_moments$G_source, old_sample_moments$G_source
        ),
        .relative_matrix_change(
          sample_moments$G_target, old_sample_moments$G_target
        ),
        .relative_matrix_change(
          feature_moments$G_source, old_feature_moments$G_source
        ),
        .relative_matrix_change(
          feature_moments$G_target, old_feature_moments$G_target
        ),
        .factorized_relative_change(
          sample_moments$source_mass, old_sample_moments$source_mass
        ),
        .factorized_relative_change(
          sample_moments$target_mass, old_sample_moments$target_mass
        ),
        .factorized_relative_change(
          feature_moments$source_mass, old_feature_moments$source_mass
        ),
        .factorized_relative_change(
          feature_moments$target_mass, old_feature_moments$target_mass
        ),
        abs(sample_moments$transported_mass -
              old_sample_moments$transported_mass),
        abs(feature_moments$transported_mass -
              old_feature_moments$transported_mass)
      )
      sample_screening <- .factorized_plan_screening_residual(
        sample_plan, old_sample
      )
      feature_screening <- .factorized_plan_screening_residual(
        feature_plan, old_feature
      )
      plan_screening_residual <- max(
        sample_screening, feature_screening
      )
      screening_tolerance <- max(1e-5, 100 * tol)
      plan_update_audited <- is.finite(plan_screening_residual) &&
        plan_screening_residual <= screening_tolerance &&
        moment_residual <= screening_tolerance
      if (plan_update_audited) {
        sample_difference <- .factorized_plan_difference(
          sample_plan, old_sample, block_size
        )
        feature_difference <- .factorized_plan_difference(
          feature_plan, old_feature, block_size
        )
        sample_update <- sample_difference$l1
        feature_update <- feature_difference$l1
      } else {
        sample_update <- sample_screening
        feature_update <- feature_screening
      }
      outer_residual <- max(sample_update, feature_update)
      objective <- .fugw_factorized_objective(
        sample_plan,
        feature_plan,
        M,
        feature_moments,
        source_factors$left,
        target_factors$left,
        wx,
        wy,
        reg_marginals,
        epsilon,
        weights$feature_weight,
        weights$structure_weight
      )
      objective_trace <- c(objective_trace, objective$fugw_cost)
      sample_update_trace <- c(sample_update_trace, sample_update)
      feature_update_trace <- c(feature_update_trace, feature_update)
      moment_trace <- c(moment_trace, moment_residual)
      plan_screening_trace <- c(
        plan_screening_trace, plan_screening_residual
      )
      plan_update_audited_trace <- c(
        plan_update_audited_trace, plan_update_audited
      )
      objective_check_certified <- c(
        objective_check_certified, post_rescale_candidate
      )
      iterate_stationary <- plan_update_audited &&
        outer_residual <= tol && moment_residual <= tol
      if (iterate_stationary && current_inner_certified) {
        final_block_audits <- .fugw_final_block_audits(
          sample_plan,
          feature_plan,
          sample_moments,
          feature_moments,
          source_factors$left,
          target_factors$left,
          feature_terms,
          wx,
          wy,
          reg_marginals,
          epsilon,
          weights$feature_weight,
          weights$structure_weight,
          tol_ot,
          block_size,
          sample_scale,
          feature_scale,
          exact = isTRUE(cost_provenance(Cx)$exact) &&
            isTRUE(cost_provenance(Cy)$exact) &&
            isTRUE(cost_provenance(M)$exact)
        )
        current_equal_mass_residual <- abs(
          transport_plan_mass(sample_plan) -
            transport_plan_mass(feature_plan)
        )
        current_equal_mass_tolerance <- max(1e-10, 10 * tol_ot) *
          max(
            1,
            transport_plan_mass(sample_plan),
            transport_plan_mass(feature_plan)
          )
        current_monotonicity_tolerance <- max(1e-10, 100 * tol) *
          max(1, max(abs(objective_trace), 0))
        current_certified_pairs <- if (length(objective_trace) > 1L) {
          head(objective_check_certified, -1L) &
            tail(objective_check_certified, -1L)
        } else logical()
        current_monotonicity_violations <- if (length(objective_trace) > 1L) {
          which(
            diff(objective_trace) > current_monotonicity_tolerance &
              current_certified_pairs
          )
        } else integer()
        if (isTRUE(final_block_audits$certified) &&
            current_equal_mass_residual <= current_equal_mass_tolerance &&
            !length(current_monotonicity_violations) &&
            is.finite(objective$fugw_cost)) {
          converged <- TRUE
          break
        }
      }
    }
  }

  objective <- .fugw_factorized_objective(
    sample_plan,
    feature_plan,
    M,
    feature_moments,
    source_factors$left,
    target_factors$left,
    wx,
    wy,
    reg_marginals,
    epsilon,
    weights$feature_weight,
    weights$structure_weight
  )
  coupling_discrepancy <- .factorized_plan_difference(
    sample_plan, feature_plan, block_size
  )
  final_block_audits <- .fugw_final_block_audits(
    sample_plan,
    feature_plan,
    sample_moments,
    feature_moments,
    source_factors$left,
    target_factors$left,
    feature_terms,
    wx,
    wy,
    reg_marginals,
    epsilon,
    weights$feature_weight,
    weights$structure_weight,
    tol_ot,
    block_size,
    sample_scale,
    feature_scale,
    exact = isTRUE(cost_provenance(Cx)$exact) &&
      isTRUE(cost_provenance(Cy)$exact) &&
      isTRUE(cost_provenance(M)$exact)
  )
  equal_mass_residual <- abs(
    transport_plan_mass(sample_plan) - transport_plan_mass(feature_plan)
  )
  equal_mass_tolerance <- max(1e-10, 10 * tol_ot) * max(
    1,
    transport_plan_mass(sample_plan),
    transport_plan_mass(feature_plan)
  )
  equal_mass_certified <- is.finite(equal_mass_residual) &&
    equal_mass_residual <= equal_mass_tolerance
  geometry <- .factorized_geometry_certificate(
    Cx, Cy, wx, wy, block_size = min(block_size, 64L)
  )
  feature_provenance <- cost_provenance(M)
  objective_tolerance <- max(1e-8, 100 * tol) *
    max(1, abs(objective$fugw_cost))
  monotonicity_tolerance <- max(1e-10, 100 * tol) *
    max(1, max(abs(objective_trace), 0))
  monotonicity_violations <- if (length(objective_trace) > 1L) {
    which(diff(objective_trace) > monotonicity_tolerance)
  } else integer()
  certified_objective_pairs <- if (length(objective_trace) > 1L) {
    head(objective_check_certified, -1L) &
      tail(objective_check_certified, -1L)
  } else logical()
  certification_monotonicity_violations <- if (
      length(objective_trace) > 1L) {
    which(
      diff(objective_trace) > monotonicity_tolerance &
        certified_objective_pairs
    )
  } else integer()
  inner_trace <- do.call(rbind, inner_trace)
  historical_inner_failures <- inner_trace[!inner_trace$certified, , drop = FALSE]
  last_inner_solve_certified <- isTRUE(current_inner_certified)
  final_inner_certified <- last_inner_solve_certified &&
    isTRUE(final_block_audits$certified)
  final_inner_residual <- max(
    last_sample_fit$inner_residual,
    last_feature_fit$inner_residual,
    final_block_audits$residual
  )
  final_inner_tolerance <- max(1e-7, 100 * tol_ot)
  iterate_stationary <- is.finite(outer_residual) &&
    is.finite(moment_residual) &&
    outer_residual <= tol && moment_residual <= tol
  objective_numerically_valid <- all(is.finite(unlist(objective)))
  numerical_finite <- objective_numerically_valid &&
    isTRUE(final_block_audits$sample$numerical_finite) &&
    isTRUE(final_block_audits$feature$numerical_finite) &&
    is.finite(final_inner_residual)
  final_stationarity_certified <- iterate_stationary &&
    final_inner_certified &&
    equal_mass_certified &&
    !length(certification_monotonicity_violations) &&
    numerical_finite
  converged <- final_stationarity_certified
  certification_failures <- names(which(!c(
    consecutive_plan_and_moment_stationarity = iterate_stationary,
    final_inner_solve = last_inner_solve_certified,
    final_sample_block_kkt = isTRUE(final_block_audits$sample$certified),
    final_feature_block_kkt = isTRUE(final_block_audits$feature$certified),
    equal_coupling_mass = equal_mass_certified,
    objective_monotonicity = !length(certification_monotonicity_violations),
    numerical_finiteness = numerical_finite
  )))
  sample_output <- if (identical(plan, "dense")) {
    transport_plan_materialize(sample_plan)
  } else sample_plan
  feature_output <- if (identical(plan, "dense")) {
    transport_plan_materialize(feature_plan)
  } else feature_plan
  out <- list(
    plan = sample_output,
    plans = list(sample = sample_output, feature = feature_output),
    implicit_plans = list(sample = sample_plan, feature = feature_plan),
    formulation = "fugw_factorized_joint_kl",
    backend = "cpp_blocked_affine_bilinear_ti",
    fugw_cost = objective$fugw_cost,
    linear_cost = weights$structure_weight * objective$structure_unweighted,
    feature_cost = weights$feature_weight * objective$feature_unweighted,
    regularization_cost = objective$regularization,
    source_marginal_divergence = objective$source_marginal_divergence,
    target_marginal_divergence = objective$target_marginal_divergence,
    plan_divergence = objective$plan_divergence,
    iterations = iteration_done,
    error = max(outer_residual, moment_residual),
    outer_plan_residual = outer_residual,
    moment_residual = moment_residual,
    sample_update_trace = sample_update_trace,
    feature_update_trace = feature_update_trace,
    moment_trace = moment_trace,
    plan_screening_trace = plan_screening_trace,
    plan_update_audited_trace = plan_update_audited_trace,
    outer_plan_residual_exact = isTRUE(tail(
      plan_update_audited_trace, 1L
    )),
    objective_trace = objective_trace,
    objective_check_certified = objective_check_certified,
    objective_monotone = !length(monotonicity_violations),
    objective_monotonicity_violations = monotonicity_violations,
    certification_monotonicity_violations =
      certification_monotonicity_violations,
    objective_monotonicity_tolerance = monotonicity_tolerance,
    inner_trace = inner_trace,
    historical_inner_failures = historical_inner_failures,
    historical_inner_failure_count = nrow(historical_inner_failures),
    historical_all_inner_certified = historical_all_inner_certified,
    historical_max_inner_residual = historical_max_inner_residual,
    coupling_discrepancy = coupling_discrepancy,
    equal_mass_residual = equal_mass_residual,
    equal_mass_tolerance = equal_mass_tolerance,
    equal_mass_certified = equal_mass_certified,
    moments = list(sample = sample_moments, feature = feature_moments),
    potentials = list(
      sample = list(
        source_bar = .factorized_plan_state(sample_plan)$source_bar,
        target_bar = .factorized_plan_state(sample_plan)$target_bar
      ),
      feature = list(
        source_bar = .factorized_plan_state(feature_plan)$source_bar,
        target_bar = .factorized_plan_state(feature_plan)$target_bar
      )
    ),
    inner_uot_certified = final_inner_certified,
    final_inner_solve_certified = last_inner_solve_certified,
    final_block_optimality_certified = isTRUE(final_block_audits$certified),
    final_stationarity_certified = final_stationarity_certified,
    certification_failures = certification_failures,
    geometry_exact = geometry$exact,
    geometry_certified = geometry$certified,
    geometry_status = geometry$status,
    geometry_relative_error = geometry$relative_error,
    geometry_heldout_error = geometry$heldout_error,
    support_tail_bound = 0,
    global_optimality = "not_claimed",
    feature_weight = weights$feature_weight,
    structure_weight = weights$structure_weight,
    alpha = weights$solver_alpha,
    regularization = epsilon,
    reg_marginals = reg_marginals,
    rescale_plan = isTRUE(rescale_plan),
    block_size = block_size,
    max_tile_elements = min(n_source, block_size) * min(n_target, block_size),
    full_matrix_elements = as.double(n_source) * as.double(n_target),
    plan_representation = transport_plan_representation(sample_output)
  )
  out$initialization <- if (is.null(init)) {
    list(mode = "independent_product", transferred = FALSE)
  } else {
    list(mode = "transferred_summary", transferred = TRUE,
         provenance = init$provenance)
  }
  ans <- .attach_solver_diagnostics(
    out,
    residual = max(outer_residual, moment_residual),
    converged = final_stationarity_certified,
    iterations = iteration_done,
    max_iter = max_iter,
    plan = sample_output,
    inner_residual = final_inner_residual,
    max_inner_residual = final_inner_residual,
    inner_iterations = total_inner_iterations,
    inner_converged = final_inner_certified,
    inner_status = if (final_inner_certified) {
      "converged_post_rescale_two_sided"
    } else {
      "certificate_failure"
    },
    feasibility = "unbalanced",
    feasibility_tol = final_inner_tolerance,
    objective_recomputed = objective$fugw_cost,
    objective_tolerance = objective_tolerance,
    objective_components = list(
      linear_cost = weights$structure_weight * objective$structure_unweighted,
      feature_cost = weights$feature_weight * objective$feature_unweighted,
      regularization_cost = objective$regularization
    ),
    numerical_ok = numerical_finite
  )
  objective_certificate <- isTRUE(ans$objective_consistent) &&
    isTRUE(ans$objective_components_consistent)
  if (!objective_certificate) {
    certification_failures <- unique(c(
      certification_failures, "objective_consistency"
    ))
  }
  fully_certified <- final_stationarity_certified && objective_certificate
  ans$final_stationarity_certified <- fully_certified
  ans$certification_failures <- certification_failures
  ans$termination_reason <- .termination_reason_from_result(ans, max_iter)
  if (isTRUE(ans$converged)) {
    ans$status <- if (!isTRUE(geometry$certified)) {
      "converged_uncertified_geometry"
    } else if (geometry$exact && isTRUE(feature_provenance$exact)) {
      "converged_stationary"
    } else {
      "converged_approximate_stationary"
    }
    ans$warning_payload <- if (!isTRUE(geometry$certified)) {
      list(
        code = "geometry_uncertified",
        message = paste0(
          "The represented FUGW objective is stationary, but the requested ",
          "source or target geometry is not certified."
        )
      )
    } else NULL
  } else if (iterate_stationary && !fully_certified &&
             !(ans$status %in% c("numerical_failure", "objective_mismatch"))) {
    ans$status <- "certificate_failure"
    ans$termination_reason <- "certificate_failure"
    ans$warning_payload <- list(
      code = "certificate_failure",
      message = paste0(
        "Final Moment-FUGW iterate is not certified: ",
        paste(certification_failures, collapse = ", "), "."
      )
    )
  }
  certification_class <- if (!isTRUE(ans$converged)) {
    "uncertified"
  } else if (!isTRUE(geometry$certified)) {
    "stationary_geometry_uncertified"
  } else if (geometry$exact && isTRUE(feature_provenance$exact)) {
    "stationary"
  } else {
    "approximate_stationary"
  }
  ans$certificate <- list(
    inner_uot = list(
      certified = final_inner_certified,
      final_solve_certified = last_inner_solve_certified,
      historical_all_certified = historical_all_inner_certified,
      historical_failure_count = nrow(historical_inner_failures),
      historical_failures = historical_inner_failures,
      last_sample = if (isTRUE(certify)) last_sample_fit$uot_certificate else NULL,
      last_feature = if (isTRUE(certify)) last_feature_fit$uot_certificate else NULL
    ),
    outer_stationarity = list(
      certified = fully_certified,
      classification = certification_class,
      scope = "single_level",
      plan_update_residual = outer_residual,
      plan_update_residual_exact = isTRUE(tail(
        plan_update_audited_trace, 1L
      )),
      plan_screening_residual = tail(plan_screening_trace, 1L),
      plan_screening_tolerance = max(1e-5, 100 * tol),
      moment_residual = moment_residual,
      tolerance = tol,
      objective_monotone = !length(monotonicity_violations),
      objective_monotonicity_violations = monotonicity_violations,
      certification_monotonicity_violations =
        certification_monotonicity_violations,
      equal_mass_residual = equal_mass_residual,
      equal_mass_tolerance = equal_mass_tolerance,
      equal_mass_certified = equal_mass_certified,
      objective_consistent = objective_certificate,
      numerical_finite = numerical_finite,
      failure_reasons = certification_failures,
      final_blocks = final_block_audits,
      interpretation = paste0(
        "Consecutive plan and moment changes are combined with post-rescale ",
        "two-sided block KKT audits. This is not a global-optimality ",
        "certificate."
      )
    ),
    geometry = geometry,
    feature = list(
      exact = isTRUE(feature_provenance$exact),
      provenance = feature_provenance
    ),
    support = list(
      mode = "full_implicit",
      complete = TRUE,
      omitted_kernel_mass_bound = 0,
      adaptive_refinement = FALSE
    ),
    global_optimality = "not_claimed"
  )
  ans$runtime_provenance$memory_contract <- list(
    prohibited_dense_allocations = c("Cx", "Cy", "M", "P", "Q"),
    dense_allocations_requested = identical(plan, "dense"),
    max_tile_elements = out$max_tile_elements,
    full_matrix_elements = out$full_matrix_elements,
    source_structure_rank = ncol(source_factors$left),
    target_structure_rank = ncol(target_factors$left),
    feature_rank = ncol(feature_terms$left),
    factor_tile_copy_policy = "strided_blas_no_factor_tile_copy",
    combined_stats_moments = TRUE,
    moment_workspace_elements = max(
      last_sample_fit$moment_workspace_elements,
      last_feature_fit$moment_workspace_elements
    ),
    moment_workspace_bound = paste0(
      "tile_rows * tile_columns + tile_rows * target_rank + ",
      "source_rank * target_rank + 2 * (n_source + n_target) doubles"
    ),
    moment_workspace_scope = paste0(
      "Peak live native doubles for the reusable plan tile, weighted-target ",
      "tile, log-measure and marginal vectors, and returned cross moment."
    ),
    final_kkt_full_support = TRUE,
    final_kkt_workspace_elements = max(
      final_block_audits$sample$workspace_elements,
      final_block_audits$feature$workspace_elements
    ),
    final_kkt_workspace_bound = paste0(
      "2 * tile_rows * tile_columns + ",
      "2 * (n_source + n_target) doubles"
    )
  )
  ans$runtime_provenance$screening_contract <- list(
    method = "gauge_normalized_potentials_cost_marginals_and_moments",
    tolerance = max(1e-5, 100 * tol),
    checks = length(plan_update_audited_trace),
    exact_plan_update_audits = sum(plan_update_audited_trace),
    screened_plan_updates = sum(!plan_update_audited_trace),
    stationarity_requires_exact_plan_update = TRUE,
    full_final_two_sided_kkt_mandatory = TRUE
  )
  .attach_objective_weight_contract(
    ans,
    weights,
    feature_unweighted = objective$feature_unweighted,
    structure_unweighted = objective$structure_unweighted,
    regularization = objective$regularization
  )
}

#' Factorized FUGW Objective Value
#'
#' @param ... Arguments forwarded to [fugw_factorized()].
#' @return Numeric scalar matrix-free FUGW objective.
#' @export
fugw_factorized2 <- function(...) {
  fugw_factorized(...)$fugw_cost
}
