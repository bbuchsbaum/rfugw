.symmetry_atol <- 1e-10
.symmetry_rtol <- 1e-12

.resolve_symmetric <- function(
    symmetric, C1, C2,
    atol = .symmetry_atol,
    rtol = .symmetry_rtol) {
  detected <- .is_symmetric_cost(C1, atol, rtol) &&
    .is_symmetric_cost(C2, atol, rtol)
  if (is.null(symmetric)) {
    return(detected)
  }
  if (length(symmetric) != 1L || is.na(symmetric)) {
    stop("`symmetric` must be TRUE, FALSE, or NULL.", call. = FALSE)
  }
  symmetric <- isTRUE(symmetric)
  if (symmetric && !detected) {
    stop(
      paste0(
        "`symmetric = TRUE` requires C1 and C2 to be symmetric within ",
        "atol + rtol * scale (atol=1e-10, rtol=1e-12)."
      ),
      call. = FALSE
    )
  }
  symmetric
}

.validate_alpha <- function(alpha) {
  if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be a finite number in [0, 1].", call. = FALSE)
  }
  as.numeric(alpha)
}

.resolve_objective_weights <- function(
    alpha,
    alpha_was_missing,
    feature_weight = NULL,
    structure_weight = NULL,
    convention = c("fgw_share", "fugw_coefficients")) {
  convention <- match.arg(convention)
  aliases_supplied <- !is.null(feature_weight) || !is.null(structure_weight)
  if (xor(is.null(feature_weight), is.null(structure_weight))) {
    stop(
      "`feature_weight` and `structure_weight` must be supplied together.",
      call. = FALSE
    )
  }
  if (aliases_supplied && !isTRUE(alpha_was_missing)) {
    stop(
      "Supply either `alpha` or `feature_weight`/`structure_weight`, not both.",
      call. = FALSE
    )
  }

  if (!aliases_supplied) {
    alpha <- .validate_alpha(alpha)
    if (identical(convention, "fgw_share")) {
      return(list(
        solver_alpha = alpha,
        feature_weight = 1 - alpha,
        structure_weight = alpha,
        requested_feature_weight = 1 - alpha,
        requested_structure_weight = alpha,
        alpha = alpha,
        alpha_convention = "structure_share",
        weights_normalized = FALSE
      ))
    }
    return(list(
      solver_alpha = alpha,
      feature_weight = alpha,
      structure_weight = 1,
      requested_feature_weight = alpha,
      requested_structure_weight = 1,
      alpha = alpha,
      alpha_convention = "feature_coefficient_with_unit_structure",
      weights_normalized = FALSE
    ))
  }

  feature_weight <- .validate_nonneg_scalar(feature_weight, "feature_weight")
  structure_weight <- .validate_nonneg_scalar(structure_weight, "structure_weight")
  total <- feature_weight + structure_weight
  if (!is.finite(total) || total <= 0) {
    stop(
      "`feature_weight` and `structure_weight` cannot both be zero.",
      call. = FALSE
    )
  }
  if (identical(convention, "fgw_share")) {
    effective_feature <- feature_weight / total
    effective_structure <- structure_weight / total
    return(list(
      solver_alpha = effective_structure,
      feature_weight = effective_feature,
      structure_weight = effective_structure,
      requested_feature_weight = feature_weight,
      requested_structure_weight = structure_weight,
      alpha = effective_structure,
      alpha_convention = "structure_share",
      weights_normalized = !isTRUE(all.equal(total, 1))
    ))
  }
  list(
    solver_alpha = feature_weight,
    feature_weight = feature_weight,
    structure_weight = structure_weight,
    requested_feature_weight = feature_weight,
    requested_structure_weight = structure_weight,
    alpha = feature_weight,
    alpha_convention = "feature_coefficient_with_explicit_structure_coefficient",
    weights_normalized = FALSE
  )
}

.attach_objective_weight_contract <- function(
    out,
    weights,
    feature_unweighted,
    structure_unweighted,
    regularization = 0) {
  weighted_feature <- weights$feature_weight * feature_unweighted
  weighted_structure <- weights$structure_weight * structure_unweighted
  out$alpha <- weights$alpha
  out$alpha_convention <- weights$alpha_convention
  out$feature_weight <- weights$feature_weight
  out$structure_weight <- weights$structure_weight
  out$requested_feature_weight <- weights$requested_feature_weight
  out$requested_structure_weight <- weights$requested_structure_weight
  out$weights_normalized <- weights$weights_normalized
  out$objective_decomposition <- list(
    feature_unweighted = as.numeric(feature_unweighted),
    structure_unweighted = as.numeric(structure_unweighted),
    feature_weighted = as.numeric(weighted_feature),
    structure_weighted = as.numeric(weighted_structure),
    regularization = as.numeric(regularization),
    total = as.numeric(weighted_feature + weighted_structure + regularization)
  )
  out
}

.runtime_environment_specs <- function() {
  data.frame(
    name = c(
      "RFUGW_ENABLE_WARM_START", "RFUGW_ADAPTIVE_INNER_TOL",
      "RFUGW_ADAPTIVE_INNER_TOL_STAGE1", "RFUGW_ADAPTIVE_INNER_TOL_STAGE2",
      "RFUGW_SINKHORN_TARGET_TOL_ABS_MULT", "RFUGW_SINKHORN_TARGET_TOL_MULT",
      "RFUGW_SINKHORN_TARGET_TOL_CAP_D", "RFUGW_SINKHORN_TARGET_TOL_CAP_F",
      "RFUGW_SINKHORN_REL_TOL_MULT", "RFUGW_SINKHORN_REL_TOL_FLOOR_D",
      "RFUGW_SINKHORN_REL_TOL_FLOOR_F", "RFUGW_SINKHORN_COL_REL_TOL_MULT",
      "RFUGW_SINKHORN_COL_REL_TOL_FLOOR_D",
      "RFUGW_SINKHORN_COL_REL_TOL_FLOOR_F",
      "RFUGW_MATVEC_BLOCKED_MIN_WORK_D", "RFUGW_MATVEC_GEMV_MIN_WORK_D",
      "RFUGW_MATVEC_BLOCKED_MIN_WORK_F", "RFUGW_MATVEC_GEMV_MIN_WORK_F",
      "RFUGW_AUTOTUNE_MATVEC", "RFUGW_SAMPLED_MIXED",
      "RFUGW_SEMIRELAXED_MIXED", "RFUGW_PIN_BLAS_THREADS",
      "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
      "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"
    ),
    type = c(
      rep("bool", 2), rep("positive_double", 12), rep("nonnegative_integer", 4),
      rep("bool", 4), rep("positive_integer", 5)
    ),
    default = c(
      "0", "0", "1e-6", "5e-7", "100", "52", "1.5e-4", "2e-4",
      "200", "4e-5", "2e-4", "1000", "1.5e-4", "2e-4",
      "25000", "120000", "160000", "25000", "1", "0", "1", "1",
      NA, NA, NA, NA, NA
    ),
    impact = c(
      "FUGW inner warm starts", "FUGW adaptive inner tolerance",
      "FUGW early-stage inner tolerance", "FUGW middle-stage inner tolerance",
      "unbalanced Sinkhorn absolute target", "unbalanced Sinkhorn target multiplier",
      "double target tolerance cap", "float target tolerance cap",
      "unbalanced relative update multiplier", "double relative update floor",
      "float relative update floor", "column-relative multiplier",
      "double column-relative floor", "float column-relative floor",
      "double blocked matvec dispatch", "double BLAS matvec dispatch",
      "float blocked matvec dispatch", "float BLAS matvec dispatch",
      "matvec threshold autotuning", "sampled solver precision",
      "semirelaxed solver precision", "batch BLAS thread pinning",
      "OpenMP thread request", "OpenBLAS thread request", "MKL thread request",
      "Accelerate thread request", "BLIS thread request"
    ),
    stringsAsFactors = FALSE
  )
}

.parse_runtime_environment_value <- function(raw, type, fallback) {
  if (is.na(raw) || !nzchar(raw)) {
    return(list(value = fallback, valid = TRUE, source = "default"))
  }
  valid <- switch(
    type,
    bool = grepl("^[01]$", raw),
    positive_double = grepl(
      "^[+]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?$",
      raw
    ) && is.finite(suppressWarnings(as.numeric(raw))) && as.numeric(raw) > 0,
    nonnegative_integer = grepl("^[0-9]+$", raw) &&
      is.finite(suppressWarnings(as.numeric(raw))),
    positive_integer = grepl("^[0-9]+$", raw) &&
      is.finite(suppressWarnings(as.numeric(raw))) && as.numeric(raw) > 0,
    FALSE
  )
  if (!isTRUE(valid)) {
    return(list(value = fallback, valid = FALSE, source = "invalid_environment_fallback"))
  }
  value <- if (identical(type, "bool")) {
    identical(raw, "1")
  } else {
    as.numeric(raw)
  }
  list(value = value, valid = TRUE, source = "environment")
}

.runtime_env_bool <- function(name, default) {
  raw <- Sys.getenv(name, unset = NA_character_)
  fallback <- if (isTRUE(default)) "1" else "0"
  parsed <- .parse_runtime_environment_value(raw, "bool", fallback)
  if (is.logical(parsed$value)) parsed$value else identical(parsed$value, "1")
}

.sinkhorn_dynamic_range <- function(cost, regularization) {
  regularization <- .validate_positive_scalar(regularization, "regularization")
  cost <- as.numeric(cost)
  if (!length(cost) || any(!is.finite(cost))) {
    stop("Sinkhorn dispatch cost proxy must be finite and nonempty.", call. = FALSE)
  }
  scaled <- cost / regularization
  max(max(abs(scaled)), diff(range(scaled)))
}

.select_sinkhorn_method <- function(
    requested,
    cost,
    regularization,
    precision = "double",
    context = "Sinkhorn") {
  requested <- match.arg(requested, c("scaling", "log", "auto"))
  threshold <- if (identical(precision, "mixed")) 50 else 500
  metric <- .sinkhorn_dynamic_range(cost, regularization)
  unsafe <- metric > threshold
  if (identical(requested, "scaling") && unsafe) {
    stop(
      sprintf(
        "%s scaling-domain request is outside its certified regime: dynamic-range metric %.6g exceeds %.6g. Use `method = \"log\"` or rescale the problem.",
        context, metric, threshold
      ),
      call. = FALSE
    )
  }
  effective <- if (identical(requested, "auto")) {
    if (unsafe) "log" else "scaling"
  } else {
    requested
  }
  list(
    requested = requested,
    effective = effective,
    metric = metric,
    threshold = threshold,
    reason = if (!identical(requested, "auto")) {
      "explicit_request"
    } else if (unsafe) {
      "dynamic_range_exceeds_scaling_threshold"
    } else {
      "dynamic_range_within_scaling_threshold"
    },
    transition = if (identical(requested, effective)) {
      "none"
    } else {
      paste0("auto_to_", effective)
    }
  )
}

.logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

.sinkhorn_unbalanced_log <- function(
    M,
    a,
    b,
    epsilon,
    rho,
    max_iter,
    tol,
    init_plan = NULL,
    base_measure = NULL) {
  active_i <- which(a > 0)
  active_j <- which(b > 0)
  aa <- a[active_i]
  bb <- b[active_j]
  MM <- M[active_i, active_j, drop = FALSE]
  if (is.null(base_measure)) {
    cc <- aa %o% bb
  } else {
    cc <- base_measure[active_i, active_j, drop = FALSE]
  }
  if (any(!is.finite(cc)) || any(cc < 0) || !any(cc > 0)) {
    stop("`base_measure` must be finite, nonnegative, and positive on support.", call. = FALSE)
  }
  logK <- log(cc)
  logK[cc == 0] <- -Inf
  logK <- logK - MM / epsilon
  tau_x <- rho[[1]] / (rho[[1]] + epsilon)
  tau_y <- rho[[2]] / (rho[[2]] + epsilon)
  logu <- rep(0, length(aa))
  logv <- rep(0, length(bb))
  if (!is.null(init_plan) && length(init_plan)) {
    init_active <- init_plan[active_i, active_j, drop = FALSE]
    target_rows <- rowSums(init_active)
    base_rows <- apply(logK, 1L, .logsumexp)
    usable <- target_rows > 0 & is.finite(base_rows)
    logu[usable] <- log(target_rows[usable]) - base_rows[usable]
  }

  err <- Inf
  it <- 0L
  for (k in seq_len(max_iter)) {
    old_logu <- logu
    old_logv <- logv
    logKv <- apply(sweep(logK, 2L, logv, "+"), 1L, .logsumexp)
    logu <- tau_x * (log(aa) - logKv)
    logKtu <- apply(sweep(logK, 1L, logu, "+"), 2L, .logsumexp)
    logv <- tau_y * (log(bb) - logKtu)
    err <- max(abs(logu - old_logu), abs(logv - old_logv))
    it <- k
    if (is.finite(err) && err <= tol) break
  }
  log_plan <- sweep(sweep(logK, 1L, logu, "+"), 2L, logv, "+")
  active_plan <- exp(log_plan)
  plan <- matrix(0, nrow(M), ncol(M))
  plan[active_i, active_j] <- active_plan
  list(
    plan = plan,
    ot_dist = sum(M * plan),
    iterations = as.integer(it),
    error = as.numeric(err),
    inner_iterations = as.integer(it),
    inner_residual = as.numeric(err),
    max_inner_residual = as.numeric(err),
    inner_converged = is.finite(err) && err <= tol,
    inner_status = if (is.finite(err) && err <= tol) "converged" else "max_iter",
    log_potentials = list(source = logu, target = logv),
    active_source = active_i,
    active_target = active_j
  )
}

.automatic_control_provenance <- function(out) {
  keys <- list(
    precision = c("requested_precision", "effective_precision", "backend_transition"),
    outer_tolerance = c("requested_tol", "effective_tol", NA_character_),
    inner_tolerance = c("requested_inner_tol", "effective_inner_tol", NA_character_),
    sinkhorn_method = c(
      "requested_sinkhorn_method", "effective_sinkhorn_method",
      "sinkhorn_backend_transition"
    ),
    structure_rank = c("structure_rank_requested", "structure_rank", NA_character_),
    sample_budget = c("nb_samples_grad_requested", "nb_samples_grad", NA_character_)
  )
  controls <- Map(function(name, fields) {
    requested <- out[[fields[[1]]]]
    effective <- out[[fields[[2]]]]
    if (is.null(requested) && is.null(effective)) return(NULL)
    transition <- if (!is.na(fields[[3]])) out[[fields[[3]]]] else NULL
    changed <- !is.null(requested) && !is.null(effective) &&
      !isTRUE(all.equal(requested, effective))
    list(
      requested = requested,
      effective = effective,
      source = if (isTRUE(changed)) "automatic_policy" else "user_argument_or_default",
      requested_source = out$control_sources[[name]] %||% "untracked",
      transition = transition %||% if (isTRUE(changed)) "adjusted" else "none"
    )
  }, names(keys), keys)
  names(controls) <- names(keys)
  for (name in names(controls)) {
    if (!is.null(controls[[name]]) &&
        identical(controls[[name]]$source, "user_argument_or_default")) {
      controls[[name]]$source <- controls[[name]]$requested_source
    }
  }
  Filter(Negate(is.null), controls)
}

.solver_runtime_provenance <- function(out = list()) {
  specs <- .runtime_environment_specs()
  raw <- Sys.getenv(specs$name, unset = NA_character_)
  parsed <- lapply(seq_len(nrow(specs)), function(i) {
    .parse_runtime_environment_value(raw[[i]], specs$type[[i]], specs$default[[i]])
  })
  native <- if (exists("cpp_runtime_controls", mode = "function")) {
    cpp_runtime_controls()
  } else {
    list(native_controls = "unavailable")
  }
  controls <- lapply(seq_len(nrow(specs)), function(i) {
    list(
      name = specs$name[[i]],
      raw = if (is.na(raw[[i]])) NULL else raw[[i]],
      parsed = parsed[[i]]$value,
      default = if (is.na(specs$default[[i]])) NULL else specs$default[[i]],
      source = parsed[[i]]$source,
      valid = parsed[[i]]$valid,
      impact = specs$impact[[i]]
    )
  })
  names(controls) <- specs$name
  invalid <- specs$name[!vapply(parsed, `[[`, logical(1), "valid")]
  list(
    schema_version = 1L,
    environment = controls,
    native_effective = native,
    automatic = .automatic_control_provenance(out),
    warnings = if (length(invalid)) {
      sprintf(
        "Invalid %s; native/default fallback is recorded.",
        paste(invalid, collapse = ", ")
      )
    } else {
      character()
    }
  )
}

.validate_positive_scalar <- function(x, name) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x <= 0) {
    stop(sprintf("`%s` must be a finite positive number.", name), call. = FALSE)
  }
  as.numeric(x)
}

.validate_nonneg_scalar <- function(x, name) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x < 0) {
    stop(sprintf("`%s` must be a finite nonnegative number.", name), call. = FALSE)
  }
  as.numeric(x)
}

.validate_count <- function(x, name, min = 1L) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x != trunc(x) || x < min || x > .Machine$integer.max) {
    stop(sprintf("`%s` must be an integer >= %d.", name, min), call. = FALSE)
  }
  as.integer(x)
}

.validate_finite_matrix <- function(x, name, square = FALSE) {
  .assert_matrix(x, name)
  if (nrow(x) == 0L || ncol(x) == 0L) {
    stop(sprintf("`%s` must be nonempty.", name), call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(sprintf("`%s` must be finite.", name), call. = FALSE)
  }
  if (square && nrow(x) != ncol(x)) {
    stop(sprintf("`%s` must be square.", name), call. = FALSE)
  }
  x
}

.validate_balanced_plan <- function(G, p, q, name = "G0", tol = 1e-8) {
  .assert_matrix(G, name)
  if (nrow(G) != length(p) || ncol(G) != length(q)) {
    stop(sprintf("`%s` must have shape length(p) x length(q).", name), call. = FALSE)
  }
  if (any(!is.finite(G)) || any(G < 0)) {
    stop(sprintf("`%s` must be finite and nonnegative.", name), call. = FALSE)
  }
  if (max(abs(rowSums(G) - p)) > tol || max(abs(colSums(G) - q)) > tol) {
    stop(sprintf("`%s` must satisfy the marginal constraints.", name), call. = FALSE)
  }
  G
}

.empty_feature_cost <- function() {
  matrix(numeric(0), nrow = 0, ncol = 0)
}

.validate_optional_init_plan <- function(G, ns, nt, name = "init_plan") {
  if (is.null(G) || (is.matrix(G) && length(G) == 0L)) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }
  .assert_matrix(G, name)
  if (nrow(G) != ns || ncol(G) != nt) {
    stop(sprintf("`%s` must have shape %d x %d.", name, ns, nt), call. = FALSE)
  }
  if (any(!is.finite(G)) || any(G < 0)) {
    stop(sprintf("`%s` must be finite and nonnegative.", name), call. = FALSE)
  }
  G
}

.plan_marginal_residuals <- function(plan, p = NULL, q = NULL) {
  if (is.null(plan) || !is.matrix(plan) || length(plan) == 0L) {
    return(list(row_residual = NA_real_, col_residual = NA_real_, mass = NA_real_))
  }
  list(
    row_residual = if (is.null(p)) NA_real_ else max(abs(rowSums(plan) - p)),
    col_residual = if (is.null(q)) NA_real_ else max(abs(colSums(plan) - q)),
    mass = sum(plan)
  )
}

.precision_tolerance_floor <- 1e-6

.resolve_fgw_precision_policy <- function(
    requested_precision,
    requested_tol,
    requested_inner_tol,
    ns,
    nt,
    sinkhorn_method,
    solver) {
  tight_request <- requested_tol < .precision_tolerance_floor ||
    requested_inner_tol < .precision_tolerance_floor
  if (identical(requested_precision, "mixed") && tight_request) {
    effective_precision <- "strict_double"
    transition <- "mixed_to_strict_double_for_tight_tolerance"
  } else if (identical(requested_precision, "strict_double")) {
    effective_precision <- "strict_double"
    transition <- "none"
  } else if (
    identical(requested_precision, "double") && !tight_request &&
      identical(sinkhorn_method, "scaling") && identical(solver, "PGD") &&
      ns >= 32L && nt >= 32L
  ) {
    effective_precision <- "mixed_accelerated"
    transition <- "double_to_mixed_accelerated"
  } else {
    effective_precision <- requested_precision
    transition <- "none"
  }
  list(
    requested_tol = requested_tol,
    effective_tol = requested_tol,
    requested_inner_tol = requested_inner_tol,
    effective_inner_tol = requested_inner_tol,
    requested_precision = requested_precision,
    effective_precision = effective_precision,
    backend_transition = transition,
    automatic_backend_transition = !identical(transition, "none")
  )
}

.termination_reason_from_result <- function(out, max_iter) {
  if (isTRUE(out$converged)) {
    return("tolerance")
  }
  if (out$status %in% c(
    "inner_failure", "lp_failure", "numerical_failure",
    "infeasible", "objective_mismatch", "stagnation"
  )) {
    return(out$status)
  }
  if (isTRUE(out$iterations >= max_iter)) {
    return("max_iter")
  }
  if (is.finite(out$residual)) {
    return("stagnation")
  }
  out$status %||% "failure"
}

.attach_solver_diagnostics <- function(out,
                                       residual,
                                       converged,
                                       iterations,
                                       max_iter,
                                       p = NULL,
                                       q = NULL,
                                       plan = NULL,
                                       inner_residual = NA_real_,
                                       max_inner_residual = NA_real_,
                                       inner_iterations = NA_integer_,
                                       inner_converged = NA,
                                       inner_status = NULL,
                                       feasibility = c(
                                         "auto", "none", "balanced", "partial",
                                         "semirelaxed", "unbalanced"
                                       ),
                                       feasibility_tol = NA_real_,
                                       mass_target = NULL,
                                       objective_recomputed = NA_real_,
                                       objective_tolerance = 1e-8,
                                       objective_components = list(),
                                       numerical_ok = TRUE,
                                       lp_ok = TRUE) {
  feasibility <- match.arg(feasibility)
  if (is.null(plan)) {
    plan <- out$plan
    if (is.null(plan) && !is.null(out$pi_samp)) {
      plan <- out$pi_samp
    }
  }
  if (identical(feasibility, "auto")) {
    feasibility <- if (!is.null(p) && !is.null(q)) "balanced" else "none"
  }
  if (!is.finite(feasibility_tol)) {
    feasibility_tol <- 1e-8
  }
  residuals <- .plan_marginal_residuals(plan, p, q)
  if (identical(feasibility, "partial") && !is.null(plan) &&
      !is.null(p) && !is.null(q) && !is.null(mass_target)) {
    residuals$row_residual <- max(c(rowSums(plan) - p, 0))
    residuals$col_residual <- max(c(colSums(plan) - q, 0))
  }
  plan_ok <- is.null(plan) || (is.matrix(plan) && all(is.finite(plan)) && all(plan >= 0))
  value_fields <- c("objective", "ot_dist", "fgw_dist", "gw_dist", "fugw_cost", "ucoot_cost",
                    "srgw_dist", "srfgw_dist", "partial_gw_dist", "partial_fgw_dist")
  value_ok <- TRUE
  for (nm in value_fields) {
    if (!is.null(out[[nm]]) && !all(is.finite(out[[nm]]))) {
      value_ok <- FALSE
    }
  }
  has_inner_certificate <- length(inner_converged) == 1L && !is.na(inner_converged)
  finite_inner_residuals <- !has_inner_certificate || (
    length(inner_residual) == 1L && is.finite(inner_residual) &&
      length(max_inner_residual) == 1L && is.finite(max_inner_residual)
  )
  inner_ok <- !has_inner_certificate ||
    (isTRUE(inner_converged) && finite_inner_residuals)

  mass_residual <- if (is.null(mass_target) || is.null(plan)) {
    NA_real_
  } else {
    abs(sum(plan) - mass_target)
  }
  mass_target_reported <- !is.null(mass_target)
  mass_certified <- identical(feasibility, "partial") &&
    mass_target_reported && is.finite(mass_residual) &&
    mass_residual <= feasibility_tol
  feasibility_residual <- switch(
    feasibility,
    none = 0,
    balanced = max(residuals$row_residual, residuals$col_residual),
    partial = max(residuals$row_residual, residuals$col_residual, mass_residual),
    semirelaxed = residuals$row_residual,
    unbalanced = if (has_inner_certificate) max(inner_residual, max_inner_residual) else Inf
  )
  feasibility_ok <- isTRUE(plan_ok) && is.finite(feasibility_residual) &&
    feasibility_residual <= feasibility_tol

  reported_names <- c(
    "ot_dist", "fgw_dist", "gw_dist", "fugw_cost", "ucoot_cost",
    "srgw_dist", "srfgw_dist", "partial_gw_dist", "partial_fgw_dist", "objective"
  )
  reported_objective <- NA_real_
  for (nm in reported_names) {
    if (!is.null(out[[nm]])) {
      reported_objective <- as.numeric(out[[nm]])[1]
      break
    }
  }
  has_objective_check <- length(objective_recomputed) == 1L && !is.na(objective_recomputed)
  objective_residual <- if (has_objective_check && is.finite(reported_objective)) {
    abs(reported_objective - objective_recomputed)
  } else {
    NA_real_
  }
  objective_consistent <- if (has_objective_check) {
    is.finite(objective_recomputed) && is.finite(reported_objective) &&
      is.finite(objective_residual) && objective_residual <= objective_tolerance
  } else {
    isTRUE(value_ok)
  }
  component_consistency <- logical()
  if (length(objective_components)) {
    for (nm in names(objective_components)) {
      recomputed <- as.numeric(objective_components[[nm]])[1]
      reported <- as.numeric(out[[nm]])[1]
      component_residual <- abs(reported - recomputed)
      component_ok <- length(reported) == 1L && is.finite(reported) &&
        is.finite(recomputed) && is.finite(component_residual) &&
        component_residual <= objective_tolerance
      out[[paste0(nm, "_recomputed")]] <- recomputed
      out[[paste0(nm, "_residual")]] <- component_residual
      out[[paste0(nm, "_consistent")]] <- component_ok
      component_consistency[[nm]] <- component_ok
    }
  }
  objective_components_consistent <- !length(component_consistency) || all(component_consistency)
  objective_ok <- isTRUE(value_ok) && isTRUE(objective_consistent) &&
    isTRUE(objective_components_consistent)

  if (!isTRUE(numerical_ok) || !isTRUE(plan_ok) || !isTRUE(value_ok)) {
    status <- "numerical_failure"
  } else if (!isTRUE(lp_ok)) {
    status <- "lp_failure"
  } else if (!inner_ok && isTRUE(converged)) {
    status <- "inner_failure"
  } else if (!feasibility_ok && isTRUE(converged)) {
    status <- "infeasible"
  } else if (!objective_ok && isTRUE(converged)) {
    status <- "objective_mismatch"
  } else if (isTRUE(converged)) {
    status <- "converged"
  } else {
    status <- "max_iter"
  }
  out$residual <- as.numeric(residual)
  if (is.null(out$error)) {
    out$error <- out$residual
  }
  out$status <- status
  out$converged <- identical(status, "converged")
  out$iterations <- as.integer(iterations)
  out$max_iter <- as.integer(max_iter)
  out$row_residual <- residuals$row_residual
  out$col_residual <- residuals$col_residual
  out$mass <- residuals$mass
  out$mass_residual <- mass_residual
  out$mass_target <- if (mass_target_reported) as.numeric(mass_target) else NA_real_
  out$mass_certified <- if (identical(feasibility, "partial")) {
    isTRUE(mass_certified)
  } else {
    NA
  }
  out$mass_certification <- if (!identical(feasibility, "partial")) {
    "not_applicable"
  } else if (!mass_target_reported) {
    "not_certified_missing_target"
  } else if (isTRUE(mass_certified)) {
    "certified"
  } else {
    "not_certified_mass_residual"
  }
  out$feasibility <- feasibility
  out$feasibility_residual <- feasibility_residual
  out$feasibility_tolerance <- feasibility_tol
  out$feasible <- feasibility_ok
  out$inner_residual <- as.numeric(inner_residual)
  out$max_inner_residual <- as.numeric(max_inner_residual)
  out$inner_iterations <- as.integer(inner_iterations)
  out$inner_converged <- if (has_inner_certificate) isTRUE(inner_converged) else NA
  out$inner_status <- inner_status %||% if (!has_inner_certificate) {
    "not_reported"
  } else if (isTRUE(inner_converged)) {
    "converged"
  } else {
    "failed"
  }
  out$objective_recomputed <- as.numeric(objective_recomputed)
  out$objective_residual <- objective_residual
  out$objective_tolerance <- objective_tolerance
  out$objective_consistent <- objective_consistent
  out$objective_components_consistent <- objective_components_consistent
  out$warning_payload <- if (identical(status, "converged")) {
    NULL
  } else {
    list(
      code = status,
      message = sprintf("Solver result is not certified: %s.", status)
    )
  }
  if (is.null(out$formulation)) {
    out$formulation <- if (!is.null(out$fugw_cost)) {
      "fugw_kl"
    } else if (!is.null(out$ucoot_cost)) {
      "ucoot"
    } else if (!is.null(out$gw_dist) && is.null(out$fgw_dist)) {
      "gw"
    } else if (!is.null(out$srgw_dist)) {
      "semirelaxed_gw"
    } else if (!is.null(out$srfgw_dist)) {
      "semirelaxed_fgw"
    } else {
      "fgw"
    }
  }
  if (is.null(out$regularization)) {
    out$regularization <- NA_real_
  }
  if (is.null(out$backend)) {
    out$backend <- "cpp"
  }
  out$runtime_provenance <- .solver_runtime_provenance(out)
  class(out) <- unique(c("rfugw_result", class(out)))
  out
}

.reject_unused_dots <- function(...) {
  extras <- list(...)
  if (length(extras)) {
    nms <- names(extras)
    if (is.null(nms) || any(!nzchar(nms))) {
      stop("Unused unnamed arguments are not accepted.", call. = FALSE)
    }
    stop(
      sprintf(
        "Unsupported argument(s): %s. See inst/solver-contract.md.",
        paste(nms, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.parse_sampled_budget <- function(nb_samples_grad, ns, nt) {
  ns <- .validate_count(ns, "ns")
  nt <- .validate_count(nt, "nt")

  if (length(nb_samples_grad) == 1L) {
    nb <- .validate_count(nb_samples_grad, "nb_samples_grad")
    if (nb > ns) {
      nb_p <- ns
      nb_q <- max(1L, as.integer(nb %/% ns))
    } else {
      nb_p <- nb
      nb_q <- 1L
    }
  } else {
    if (length(nb_samples_grad) != 2L || !is.numeric(nb_samples_grad) ||
        any(!is.finite(nb_samples_grad)) ||
        any(nb_samples_grad != trunc(nb_samples_grad)) ||
        any(nb_samples_grad < 1) ||
        any(nb_samples_grad > .Machine$integer.max)) {
      stop("`nb_samples_grad` must be an integer >=1 or a length-2 integer vector.", call. = FALSE)
    }
    nb <- as.integer(nb_samples_grad)
    nb_p <- nb[[1]]
    nb_q <- nb[[2]]
  }

  clamped <- FALSE
  msgs <- character()
  if (nb_p > ns) {
    msgs <- c(msgs, sprintf("%d source samples but ns = %d", nb_p, ns))
    nb_p <- ns
    clamped <- TRUE
  }
  if (nb_q > nt) {
    msgs <- c(msgs, sprintf("%d target samples but nt = %d", nb_q, nt))
    nb_q <- nt
    clamped <- TRUE
  }
  if (length(msgs)) {
    warning(
      paste0("`nb_samples_grad` requested ", paste(msgs, collapse = "; "), "; clamping."),
      call. = FALSE
    )
  }

  list(nb_p = as.integer(nb_p), nb_q = as.integer(nb_q), clamped = clamped)
}

.parse_lowrank_rank <- function(rank, ns, nt) {
  rmax <- min(.validate_count(ns, "ns"), .validate_count(nt, "nt"))
  if (is.null(rank)) {
    return(as.integer(rmax))
  }
  r <- .validate_count(rank, "rank")
  if (r > rmax) {
    warning(
      sprintf("`rank` = %d exceeds min(ns, nt) = %d; clamping.", r, rmax),
      call. = FALSE
    )
    r <- rmax
  }
  as.integer(r)
}
