#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_arg <- function(prefix, default = NULL) {
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

input_path <- value_arg("--input=")
output_path <- value_arg("--output=", ".gate/pot-oracle/rfugw-comparison.json")
scope <- value_arg("--scope=", "pr")
if (is.null(input_path) || !file.exists(input_path)) {
  stop("Use --input=<live POT JSON receipt>.", call. = FALSE)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("The POT comparator requires package `jsonlite`.", call. = FALSE)
}
if (!requireNamespace("rfugw", quietly = TRUE)) {
  stop("Install rfugw before running the live POT comparator.", call. = FALSE)
}

oracle <- jsonlite::fromJSON(input_path, simplifyVector = FALSE)
if (!identical(oracle$schema_version, "rfugw-pot-oracle-v1")) {
  stop("Unsupported POT oracle schema: ", oracle$schema_version, call. = FALSE)
}

as_vector <- function(x) as.numeric(unlist(x, recursive = TRUE, use.names = FALSE))
as_matrix <- function(x) {
  if (is.matrix(x)) return(unname(x))
  rows <- lapply(x, as_vector)
  if (!length(rows)) return(matrix(numeric(), 0L, 0L))
  unname(do.call(rbind, rows))
}
scalar <- function(x) as.numeric(x)[[1L]]
integer_scalar <- function(x) as.integer(x)[[1L]]
logical_scalar <- function(x) isTRUE(x)

results <- list()
exceptions <- list()

record_metric <- function(
    case, metric, actual, expected, atol, rtol, strict = TRUE, note = NULL) {
  actual <- as.numeric(actual)
  expected <- as.numeric(expected)
  shape_ok <- identical(length(actual), length(expected))
  finite_ok <- length(actual) && length(expected) &&
    all(is.finite(actual)) && all(is.finite(expected))
  if (shape_ok && finite_ok) {
    scale <- pmax(abs(actual), abs(expected), .Machine$double.eps)
    difference <- abs(actual - expected)
    bound <- atol + rtol * scale
    ok <- all(difference <= bound)
    max_abs <- max(difference)
    max_ratio <- max(difference / pmax(bound, .Machine$double.eps))
  } else {
    ok <- FALSE
    max_abs <- Inf
    max_ratio <- Inf
  }
  results[[length(results) + 1L]] <<- list(
    case = case,
    metric = metric,
    strict = isTRUE(strict),
    ok = isTRUE(ok),
    max_abs = max_abs,
    max_tolerance_ratio = max_ratio,
    atol = atol,
    rtol = rtol,
    actual_length = length(actual),
    expected_length = length(expected),
    note = note
  )
  invisible(ok)
}

record_boolean <- function(case, metric, ok, strict = TRUE, note = NULL) {
  results[[length(results) + 1L]] <<- list(
    case = case,
    metric = metric,
    strict = isTRUE(strict),
    ok = isTRUE(ok),
    max_abs = if (isTRUE(ok)) 0 else Inf,
    max_tolerance_ratio = if (isTRUE(ok)) 0 else Inf,
    atol = 0,
    rtol = 0,
    actual_length = 1L,
    expected_length = 1L,
    note = note
  )
  invisible(ok)
}

case_id <- function(case) sprintf("%s-seed-%s", case$family, case$seed)

diagnostic_text <- function(x) {
  warning_count <- length(x$warnings %||% list())
  stdout <- x$stdout %||% ""
  stderr <- x$stderr %||% ""
  c(warnings = warning_count, stdout = nzchar(stdout), stderr = nzchar(stderr))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

check_diagnostics <- function(case) {
  id <- case_id(case)
  ignored <- switch(
    case$family,
    partial_gromov = "entropic_partial_fgw",
    translation_invariant_uot = "translation_invariant",
    sampled_gromov_quality = names(case$diagnostics),
    character()
  )
  for (name in names(case$diagnostics)) {
    if (name %in% ignored) next
    diagnostic <- diagnostic_text(case$diagnostics[[name]])
    record_boolean(
      id,
      paste0("pot_diagnostics_clean_", name),
      all(diagnostic == 0),
      strict = TRUE,
      note = paste(names(diagnostic), diagnostic, collapse = "; ")
    )
  }
}

run_balanced_linear <- function(case) {
  id <- case_id(case)
  x <- case$inputs
  p <- case$params
  y <- case$outputs
  cost <- as_matrix(x$cost)
  source <- as_vector(x$source)
  target <- as_vector(x$target)
  epsilon <- scalar(p$epsilon)
  rho <- as_vector(p$rho)
  max_iter <- integer_scalar(p$max_iter)
  tol <- scalar(p$tol)
  exact <- rfugw::ot_emd(cost, source, target, max_iter = max_iter, tol = tol)
  scaling <- rfugw::ot_sinkhorn(
    cost, source, target, epsilon = epsilon, method = "scaling",
    max_iter = max_iter, tol = tol
  )
  logarithmic <- rfugw::ot_sinkhorn(
    cost, source, target, epsilon = epsilon, method = "log",
    max_iter = max_iter, tol = tol
  )
  unbalanced <- rfugw::ot_sinkhorn_unbalanced(
    cost, source, target, epsilon = epsilon, rho = rho,
    max_iter = max_iter, tol = tol
  )
  a <- scalar(case$tolerance$atol); r <- scalar(case$tolerance$rtol)
  record_metric(id, "exact_plan", exact$plan, as_matrix(y$exact_plan), a, r)
  record_metric(id, "exact_linear_cost", exact$ot_dist, y$exact_linear_cost, a, r)
  record_metric(id, "scaling_plan", scaling$plan, as_matrix(y$scaling_plan), a, r)
  record_metric(id, "scaling_linear_cost", scaling$ot_dist, y$scaling_linear_cost, a, r)
  record_metric(id, "scaling_regularized_objective", scaling$regularized_objective,
                y$scaling_regularized_objective, a, r)
  record_metric(id, "log_plan", logarithmic$plan, as_matrix(y$log_plan), a, r)
  record_metric(id, "log_regularized_objective", logarithmic$regularized_objective,
                y$log_regularized_objective, a, r)
  record_metric(id, "unbalanced_plan", unbalanced$plan,
                as_matrix(y$unbalanced_plan), a, r)
  record_metric(id, "unbalanced_mass", unbalanced$mass,
                y$unbalanced_mass, a, r)
  record_metric(id, "unbalanced_regularized_objective",
                unbalanced$regularized_objective,
                y$unbalanced_regularized_objective, a, r)
}

run_partial_linear <- function(case) {
  id <- case_id(case)
  x <- case$inputs; p <- case$params; y <- case$outputs
  cost <- as_matrix(x$cost); source <- as_vector(x$source); target <- as_vector(x$target)
  mass <- scalar(p$mass); epsilon <- scalar(p$epsilon)
  penalty <- scalar(p$discard_penalty); max_iter <- integer_scalar(p$max_iter)
  tol <- scalar(p$tol)
  exact <- rfugw::ot_partial_emd(
    cost, source, target, mass = mass, max_iter = max_iter, tol = tol
  )
  scaling <- rfugw::ot_partial_sinkhorn(
    cost, source, target, mass = mass, epsilon = epsilon,
    method = "scaling", max_iter = max_iter, tol = tol
  )
  logarithmic <- rfugw::ot_partial_sinkhorn(
    cost, source, target, mass = mass, epsilon = epsilon,
    method = "log", max_iter = max_iter, tol = tol
  )
  penalized <- rfugw::ot_partial_penalized(
    cost, source, target, discard_penalty = penalty,
    max_iter = max_iter, tol = tol
  )
  a <- scalar(case$tolerance$atol); r <- scalar(case$tolerance$rtol)
  record_metric(id, "exact_plan", exact$plan, as_matrix(y$exact_plan), a, r)
  record_metric(id, "exact_linear_cost", exact$ot_dist, y$exact_linear_cost, a, r)
  record_metric(id, "scaling_plan", scaling$plan, as_matrix(y$scaling_plan), a, r)
  record_metric(id, "scaling_objective", scaling$regularized_objective,
                y$scaling_objective, a, r)
  record_metric(id, "log_plan", logarithmic$plan, as_matrix(y$log_plan), a, r)
  record_metric(id, "log_objective", logarithmic$regularized_objective,
                y$log_objective, a, r)
  record_metric(id, "penalized_plan", penalized$plan,
                as_matrix(y$penalized_plan), a, r)
  record_metric(id, "penalized_mass", penalized$transported_mass,
                y$penalized_mass, a, r)
  record_metric(id, "penalized_objective", penalized$penalized_partial_objective,
                y$penalized_objective, a, r)
}

run_balanced_gromov <- function(case) {
  id <- case_id(case)
  x <- case$inputs; p <- case$params; y <- case$outputs
  c1 <- as_matrix(x$source_cost); c2 <- as_matrix(x$target_cost)
  feature <- as_matrix(x$feature_cost); source <- as_vector(x$source)
  target <- as_vector(x$target); init <- as_matrix(x$init_plan)
  symmetric <- logical_scalar(p$symmetric); epsilon <- scalar(p$epsilon)
  alpha <- scalar(p$alpha); max_iter <- integer_scalar(p$max_iter); tol <- scalar(p$tol)
  gw_exact <- rfugw::gromov_wasserstein(
    c1, c2, source, target, symmetric = symmetric, G0 = init,
    max_iter = max_iter, tol_rel = tol, tol_abs = tol,
    lp_solver = "cpp_transport"
  )
  fgw_exact <- rfugw::fgw_exact_cg(
    feature, c1, c2, source, target, alpha = alpha,
    symmetric = symmetric, G0 = init, max_iter = max_iter,
    tol_rel = tol, tol_abs = tol, lp_solver = "cpp_transport"
  )
  gw_entropic <- rfugw::entropic_gromov_wasserstein(
    c1, c2, source, target, epsilon = epsilon, symmetric = symmetric,
    G0 = init, max_iter = max_iter, tol = tol, solver = "PGD",
    sinkhorn_max_iter = 5000L, sinkhorn_tol = 1e-12,
    sinkhorn_method = "log", precision = "double"
  )
  fgw_entropic <- rfugw::fgw_entropic(
    feature, c1, c2, source, target, alpha = alpha,
    epsilon = epsilon, symmetric = symmetric, init_plan = init,
    max_iter = max_iter, tol = tol, solver = "PGD",
    sinkhorn_max_iter = 5000L, sinkhorn_tol = 1e-12,
    sinkhorn_method = "log", precision = "double"
  )
  a <- scalar(case$tolerance$atol); r <- scalar(case$tolerance$rtol)
  record_metric(id, "gw_exact_plan", gw_exact$plan, as_matrix(y$gw_exact_plan), a, r)
  record_metric(id, "gw_exact_objective", gw_exact$gw_dist, y$gw_exact_objective, a, r)
  record_metric(id, "fgw_exact_plan", fgw_exact$plan, as_matrix(y$fgw_exact_plan), a, r)
  record_metric(id, "fgw_exact_objective", fgw_exact$fgw_dist, y$fgw_exact_objective, a, r)
  record_metric(id, "gw_entropic_plan", gw_entropic$plan,
                as_matrix(y$gw_entropic_plan), a, r)
  record_metric(id, "gw_entropic_objective", gw_entropic$gw_dist,
                y$gw_entropic_objective, a, r)
  record_metric(id, "fgw_entropic_plan", fgw_entropic$plan,
                as_matrix(y$fgw_entropic_plan), a, r)
  record_metric(id, "fgw_entropic_objective", fgw_entropic$fgw_dist,
                y$fgw_entropic_objective, a, r)
}

run_partial_gromov <- function(case) {
  id <- case_id(case)
  x <- case$inputs; p <- case$params; y <- case$outputs
  c1 <- as_matrix(x$source_cost); c2 <- as_matrix(x$target_cost)
  feature <- as_matrix(x$feature_cost); source <- as_vector(x$source)
  target <- as_vector(x$target); init <- as_matrix(x$init_plan)
  mass <- scalar(p$mass); alpha <- scalar(p$alpha); epsilon <- scalar(p$epsilon)
  symmetric <- logical_scalar(p$symmetric); max_iter <- integer_scalar(p$max_iter)
  tol <- scalar(p$tol)
  pgw <- rfugw::partial_gromov_wasserstein(
    c1, c2, source, target, m = mass, symmetric = symmetric, G0 = init,
    numItermax = max_iter, tol = tol, log = TRUE,
    lp_solver = "cpp_transport"
  )
  pfgw <- rfugw::partial_fused_gromov_wasserstein(
    feature, c1, c2, source, target, m = mass, alpha = alpha,
    symmetric = symmetric, G0 = init, numItermax = max_iter,
    tol = tol, log = TRUE, lp_solver = "cpp_transport"
  )
  epgw <- rfugw::entropic_partial_gromov_wasserstein(
    c1, c2, source, target, reg = epsilon, m = mass,
    symmetric = symmetric, G0 = init, numItermax = 300L,
    tol = 1e-8, log = TRUE
  )
  epfgw <- rfugw::entropic_partial_fused_gromov_wasserstein(
    feature, c1, c2, source, target, reg = epsilon, m = mass,
    alpha = alpha, symmetric = symmetric, G0 = init,
    numItermax = 300L, tol = 1e-8, log = TRUE
  )
  a <- scalar(case$tolerance$atol); r <- scalar(case$tolerance$rtol)
  record_metric(id, "partial_gw_plan", pgw$plan, as_matrix(y$partial_gw_plan), a, r)
  record_metric(id, "partial_gw_objective", pgw$partial_gw_dist,
                y$partial_gw_objective, a, r)
  record_metric(id, "partial_fgw_plan", pfgw$plan,
                as_matrix(y$partial_fgw_plan), a, r)
  record_metric(id, "partial_fgw_objective", pfgw$partial_fgw_dist,
                y$partial_fgw_objective, a, r)
  record_metric(id, "entropic_partial_gw_plan", epgw$plan,
                as_matrix(y$entropic_partial_gw_plan), a, r)
  record_metric(id, "entropic_partial_gw_objective", epgw$partial_gw_dist,
                y$entropic_partial_gw_objective, a, r)
  record_metric(
    id, "entropic_partial_fgw_plan_monitor", epfgw$plan,
    as_matrix(y$entropic_partial_fgw_plan_monitor), a, r,
    strict = FALSE,
    note = case$exception$reason
  )
  record_metric(
    id, "entropic_partial_fgw_independent_objective",
    epfgw$partial_fgw_dist,
    rfugw::ot_fgw_square(feature, c1, c2, epfgw, alpha = alpha),
    1e-9, 1e-8, strict = TRUE,
    note = "rfugw square-loss objective is authoritative for this known POT gradient exception"
  )
  exceptions[[length(exceptions) + 1L]] <<- c(
    list(case = id), case$exception
  )
}

run_semirelaxed <- function(case) {
  id <- case_id(case)
  x <- case$inputs; p <- case$params; y <- case$outputs
  c1 <- as_matrix(x$source_cost); c2 <- as_matrix(x$target_cost)
  feature <- as_matrix(x$feature_cost); source <- as_vector(x$source)
  init <- as_matrix(x$init_plan); symmetric <- logical_scalar(p$symmetric)
  epsilon <- scalar(p$epsilon); alpha <- scalar(p$alpha)
  max_iter <- integer_scalar(p$max_iter); tol <- scalar(p$tol)
  gw_exact <- rfugw::semirelaxed_gromov_wasserstein(
    c1, c2, source, symmetric = symmetric, G0 = init,
    max_iter = max_iter, tol_rel = tol, tol_abs = tol
  )
  fgw_exact <- rfugw::semirelaxed_fused_gromov_wasserstein(
    feature, c1, c2, source, symmetric = symmetric, alpha = alpha,
    G0 = init, max_iter = max_iter, tol_rel = tol, tol_abs = tol
  )
  gw_entropic <- rfugw::entropic_semirelaxed_gromov_wasserstein(
    c1, c2, source, epsilon = epsilon, symmetric = symmetric,
    G0 = init, max_iter = max_iter, tol = tol,
    precision = "double", backend = "cpp"
  )
  fgw_entropic <- rfugw::entropic_semirelaxed_fused_gromov_wasserstein(
    feature, c1, c2, source, epsilon = epsilon, alpha = alpha,
    symmetric = symmetric, G0 = init, max_iter = max_iter,
    tol = tol, precision = "double", backend = "cpp"
  )
  a <- scalar(case$tolerance$atol); r <- scalar(case$tolerance$rtol)
  values <- list(
    gw_exact = list(gw_exact$plan, gw_exact$srgw_dist),
    fgw_exact = list(fgw_exact$plan, fgw_exact$srfgw_dist),
    gw_entropic = list(gw_entropic$plan, gw_entropic$srgw_dist),
    fgw_entropic = list(fgw_entropic$plan, fgw_entropic$srfgw_dist)
  )
  for (name in names(values)) {
    record_metric(id, paste0(name, "_plan"), values[[name]][[1L]],
                  as_matrix(y[[paste0(name, "_plan")]]), a, r)
    record_metric(id, paste0(name, "_objective"), values[[name]][[2L]],
                  y[[paste0(name, "_objective")]], a, r)
  }
}

run_ti_uot <- function(case) {
  id <- case_id(case)
  x <- case$inputs; p <- case$params; y <- case$outputs
  fit <- rfugw::ot_sinkhorn_unbalanced_ti(
    as_matrix(x$cost), as_vector(x$source), as_vector(x$target),
    epsilon = scalar(p$epsilon), rho = as_vector(p$rho),
    max_iter = integer_scalar(p$max_iter), tol = scalar(p$tol),
    plan = "dense"
  )
  a <- scalar(case$tolerance$atol); r <- scalar(case$tolerance$rtol)
  record_metric(id, "reference_plan", rfugw::rfugw_plan(fit),
                as_matrix(y$reference_plan), a, r)
  record_metric(id, "reference_objective", fit$regularized_objective,
                y$reference_objective, a, r)
  specialized_plan <- as_matrix(y$specialized_plan_monitor)
  record_metric(
    id, "specialized_pot_plan_monitor", rfugw::rfugw_plan(fit),
    specialized_plan, a, r, strict = FALSE, note = case$exception$reason
  )
  record_boolean(
    id, "specialized_pot_objective_not_better_than_certified_reference",
    scalar(y$specialized_objective_monitor) + 1e-12 >= fit$regularized_objective,
    strict = TRUE,
    note = "the monitored specialized POT path must not undercut the independently certified reference objective"
  )
  exceptions[[length(exceptions) + 1L]] <<- c(list(case = id), case$exception)
}

run_sinkhorn_divergence <- function(case) {
  id <- case_id(case)
  x <- case$inputs; p <- case$params; y <- case$outputs
  fit <- rfugw::ot_sinkhorn_divergence(
    source = as_matrix(x$source_points),
    target = as_matrix(x$target_points),
    p = integer_scalar(p$power),
    source_weights = as_vector(x$source),
    target_weights = as_vector(x$target),
    epsilon = scalar(p$epsilon), method = "log",
    max_iter = integer_scalar(p$max_iter), tol = scalar(p$tol)
  )
  a <- scalar(case$tolerance$atol); r <- scalar(case$tolerance$rtol)
  for (name in c("cross", "source_self", "target_self")) {
    record_metric(
      id, paste0(name, "_plan"), fit$component_solves[[name]]$plan,
      as_matrix(y$plans[[name]]), a, r
    )
    record_metric(
      id, paste0(name, "_value"), fit$component_values[[name]],
      y$component_values[[name]], a, r
    )
  }
  record_metric(id, "divergence", rfugw::rfugw_value(fit),
                y$divergence, a, r)
}

run_fugw <- function(case) {
  id <- case_id(case)
  x <- case$inputs; p <- case$params; y <- case$outputs
  fit <- rfugw::fugw_kl(
    Cx = as_matrix(x$source_cost), Cy = as_matrix(x$target_cost),
    wx = as_vector(x$source), wy = as_vector(x$target),
    reg_marginals = as_vector(p$rho), epsilon = scalar(p$epsilon),
    alpha = scalar(p$alpha), M = as_matrix(x$feature_cost),
    max_iter = integer_scalar(p$max_iter), tol = scalar(p$tol),
    max_iter_ot = integer_scalar(p$max_iter_ot), tol_ot = scalar(p$tol_ot),
    precision = "double"
  )
  a <- scalar(case$tolerance$atol); r <- scalar(case$tolerance$rtol)
  record_metric(id, "sample_plan", fit$pi_samp, as_matrix(y$sample_plan), a, r)
  record_metric(id, "feature_plan", fit$pi_feat, as_matrix(y$feature_plan), a, r)
  record_metric(id, "objective", fit$fugw_cost, y$objective, a, r)
}

run_ucoot <- function(case) {
  id <- case_id(case)
  x <- case$inputs; p <- case$params; y <- case$outputs
  fit <- rfugw::unbalanced_co_optimal_transport(
    X = as_matrix(x$source), Y = as_matrix(x$target),
    reg_marginals = as_vector(p$rho), epsilon = as_vector(p$epsilon),
    divergence = "kl", unbalanced_solver = "sinkhorn",
    alpha = c(0, 0), max_iter = integer_scalar(p$max_iter),
    tol = scalar(p$tol), max_iter_ot = integer_scalar(p$max_iter_ot),
    tol_ot = scalar(p$tol_ot), log = TRUE
  )
  a <- scalar(case$tolerance$atol); r <- scalar(case$tolerance$rtol)
  record_metric(id, "sample_plan", fit$pi_samp, as_matrix(y$sample_plan), a, r)
  record_metric(id, "feature_plan", fit$pi_feat, as_matrix(y$feature_plan), a, r)
  record_metric(id, "objective", fit$ucoot_cost, y$objective, a, r)
}

run_sampled_quality <- function(case) {
  id <- case_id(case)
  x <- case$inputs; p <- case$params; y <- case$outputs
  c1 <- as_matrix(x$source_cost); c2 <- as_matrix(x$target_cost)
  source <- as_vector(x$source); target <- as_vector(x$target)
  init <- as_matrix(x$init_plan); epsilon <- scalar(p$epsilon)
  max_iter <- integer_scalar(p$max_iter)
  random_states <- as_vector(p$random_states)
  dense <- rfugw::entropic_gromov_wasserstein(
    c1, c2, source, target, epsilon = epsilon, G0 = init,
    max_iter = max_iter, tol = 1e-9, solver = "PGD",
    sinkhorn_max_iter = 5000L, sinkhorn_tol = 1e-12,
    sinkhorn_method = "log", precision = "double"
  )
  sampled <- function(budget, random_state) {
    rfugw::sampled_gromov_wasserstein(
      c1, c2, source, target, nb_samples_grad = budget,
      epsilon = epsilon, max_iter = max_iter, random_state = random_state,
      sinkhorn_max_iter = 5000L, sinkhorn_tol = 1e-12, log = TRUE
    )
  }
  tiny <- lapply(random_states, function(state) {
    sampled(as_vector(p$tiny_budget), state)
  })
  high <- sampled(as_vector(p$rfugw_high_budget), random_states[[1L]])
  objective <- function(fit) rfugw::ot_gw_square(c1, c2, fit$plan)
  dense_objective <- objective(dense)
  tiny_objectives <- vapply(tiny, objective, numeric(1))
  high_objective <- objective(high)
  tiny_distances <- vapply(
    tiny, function(fit) sqrt(sum((fit$plan - dense$plan)^2)), numeric(1)
  )
  high_distance <- sqrt(sum((high$plan - dense$plan)^2))
  a <- scalar(case$tolerance$atol); r <- scalar(case$tolerance$rtol)
  record_metric(id, "dense_plan", dense$plan, as_matrix(y$dense_plan), a, r)
  record_metric(id, "dense_objective", dense_objective,
                y$independent_objectives$dense, a, r)
  record_metric(
    id, "tiny_plan_pot_monitor", tiny[[1L]]$plan,
    as_matrix(y$tiny_plans_monitor[[1L]]),
    0.2, 0.2, strict = FALSE,
    note = "stochastic implementations are compared by quality, not exact plan identity"
  )
  record_metric(
    id, "high_plan_pot_monitor", high$plan, as_matrix(y$high_plan_monitor),
    0.2, 0.2, strict = FALSE,
    note = "stochastic implementations are compared by quality, not exact plan identity"
  )
  record_boolean(id, "pot_high_budget_beats_median_tiny_objective_gap",
                 logical_scalar(y$pot_high_beats_median_tiny_objective_gap),
                 strict = FALSE,
                 note = case$exception$reason)
  record_boolean(id, "pot_high_budget_beats_median_tiny_plan_distance",
                 logical_scalar(y$pot_high_beats_median_tiny_plan_distance),
                 strict = FALSE,
                 note = case$exception$reason)
  record_boolean(
    id, "rfugw_high_budget_beats_median_tiny_objective_gap",
    abs(high_objective - dense_objective) <
      stats::median(abs(tiny_objectives - dense_objective)), strict = TRUE,
    note = sprintf("high_gap=%.8g median_tiny_gap=%.8g",
                   abs(high_objective - dense_objective),
                   stats::median(abs(tiny_objectives - dense_objective)))
  )
  record_boolean(
    id, "rfugw_high_budget_beats_median_tiny_plan_distance",
    high_distance < stats::median(tiny_distances), strict = TRUE,
    note = sprintf("high=%.8g median_tiny=%.8g", high_distance,
                   stats::median(tiny_distances))
  )
  record_metric(id, "tiny_row_marginal", rowSums(tiny[[1L]]$plan), source, 2e-5, 2e-5)
  record_metric(id, "tiny_col_marginal", colSums(tiny[[1L]]$plan), target, 2e-5, 2e-5)
  record_metric(id, "high_row_marginal", rowSums(high$plan), source, 2e-5, 2e-5)
  record_metric(id, "high_col_marginal", colSums(high$plan), target, 2e-5, 2e-5)
  exceptions[[length(exceptions) + 1L]] <<- c(list(case = id), case$exception)
}

run_barycenter_semantics <- function(case) {
  id <- case_id(case)
  x <- case$inputs; p <- case$params; y <- case$outputs
  costs <- lapply(x$costs, as_matrix)
  features <- lapply(x$features, as_matrix)
  weights <- lapply(x$weights, as_vector)
  barycenter_weights <- as_vector(x$barycenter_weights)
  lambdas <- as_vector(p$lambdas)
  common <- list(
    N = length(barycenter_weights), ps = weights, p = barycenter_weights,
    lambdas = lambdas, epsilon = scalar(p$epsilon),
    max_iter = integer_scalar(p$max_iter), tol = scalar(p$tol),
    sinkhorn_max_iter = integer_scalar(p$sinkhorn_max_iter),
    sinkhorn_tol = scalar(p$sinkhorn_tol), sinkhorn_method = "log",
    precision = "double", check_every = 1L, log = TRUE
  )
  gw <- do.call(
    rfugw::entropic_gromov_barycenters,
    c(common, list(Cs = costs, init_C = as_matrix(x$init_cost)))
  )
  fused <- do.call(
    rfugw::entropic_fused_gromov_barycenters,
    c(common, list(
      Ys = features, Cs = costs, alpha = scalar(p$alpha),
      init_C = as_matrix(x$init_cost), init_Y = as_matrix(x$init_features),
      feature_cost_metric = "sqeuclidean", feature_cost_normalize = FALSE
    ))
  )
  off_diagonal <- row(gw$C) != col(gw$C)
  a <- scalar(case$tolerance$atol); r <- scalar(case$tolerance$rtol)
  record_metric(id, "gw_off_diagonal", gw$C[off_diagonal],
                as_matrix(y$gw_cost)[off_diagonal], a, r)
  record_metric(id, "fused_off_diagonal", fused$C[off_diagonal],
                as_matrix(y$fused_cost)[off_diagonal], a, r)
  pot_feature_update <- Reduce(
    `+`,
    lapply(seq_along(features), function(index) {
      lambdas[[index]] *
        (as_matrix(y$fused_couplings[[index]]) %*% features[[index]])
    })
  )
  pot_feature_update <- sweep(
    pot_feature_update, 1L, 1 / barycenter_weights, `*`
  )
  record_metric(id, "fused_features_from_pot_couplings", fused$Y,
                pot_feature_update, a, r)
  record_metric(
    id, "pot_returned_fused_features_monitor", fused$Y,
    as_matrix(y$fused_features), a, r, strict = FALSE,
    note = case$exception$reason
  )
  for (index in seq_along(gw$couplings)) {
    record_metric(id, paste0("gw_coupling_", index), gw$couplings[[index]],
                  as_matrix(y$gw_couplings[[index]]), a, r)
    record_metric(id, paste0("fused_coupling_", index), fused$couplings[[index]],
                  as_matrix(y$fused_couplings[[index]]), a, r)
  }
  record_boolean(id, "rfugw_gw_diagonal_is_zero",
                 max(abs(diag(gw$C))) <= 1e-14, strict = TRUE)
  record_boolean(id, "rfugw_fused_diagonal_is_zero",
                 max(abs(diag(fused$C))) <= 1e-14, strict = TRUE)
  record_metric(
    id, "pot_gw_diagonal_monitor", diag(gw$C), diag(as_matrix(y$gw_cost)),
    a, r, strict = FALSE, note = case$exception$reason
  )
  record_metric(
    id, "pot_fused_diagonal_monitor", diag(fused$C), diag(as_matrix(y$fused_cost)),
    a, r, strict = FALSE, note = case$exception$reason
  )
  exceptions[[length(exceptions) + 1L]] <<- c(list(case = id), case$exception)
}

runners <- list(
  balanced_linear = run_balanced_linear,
  partial_linear = run_partial_linear,
  balanced_gromov = run_balanced_gromov,
  partial_gromov = run_partial_gromov,
  semirelaxed_gromov = run_semirelaxed,
  translation_invariant_uot = run_ti_uot,
  sinkhorn_divergence = run_sinkhorn_divergence,
  fugw = run_fugw,
  ucoot = run_ucoot,
  sampled_gromov_quality = run_sampled_quality,
  barycenter_semantics = run_barycenter_semantics
)

for (case in oracle$cases) {
  id <- case_id(case)
  runner <- runners[[case$family]]
  if (is.null(runner)) {
    record_boolean(id, "known_case_family", FALSE, strict = TRUE)
    next
  }
  check_diagnostics(case)
  tryCatch(
    runner(case),
    error = function(error) {
      record_boolean(
        id, "case_execution", FALSE, strict = TRUE,
        note = conditionMessage(error)
      )
    }
  )
}

strict_results <- Filter(function(x) isTRUE(x$strict), results)
strict_failures <- Filter(function(x) isTRUE(x$strict) && !isTRUE(x$ok), results)
receipt <- list(
  schema_version = "rfugw-pot-comparison-v1",
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  scope = scope,
  oracle = list(
    schema_version = oracle$schema_version,
    profile = oracle$profile,
    channel = oracle$channel,
    runtime = oracle$runtime,
    cases_sha256 = oracle$cases_sha256,
    case_count = oracle$case_count
  ),
  rfugw = list(
    version = as.character(utils::packageVersion("rfugw")),
    session = capture.output(utils::sessionInfo())
  ),
  metric_count = length(results),
  strict_metric_count = length(strict_results),
  strict_failure_count = length(strict_failures),
  all_strict_passed = !length(strict_failures),
  results = results,
  exceptions = exceptions
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(receipt, output_path, auto_unbox = TRUE, pretty = TRUE,
                     na = "null", null = "null")

cat(sprintf(
  "POT oracle comparison: cases=%d metrics=%d strict_failures=%d output=%s\n",
  length(oracle$cases), length(results), length(strict_failures), output_path
))
if (length(strict_failures)) {
  for (failure in strict_failures) {
    cat(sprintf(
      "FAIL %s :: %s (max_abs=%s, ratio=%s)%s\n",
      failure$case, failure$metric, format(failure$max_abs),
      format(failure$max_tolerance_ratio),
      if (is.null(failure$note)) "" else paste0(" -- ", failure$note)
    ))
  }
  quit(status = 1L)
}
cat("All strict live-POT differential checks passed.\n")
