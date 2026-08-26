# POT API extensions: partial GW/FGW, non-entropic semirelaxed GW/FGW,
# entropic barycenters wrappers, sampled/low-rank sampled GW, and
# unbalanced co-optimal/across-space divergences.

.clamp01 <- function(x) {
  min(1, max(0, x))
}

.parse_pair <- function(x, name) {
  if (length(x) == 1L) {
    x <- c(x, x)
  }
  if (length(x) != 2L || any(!is.finite(x))) {
    stop(sprintf("`%s` must have length 1 or 2 with finite numeric values.", name), call. = FALSE)
  }
  as.numeric(x)
}

.validate_mass <- function(m, p, q) {
  if (is.null(m)) {
    return(min(sum(p), sum(q)))
  }
  if (length(m) != 1L || !is.numeric(m) || !is.finite(m) || m < 0) {
    stop("`m` must be one finite number >= 0.", call. = FALSE)
  }
  mmax <- min(sum(p), sum(q))
  if (m > mmax + 1e-12) {
    stop("`m` must be <= min(sum(p), sum(q)).", call. = FALSE)
  }
  as.numeric(m)
}

.validate_partial_init <- function(G0, p, q, m, ns, nt) {
  if (is.null(G0)) {
    return((p %o% q) * (m / (sum(p) * sum(q))))
  }
  .assert_matrix(G0, "G0")
  if (nrow(G0) != ns || ncol(G0) != nt) {
    stop("`G0` must have shape nrow(C1) x nrow(C2).", call. = FALSE)
  }
  if (any(!is.finite(G0)) || any(G0 < 0)) {
    stop("`G0` must be finite and nonnegative.", call. = FALSE)
  }
  if (any(rowSums(G0) - p > 1e-8) || any(colSums(G0) - q > 1e-8)) {
    stop("`G0` must satisfy rowSums(G0) <= p and colSums(G0) <= q.", call. = FALSE)
  }
  if (abs(sum(G0) - m) > 1e-6) {
    stop("`sum(G0)` must equal `m` (within tolerance).", call. = FALSE)
  }
  G0
}

.partial_penalty <- function(cost) {
  ma <- suppressWarnings(max(cost, na.rm = TRUE))
  aa <- suppressWarnings(max(abs(cost), na.rm = TRUE))
  if (!is.finite(ma)) ma <- 0
  if (!is.finite(aa) || aa <= 0) aa <- 1
  ma + aa + 1
}

.make_partial_transport_lp_solver <- function(
    a,
    b,
    m,
    nb_dummies = 1L,
    lp_solver = c("cpp_transport", "lp_matrix", "lp_transport"),
    lp_scale = 1e6) {
  lp_solver <- match.arg(lp_solver)
  if (!identical(lp_solver, "cpp_transport") && !requireNamespace("lpSolve", quietly = TRUE)) {
    stop(
      paste(
        "`lp_solver = \"lp_transport\"` and `\"lp_matrix\"` require package `lpSolve`.",
        "Install it with install.packages(\"lpSolve\"), or use the default",
        "`lp_solver = \"cpp_transport\"` which has no extra dependency.",
        "Entropic partial solvers do not need lpSolve."
      ),
      call. = FALSE
    )
  }
  nb_dummies <- .validate_count(nb_dummies, "nb_dummies")

  ns <- length(a)
  nt <- length(b)
  a_ext <- c(a, rep((sum(b) - m) / nb_dummies, nb_dummies))
  b_ext <- c(b, rep((sum(a) - m) / nb_dummies, nb_dummies))
  total <- sum(a_ext)
  if (total <= 0) {
    stop("Infeasible partial transport setup: extended total mass must be positive.", call. = FALSE)
  }

  solve_ext <- .make_transport_lp_solver(
    p = a_ext / total,
    q = b_ext / total,
    scale = lp_scale,
    solver = lp_solver
  )

  function(cost) {
    .assert_matrix(cost, "cost")
    if (nrow(cost) != ns || ncol(cost) != nt) {
      stop("`cost` has incompatible shape for partial LP direction step.", call. = FALSE)
    }
    penalty <- .partial_penalty(cost)
    cost_ext <- matrix(0, nrow = ns + nb_dummies, ncol = nt + nb_dummies)
    cost_ext[seq_len(ns), seq_len(nt)] <- cost
    cost_ext[ns + seq_len(nb_dummies), nt + seq_len(nb_dummies)] <- penalty

    G_ext <- solve_ext(cost_ext) * total
    G_ext[seq_len(ns), seq_len(nt), drop = FALSE]
  }
}

.gw_square_terms <- function(C1, C2, G, symmetric = TRUE) {
  cpp_gw_square_terms_square(C1, C2, G, symmetric = isTRUE(symmetric))
}

.partial_fgw_cg_core <- function(
    M,
    C1,
    C2,
    p,
    q,
    m,
    alpha,
    G0,
    max_iter,
    tol,
    symmetric,
    lp_solver,
    lp_scale,
    nb_dummies,
    verbose = FALSE) {
  lin_w <- 1 - alpha
  quad_w <- alpha

  solve_direction <- .make_partial_transport_lp_solver(
    a = p,
    b = q,
    m = m,
    nb_dummies = nb_dummies,
    lp_solver = lp_solver,
    lp_scale = lp_scale
  )

  gw <- .gw_square_terms(C1, C2, G0, symmetric = symmetric)
  cost <- lin_w * sum(M * G0) + quad_w * gw$loss
  grad <- gw$grad
  loss_trace <- numeric(max_iter + 1L)
  loss_trace[1] <- cost

  it <- 0L
  rel_delta <- Inf
  abs_delta <- Inf

  for (k in seq_len(max_iter)) {
    Mi <- lin_w * M + quad_w * grad
    Gc <- solve_direction(Mi)
    delta <- Gc - G0

    gw_c <- .gw_square_terms(C1, C2, Gc, symmetric = symmetric)
    grad_delta <- gw_c$grad - grad

    a_ls <- quad_w * 0.5 * sum(grad_delta * delta)
    b_ls <- lin_w * sum(M * delta) + quad_w * sum(grad * delta)
    step <- .clamp01(.solve_1d_linesearch_quad(a_ls, b_ls))

    new_cost <- cost + a_ls * step * step + b_ls * step
    G0 <- G0 + step * delta
    grad <- grad + step * grad_delta

    abs_delta <- abs(new_cost - cost)
    rel_delta <- abs_delta / (abs(new_cost) + 1e-15)
    cost <- new_cost
    loss_trace[k + 1L] <- cost
    it <- k

    if (isTRUE(verbose) && (k %% 25L == 0L || k == 1L)) {
      cat(sprintf("iter=%d cost=%.8e rel=%.3e abs=%.3e\n", k, cost, rel_delta, abs_delta))
    }

    if (rel_delta <= tol) {
      break
    }
  }

  gw_end <- .gw_square_terms(C1, C2, G0, symmetric = symmetric)
  lin_loss <- lin_w * sum(M * G0)
  quad_loss <- quad_w * gw_end$loss

  list(
    plan = G0,
    objective = lin_loss + quad_loss,
    lin_loss = lin_loss,
    quad_loss = quad_loss,
    gw_loss = gw_end$loss,
    iterations = as.integer(it),
    error = as.numeric(rel_delta),
    abs_error = as.numeric(abs_delta),
    loss_trace = loss_trace[seq_len(it + 1L)]
  )
}

.validate_variable_partial_init <- function(G0, p, q, ns, nt) {
  if (is.null(G0)) {
    return(matrix(0, ns, nt))
  }
  G0 <- .validate_finite_matrix(G0, "G0")
  if (nrow(G0) != ns || ncol(G0) != nt) {
    stop("`G0` must have shape nrow(C1) x nrow(C2).", call. = FALSE)
  }
  if (any(G0 < 0)) {
    stop("`G0` must be nonnegative.", call. = FALSE)
  }
  scale <- max(1, sum(p), sum(q))
  capacity_tol <- 1e-10 * scale
  if (any(rowSums(G0) - p > capacity_tol) ||
      any(colSums(G0) - q > capacity_tol)) {
    stop(
      "`G0` must satisfy rowSums(G0) <= p and colSums(G0) <= q.",
      call. = FALSE
    )
  }
  G0
}

.penalized_partial_fgw_terms <- function(
    M, C1, C2, G, p, q, alpha, discard_penalty, symmetric) {
  gw <- .gw_square_terms(C1, C2, G, symmetric = symmetric)
  feature_raw <- sum(M * G)
  feature_term <- (1 - alpha) * feature_raw
  structure_term <- alpha * gw$loss
  transported_mass <- sum(G)
  discarded_source <- sum(p) - transported_mass
  discarded_target <- sum(q) - transported_mass
  penalty_term <- discard_penalty * (discarded_source + discarded_target)
  list(
    objective = feature_term + structure_term + penalty_term,
    feature_raw = feature_raw,
    feature_term = feature_term,
    structure_raw = gw$loss,
    structure_term = structure_term,
    penalty_term = penalty_term,
    transported_mass = transported_mass,
    discarded_source_mass = discarded_source,
    discarded_target_mass = discarded_target,
    gw_gradient = gw$grad,
    gradient = (1 - alpha) * M + alpha * gw$grad -
      2 * discard_penalty
  )
}

.solve_penalized_partial_fgw_direction <- function(
    linear_cost, p, q, discard_penalty, max_iter, tol) {
  cost_shift <- max(0, -min(linear_cost))
  adjusted_penalty <- discard_penalty + cost_shift / 2
  out <- ot_partial_penalized(
    M = linear_cost + cost_shift,
    p = p,
    q = q,
    discard_penalty = adjusted_penalty,
    max_iter = max_iter,
    tol = tol
  )
  out$linear_cost_shift <- cost_shift
  out$adjusted_discard_penalty <- adjusted_penalty
  out$linearized_objective <- sum(linear_cost * out$plan) -
    2 * discard_penalty * sum(out$plan)
  out
}

.penalized_partial_fgw_line_search <- function(
    M, C1, C2, G, direction, alpha, discard_penalty, symmetric,
    current_terms = NULL) {
  if (is.null(current_terms)) {
    stop("`current_terms` is required for certified line search.", call. = FALSE)
  }
  delta <- direction - G
  direction_gw <- .gw_square_terms(
    C1, C2, direction, symmetric = symmetric
  )
  gradient_delta <- direction_gw$grad - current_terms$gw_gradient
  quadratic <- alpha * 0.5 * sum(gradient_delta * delta)
  linear <- (1 - alpha) * sum(M * delta) +
    alpha * sum(current_terms$gw_gradient * delta) -
    2 * discard_penalty * sum(delta)
  step <- .clamp01(.solve_1d_linesearch_quad(quadratic, linear))
  list(
    step = step,
    quadratic = quadratic,
    linear = linear,
    predicted_objective = current_terms$objective +
      quadratic * step^2 + linear * step,
    delta = delta
  )
}

.variable_partial_feasibility <- function(G, p, q) {
  scale <- max(1, sum(p), sum(q))
  tolerance <- max(1e-10 * scale, 100 * .Machine$double.eps * scale)
  row_residual <- max(c(rowSums(G) - p, 0))
  col_residual <- max(c(colSums(G) - q, 0))
  nonnegative_residual <- max(c(-G, 0))
  mass_residual <- max(c(-sum(G), sum(G) - min(sum(p), sum(q)), 0))
  residual <- max(
    row_residual, col_residual, nonnegative_residual, mass_residual
  )
  list(
    residual = residual,
    tolerance = tolerance,
    row_residual = row_residual,
    col_residual = col_residual,
    nonnegative_residual = nonnegative_residual,
    mass_residual = mass_residual,
    feasible = is.finite(residual) && residual <= tolerance
  )
}

.entropic_partial_wasserstein <- function(
    a,
    b,
    M,
    reg,
    m,
    numItermax = 1000L,
    stopThr = 1e-100,
    verbose = FALSE,
    log = FALSE,
    method = c("auto", "scaling", "log")) {
  method <- match.arg(method)
  fit <- ot_partial_sinkhorn(
    M = M,
    p = a,
    q = b,
    mass = m,
    epsilon = reg,
    method = method,
    max_iter = numItermax,
    tol = stopThr,
    check_every = 1L,
    verbose = verbose
  )
  if (isTRUE(log)) {
    return(list(
      plan = fit$plan,
      log = list(
        err = fit$convergence_trace$update_residual,
        partial_w_dist = fit$ot_dist
      ),
      result = fit
    ))
  }
  plan <- fit$plan
  attr(plan, "partial_sinkhorn_result") <- fit
  plan
}

.partial_fgw_entropic_core <- function(
    M,
    C1,
    C2,
    p,
    q,
    m,
    reg,
    alpha,
    G0,
    max_iter,
    tol,
    symmetric,
    inner_max_iter = 300L,
    inner_tol = 1e-12,
    verbose = FALSE,
    check_every = 10L,
    method = c("auto", "scaling", "log")) {
  method <- match.arg(method)
  lin_w <- 1 - alpha
  quad_w <- alpha

  G <- G0
  err <- Inf
  it <- 0L
  err_trace <- numeric()
  all_inner_converged <- TRUE
  total_inner_iterations <- 0L
  max_inner_residual <- 0
  final_inner_residual <- Inf
  final_inner_status <- "not_run"
  inner_effective_methods <- character()
  inner_dynamic_ranges <- numeric()
  inner_transitions <- character()

  for (k in seq_len(max_iter)) {
    Gprev <- G
    gw <- .gw_square_terms(C1, C2, G, symmetric = symmetric)
    M_entr <- quad_w * gw$grad + lin_w * M

    inner <- ot_partial_sinkhorn(
      M = M_entr, p = p, q = q, mass = m, epsilon = reg,
      method = method, max_iter = inner_max_iter, tol = inner_tol,
      check_every = 1L, verbose = FALSE
    )
    G <- inner$plan
    total_inner_iterations <- total_inner_iterations + inner$iterations
    final_inner_residual <- inner$residual
    max_inner_residual <- max(max_inner_residual, inner$residual)
    final_inner_status <- inner$status
    inner_effective_methods <- c(
      inner_effective_methods, inner$effective_sinkhorn_method
    )
    inner_dynamic_ranges <- c(
      inner_dynamic_ranges, inner$sinkhorn_dynamic_range
    )
    inner_transitions <- c(
      inner_transitions, inner$sinkhorn_backend_transition
    )
    if (!isTRUE(inner$converged)) {
      all_inner_converged <- FALSE
      it <- k
      break
    }

    if ((k %% check_every) == 0L || k == 1L) {
      err <- sqrt(sum((G - Gprev)^2))
      err_trace <- c(err_trace, err)
      if (isTRUE(verbose) && (k %% 50L == 0L || k == 1L)) {
        obj <- lin_w * sum(M * G) + quad_w * .gw_square_terms(C1, C2, G, symmetric = symmetric)$loss
        cat(sprintf("iter=%d err=%.8e obj=%.8e\n", k, err, obj))
      }
      if (isTRUE(err <= tol)) {
        it <- k
        break
      }
    }
    it <- k
  }

  gw_end <- .gw_square_terms(C1, C2, G, symmetric = symmetric)
  lin_loss <- lin_w * sum(M * G)
  quad_loss <- quad_w * gw_end$loss

  list(
    plan = G,
    objective = lin_loss + quad_loss,
    lin_loss = lin_loss,
    quad_loss = quad_loss,
    gw_loss = gw_end$loss,
    iterations = as.integer(it),
    error = as.numeric(err),
    err_trace = err_trace,
    lp_ok = all_inner_converged,
    inner_converged = all_inner_converged,
    inner_iterations = as.integer(total_inner_iterations),
    inner_residual = final_inner_residual,
    max_inner_residual = max_inner_residual,
    inner_status = final_inner_status,
    inner_termination_reason = if (all_inner_converged) {
      "all_partial_sinkhorn_solves_certified"
    } else {
      "partial_sinkhorn_failure"
    },
    sinkhorn_dispatch = list(
      requested = method,
      effective = if (length(unique(inner_effective_methods)) == 1L) {
        unique(inner_effective_methods)
      } else {
        paste(unique(inner_effective_methods), collapse = "+")
      },
      metric = if (length(inner_dynamic_ranges)) {
        max(inner_dynamic_ranges)
      } else {
        NA_real_
      },
      threshold = 100,
      reason = "certified_public_partial_sinkhorn_per_outer_iteration",
      transition = if (length(unique(inner_transitions)) == 1L) {
        unique(inner_transitions)
      } else {
        paste(unique(inner_transitions), collapse = "+")
      }
    )
  )
}

.partial_fgw_exact_dispatch <- function(
    M,
    C1,
    C2,
    p,
    q,
    m,
    alpha,
    G0,
    max_iter,
    tol,
    symmetric,
    lp_solver,
    lp_scale,
    nb_dummies,
    verbose = FALSE) {
  if (identical(lp_solver, "cpp_transport") &&
      exists("cpp_partial_fgw_exact_square", mode = "function")) {
    return(cpp_partial_fgw_exact_square(
      M = M,
      C1 = C1,
      C2 = C2,
      p = p,
      q = q,
      m = m,
      alpha = alpha,
      symmetric = isTRUE(symmetric),
      init_plan = if (is.null(G0)) .empty_feature_cost() else G0,
      max_iter = as.integer(max_iter),
      tol = tol,
      nb_dummies = as.integer(nb_dummies),
      lp_max_iter = 20000L,
      lp_tol = 1e-12
    ))
  }
  M_r <- if (length(M) == 0L) matrix(0, nrow(C1), nrow(C2)) else M
  out <- .partial_fgw_cg_core(
    M = M_r,
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    m = m,
    alpha = alpha,
    G0 = G0,
    max_iter = as.integer(max_iter),
    tol = tol,
    symmetric = symmetric,
    lp_solver = lp_solver,
    lp_scale = lp_scale,
    nb_dummies = as.integer(nb_dummies),
    verbose = verbose
  )
  out$lp_ok <- TRUE
  out$inner_converged <- TRUE
  out$inner_iterations <- NA_integer_
  out$inner_residual <- 0
  out$max_inner_residual <- 0
  out$inner_status <- "optimal"
  out$inner_termination_reason <- "optimal"
  out
}

.partial_fgw_entropic_dispatch <- function(
    M,
    C1,
    C2,
    p,
    q,
    m,
    reg,
    alpha,
    G0,
    max_iter,
    tol,
    symmetric,
    inner_max_iter,
    inner_tol,
    verbose = FALSE,
    check_every = 10L,
    method = c("auto", "scaling", "log")) {
  method <- match.arg(method)
  M_r <- if (length(M) == 0L) matrix(0, nrow(C1), nrow(C2)) else M
  .partial_fgw_entropic_core(
    M = M_r,
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    m = m,
    reg = reg,
    alpha = alpha,
    G0 = G0,
    max_iter = as.integer(max_iter),
    tol = tol,
    symmetric = symmetric,
    inner_max_iter = inner_max_iter,
    inner_tol = inner_tol,
    verbose = verbose,
    check_every = check_every,
    method = method
  )
}

.partial_log_with_nested_status <- function(
    out,
    value_field,
    tol,
    max_iter,
    p,
    q,
    mass_target,
    feasibility_tol,
    objective_recomputed,
    objective_components = list(),
    mass_defaulted = FALSE) {
  out[[value_field]] <- out$objective
  out$objective <- NULL
  outer_converged <- is.finite(out$error) && out$error <= tol
  ans <- .attach_solver_diagnostics(
    out,
    residual = out$error,
    converged = outer_converged,
    iterations = out$iterations,
    max_iter = max_iter,
    p = p,
    q = q,
    plan = out$plan,
    inner_residual = out$inner_residual %||% NA_real_,
    max_inner_residual = out$max_inner_residual %||% NA_real_,
    inner_iterations = out$inner_iterations %||% NA_integer_,
    inner_converged = out$inner_converged %||% NA,
    inner_status = out$inner_status,
    feasibility = "partial",
    feasibility_tol = feasibility_tol,
    mass_target = mass_target,
    objective_recomputed = objective_recomputed,
    objective_components = objective_components,
    lp_ok = out$lp_ok %||% TRUE
  )
  if (!is.na(ans$inner_converged) && !isTRUE(ans$inner_converged)) {
    ans$status <- "inner_failure"
    ans$converged <- FALSE
    ans$warning_payload <- list(
      code = "inner_failure",
      message = "Solver result is not certified: inner_failure."
    )
  }
  ans$transported_mass_target <- mass_target
  ans$transported_mass_defaulted <- isTRUE(mass_defaulted)
  ans$termination_reason <- .termination_reason_from_result(ans, max_iter)
  class(ans) <- setdiff(class(ans), "rfugw_result")
  ans
}

.attach_partial_sinkhorn_dispatch <- function(out, dispatch) {
  out$requested_sinkhorn_method <- dispatch$requested
  out$effective_sinkhorn_method <- dispatch$effective
  out$sinkhorn_backend_transition <- dispatch$transition
  out$sinkhorn_dispatch_reason <- dispatch$reason
  out$sinkhorn_dynamic_range <- dispatch$metric
  out$sinkhorn_scaling_threshold <- dispatch$threshold
  out
}

.sr_row_min_direction <- function(Mi, p) {
  if (exists("cpp_srfgw_row_min_direction", mode = "function")) {
    return(cpp_srfgw_row_min_direction(Mi, p))
  }
  row_min <- apply(Mi, 1L, min)
  mask <- Mi <= (row_min + 1e-15)
  denom <- rowSums(mask)
  denom[denom <= 0] <- 1
  mask * (p / denom)
}

.srfgw_square_terms <- function(C1, C2, G, p, symmetric = TRUE) {
  ns <- nrow(C1)
  nt <- nrow(C2)
  qG <- colSums(G)

  ones_p <- rep(1, ns)
  fC2t <- t(C2^2)

  constC <- tcrossprod(as.vector((C1^2) %*% p), rep(1, nt))
  marginal <- tcrossprod(ones_p, as.vector(qG %*% fC2t))
  tens <- constC + marginal - C1 %*% G %*% t(2 * C2)

  quad <- sum(tens * G)
  grad <- 2 * tens

  if (!symmetric) {
    C1t <- t(C1)
    C2t <- t(C2)
    constCt <- tcrossprod(as.vector((C1t^2) %*% p), rep(1, nt))
    marginal2 <- tcrossprod(ones_p, as.vector(qG %*% (C2^2)))
    tenst <- constCt + marginal2 - C1t %*% G %*% t(2 * C2t)

    quad <- 0.5 * (quad + sum(tenst * G))
    grad <- 0.5 * (grad + 2 * tenst)
  }

  list(quad = quad, grad = grad)
}

.semirelaxed_fgw_cg_core <- function(
    M,
    C1,
    C2,
    p,
    alpha,
    symmetric,
    G0,
    max_iter,
    tol_rel,
    tol_abs,
    verbose = FALSE) {
  ns <- nrow(C1)
  nt <- nrow(C2)
  lin_w <- 1 - alpha

  if (is.null(G0)) {
    G <- p %o% rep(1 / nt, nt)
  } else {
    G <- .validate_semirelaxed_init(G0, p, ns, nt)
  }

  state <- .srfgw_square_terms(C1, C2, G, p, symmetric = symmetric)
  quad_raw <- state$quad
  cost <- lin_w * sum(M * G) + alpha * quad_raw

  rel_delta <- Inf
  abs_delta <- Inf
  loss_trace <- numeric(max_iter + 1L)
  loss_trace[1] <- cost
  it <- 0L

  ones_p <- rep(1, ns)
  hC1 <- C1
  hC2 <- 2 * C2
  fC2t <- t(C2^2)

  for (k in seq_len(max_iter)) {
    Mi <- lin_w * M + alpha * state$grad
    Gc <- .sr_row_min_direction(Mi, p)
    delta <- Gc - G

    qG <- colSums(G)
    qdelta <- colSums(delta)

    dot <- hC1 %*% delta %*% t(hC2)
    dot_qG <- tcrossprod(ones_p, as.vector(qG %*% fC2t))
    dot_qdelta <- tcrossprod(ones_p, as.vector(qdelta %*% fC2t))

    a_ls <- alpha * sum((dot_qdelta - dot) * delta)
    b_ls <- sum((lin_w * M) * delta) + alpha * (
      sum((dot_qdelta - dot) * G) +
        sum((dot_qG - hC1 %*% G %*% t(hC2)) * delta)
    )

    step <- .clamp01(.solve_1d_linesearch_quad(a_ls, b_ls))
    new_cost <- cost + a_ls * step * step + b_ls * step

    G <- G + step * delta
    state <- .srfgw_square_terms(C1, C2, G, p, symmetric = symmetric)

    abs_delta <- abs(new_cost - cost)
    rel_delta <- abs_delta / (abs(new_cost) + 1e-15)
    cost <- new_cost
    loss_trace[k + 1L] <- cost
    it <- k

    if (isTRUE(verbose) && (k %% 25L == 0L || k == 1L)) {
      cat(sprintf("iter=%d cost=%.8e rel=%.3e abs=%.3e\n", k, cost, rel_delta, abs_delta))
    }

    if ((is.finite(rel_delta) && rel_delta <= tol_rel) || (is.finite(abs_delta) && abs_delta <= tol_abs)) {
      break
    }
  }

  q <- colSums(G)
  quad_raw <- state$quad
  lin_loss <- lin_w * sum(M * G)
  quad_loss <- alpha * quad_raw

  list(
    plan = G,
    q = q,
    lin_loss = lin_loss,
    quad_loss = quad_loss,
    srfgw_dist = lin_loss + quad_loss,
    srgw_dist = quad_raw,
    iterations = as.integer(it),
    error = as.numeric(rel_delta),
    abs_error = as.numeric(abs_delta),
    loss_trace = loss_trace[seq_len(it + 1L)],
    symmetric = symmetric
  )
}

.semirelaxed_fgw_exact_dispatch <- function(
    M,
    C1,
    C2,
    p,
    alpha,
    symmetric,
    G0,
    max_iter,
    tol_rel,
    tol_abs,
    verbose = FALSE) {
  ns <- nrow(C1)
  nt <- nrow(C2)
  G_init <- if (is.null(G0)) {
    .empty_feature_cost()
  } else {
    .validate_semirelaxed_init(G0, p, ns, nt)
  }
  use_mixed_precision <- (ns * nt >= 4000L) &&
    .runtime_env_bool("RFUGW_SEMIRELAXED_MIXED", TRUE)
  if (isTRUE(symmetric) && exists("cpp_semirelaxed_fgw_cg_square_fast", mode = "function")) {
    return(cpp_semirelaxed_fgw_cg_square_fast(
      M = M,
      C1 = C1,
      C2 = C2,
      p = p,
      alpha = alpha,
      init_plan = G_init,
      max_iter = as.integer(max_iter),
      tol_rel = tol_rel,
      tol_abs = tol_abs,
      verbose = isTRUE(verbose),
      use_mixed_precision = use_mixed_precision
    ))
  }
  if (exists("cpp_semirelaxed_fgw_exact_square", mode = "function")) {
    return(cpp_semirelaxed_fgw_exact_square(
      M = M,
      C1 = C1,
      C2 = C2,
      p = p,
      alpha = alpha,
      symmetric = isTRUE(symmetric),
      init_plan = G_init,
      max_iter = as.integer(max_iter),
      tol_rel = tol_rel,
      tol_abs = tol_abs
    ))
  }
  M_r <- if (length(M) == 0L) matrix(0, ns, nt) else M
  .semirelaxed_fgw_cg_core(
    M = M_r,
    C1 = C1,
    C2 = C2,
    p = p,
    alpha = alpha,
    symmetric = symmetric,
    G0 = G0,
    max_iter = as.integer(max_iter),
    tol_rel = tol_rel,
    tol_abs = tol_abs,
    verbose = verbose
  )
}

.gw_square_value <- function(C1, C2, G, p, q, symmetric = TRUE) {
  ns <- nrow(C1)
  nt <- nrow(C2)
  constC <- tcrossprod(as.vector((C1^2) %*% p), rep(1, nt)) +
    tcrossprod(rep(1, ns), as.vector(q %*% t(C2^2)))
  tens <- constC - C1 %*% G %*% t(2 * C2)
  val <- sum(tens * G)
  if (!symmetric) {
    constCt <- tcrossprod(as.vector((t(C1)^2) %*% p), rep(1, nt)) +
      tcrossprod(rep(1, ns), as.vector(q %*% (C2^2)))
    tenst <- constCt - t(C1) %*% G %*% t(2 * t(C2))
    val <- 0.5 * (val + sum(tenst * G))
  }
  val
}

.sinkhorn_balanced <- function(a, b, M, reg, max_iter = 1000L, tol = 1e-9) {
  K <- exp(-M / reg)
  u <- rep(1, length(a))
  v <- rep(1, length(b))
  err <- Inf
  for (it in seq_len(as.integer(max_iter))) {
    u_prev <- u
    Kv <- as.vector(K %*% v)
    Kv[Kv <= 0] <- 1e-300
    u <- a / Kv
    Ktu <- as.vector(t(K) %*% u)
    Ktu[Ktu <= 0] <- 1e-300
    v <- b / Ktu

    if ((it %% 10L) == 0L || it == 1L) {
      err <- max(abs(u - u_prev))
      if (err <= tol) break
    }
  }
  (u %o% v) * K
}

.kl_div_mass <- function(x, y, mass = TRUE) {
  tiny <- 1e-300
  out <- sum(x * (log(pmax(x, tiny)) - log(pmax(y, tiny))))
  if (isTRUE(mass)) {
    out <- out - sum(x) + sum(y)
  }
  out
}

.div_to_product_kl <- function(pi, a, b, pi1 = NULL, pi2 = NULL, mass = TRUE) {
  if (is.null(pi1)) pi1 <- rowSums(pi)
  if (is.null(pi2)) pi2 <- colSums(pi)
  tiny <- 1e-300
  res <- sum(pi * log(pmax(pi, tiny))) -
    sum(pi1 * log(pmax(a, tiny))) -
    sum(pi2 * log(pmax(b, tiny)))
  if (isTRUE(mass)) {
    res <- res - sum(pi1) + sum(a) * sum(b)
  }
  res
}

.div_between_product_kl <- function(mu, nu, alpha, beta) {
  m_mu <- sum(mu)
  m_nu <- sum(nu)
  m_alpha <- sum(alpha)
  m_beta <- sum(beta)
  const <- (m_mu - m_alpha) * (m_nu - m_beta)
  m_nu * .kl_div_mass(mu, alpha, mass = TRUE) +
    m_mu * .kl_div_mass(nu, beta, mass = TRUE) +
    const
}

.uot_cost_matrix_kl <- function(data, pi, tuple_p, hyperparams, reg_type) {
  X_sqr <- data$X_sqr
  Y_sqr <- data$Y_sqr
  X <- data$X
  Y <- data$Y
  Y_t <- data$Y_t
  M <- data$M

  rho_x <- hyperparams[[1]]
  rho_y <- hyperparams[[2]]
  eps <- hyperparams[[3]]
  a <- tuple_p[[1]]
  b <- tuple_p[[2]]

  pi1 <- rowSums(pi)
  pi2 <- colSums(pi)
  A <- as.vector(X_sqr %*% pi1)
  B <- as.vector(Y_sqr %*% pi2)
  uot_cost <- tcrossprod(A, rep(1, length(B))) + tcrossprod(rep(1, length(A)), B) - 2 * (X %*% pi %*% Y_t)

  if (!is.null(M)) {
    uot_cost <- uot_cost + M
  }

  if (is.finite(rho_x) && rho_x != 0) {
    uot_cost <- uot_cost + rho_x * .kl_div_mass(pi1, a, mass = FALSE)
  }
  if (is.finite(rho_y) && rho_y != 0) {
    uot_cost <- uot_cost + rho_y * .kl_div_mass(pi2, b, mass = FALSE)
  }
  if (identical(reg_type, "joint") && eps > 0) {
    uot_cost <- uot_cost + eps * .div_to_product_kl(pi, a, b, pi1 = pi1, pi2 = pi2, mass = FALSE)
  }

  uot_cost
}

.sinkhorn_unbalanced_kl <- function(
    M,
    a,
    b,
    reg,
    rho,
    c = NULL,
    max_iter = 500L,
    tol = 1e-7,
    plan_init = NULL) {
  if (!is.finite(reg) || reg <= 0) {
    stop("KL Sinkhorn unbalanced solver requires `reg > 0`.", call. = FALSE)
  }
  rho_x <- rho[[1]]
  rho_y <- rho[[2]]
  tau_x <- if (is.infinite(rho_x)) 1 else if (rho_x <= 0) 0 else rho_x / (rho_x + reg)
  tau_y <- if (is.infinite(rho_y)) 1 else if (rho_y <= 0) 0 else rho_y / (rho_y + reg)

  if (is.null(c)) {
    c <- a %o% b
  }
  K <- c * exp(-M / reg)
  K[K <= 0] <- 1e-300

  if (!is.null(plan_init)) {
    .assert_matrix(plan_init, "plan_init")
    if (!all(dim(plan_init) == dim(M))) {
      stop("`plan_init` has incompatible shape.", call. = FALSE)
    }
    pi0 <- pmax(plan_init, 1e-300)
    u <- rowSums(pi0) / pmax(rowSums(K), 1e-300)
    v <- colSums(pi0) / pmax(colSums(K), 1e-300)
    u[!is.finite(u)] <- 1
    v[!is.finite(v)] <- 1
  } else {
    u <- rep(1, nrow(M))
    v <- rep(1, ncol(M))
  }

  err <- Inf
  it <- 0L
  for (k in seq_len(as.integer(max_iter))) {
    u_prev <- u

    Kv <- as.vector(K %*% v)
    Kv[Kv <= 0] <- 1e-300
    if (tau_x == 0) {
      u <- rep(1, length(a))
    } else {
      u <- (a / Kv)^tau_x
    }

    Ktu <- as.vector(crossprod(K, u))
    Ktu[Ktu <= 0] <- 1e-300
    if (tau_y == 0) {
      v <- rep(1, length(b))
    } else {
      v <- (b / Ktu)^tau_y
    }

    if ((k %% 5L) == 0L || k == 1L) {
      err <- max(abs(u - u_prev))
      if (err <= tol) {
        it <- k
        break
      }
    }
    it <- k
  }

  plan <- (u %o% v) * K
  list(plan = plan, iterations = as.integer(it), error = as.numeric(err), potentials = list(u = u, v = v))
}

.fused_unbalanced_across_spaces_cost_kl <- function(
    M_linear,
    data,
    tuple_pxy_samp,
    tuple_pxy_feat,
    pi_samp,
    pi_feat,
    hyperparams,
    reg_type) {
  rho_x <- hyperparams[[1]]
  rho_y <- hyperparams[[2]]
  eps_samp <- hyperparams[[3]]
  eps_feat <- hyperparams[[4]]

  M_samp <- M_linear[[1]]
  M_feat <- M_linear[[2]]
  px_samp <- tuple_pxy_samp[[1]]
  py_samp <- tuple_pxy_samp[[2]]
  pxy_samp <- tuple_pxy_samp[[3]]
  px_feat <- tuple_pxy_feat[[1]]
  py_feat <- tuple_pxy_feat[[2]]
  pxy_feat <- tuple_pxy_feat[[3]]

  X_sqr <- data[[1]]
  Y_sqr <- data[[2]]
  X <- data[[3]]
  Y <- data[[4]]

  pi1_samp <- rowSums(pi_samp)
  pi2_samp <- colSums(pi_samp)
  pi1_feat <- rowSums(pi_feat)
  pi2_feat <- colSums(pi_feat)

  A_sqr <- sum((X_sqr %*% pi1_feat) * pi1_samp)
  B_sqr <- sum((Y_sqr %*% pi2_feat) * pi2_samp)
  AB <- (X %*% pi_feat %*% t(Y)) * pi_samp
  linear_cost <- A_sqr + B_sqr - 2 * sum(AB)

  ucoot_cost <- linear_cost
  if (!is.null(M_samp)) ucoot_cost <- ucoot_cost + sum(pi_samp * M_samp)
  if (!is.null(M_feat)) ucoot_cost <- ucoot_cost + sum(pi_feat * M_feat)

  if (is.finite(rho_x) && rho_x != 0) {
    ucoot_cost <- ucoot_cost + rho_x * .div_between_product_kl(pi1_samp, pi1_feat, px_samp, px_feat)
  }
  if (is.finite(rho_y) && rho_y != 0) {
    ucoot_cost <- ucoot_cost + rho_y * .div_between_product_kl(pi2_samp, pi2_feat, py_samp, py_feat)
  }

  if (identical(reg_type, "joint")) {
    if (eps_samp != 0) {
      ucoot_cost <- ucoot_cost + eps_samp * .div_between_product_kl(pi_samp, pi_feat, pxy_samp, pxy_feat)
    }
  } else {
    if (eps_samp != 0) {
      ucoot_cost <- ucoot_cost + eps_samp * .div_to_product_kl(
        pi_samp,
        px_samp,
        py_samp,
        pi1_samp,
        pi2_samp,
        mass = TRUE
      )
    }
    if (eps_feat != 0) {
      ucoot_cost <- ucoot_cost + eps_feat * .div_to_product_kl(
        pi_feat,
        px_feat,
        py_feat,
        pi1_feat,
        pi2_feat,
        mass = TRUE
      )
    }
  }

  list(linear_cost = linear_cost, ucoot_cost = ucoot_cost)
}

#' Partial Gromov-Wasserstein (square loss)
#'
#' POT-compatible exact partial GW solver using conditional gradient with
#' partial-OT direction steps.
#'
#' @inheritParams gromov_wasserstein
#' @param m Amount of mass to transport. Defaults to `min(sum(p), sum(q))`.
#' @param nb_dummies Number of dummy nodes used in the partial OT linearized step.
#' @param thres Unused POT compatibility parameter.
#' @param numItermax Maximum CG iterations.
#' @param warn Ignored; retained for POT signature compatibility.
#' @param lp_solver LP backend for the linearized partial OT step
#'   (`"cpp_transport"` default, `"lp_transport"`, `"lp_matrix"`).
#' @param lp_scale Integer scaling factor for LP marginal discretization.
#' @param tol Relative stopping tolerance on the partial objective.
#' @param log If `TRUE`, return a list with diagnostics; otherwise the plan.
#' @param verbose If `TRUE`, print CG diagnostics.
#' @return If `log = FALSE`, returns a coupling matrix. If `log = TRUE`, returns
#'   a list with `plan`, `partial_gw_dist`, `iterations`, `error`, and `loss_trace`.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
partial_gromov_wasserstein <- function(
    C1,
    C2,
    p = NULL,
    q = NULL,
    m = NULL,
    loss_fun = "square_loss",
    nb_dummies = 1L,
    G0 = NULL,
    thres = 1,
    numItermax = 10000L,
    tol = 1e-8,
    symmetric = NULL,
    warn = TRUE,
    log = FALSE,
    verbose = FALSE,
    lp_solver = c("cpp_transport", "lp_matrix", "lp_transport"),
    lp_scale = 1e6,
    ...) {
  .check_square_loss(loss_fun)
  .assert_matrix(C1, "C1")
  .assert_matrix(C2, "C2")
  ns <- nrow(C1)
  nt <- nrow(C2)
  if (ncol(C1) != ns || ncol(C2) != nt) {
    stop("`C1` and `C2` must be square.", call. = FALSE)
  }

  if (is.null(p)) p <- rep(1 / ns, ns)
  if (is.null(q)) q <- rep(1 / nt, nt)
  p <- .assert_prob(p, ns, "p")
  q <- .assert_prob(q, nt, "q")
  mass_defaulted <- is.null(m)
  m <- .validate_mass(m, p, q)
  numItermax <- .validate_count(numItermax, "numItermax")
  nb_dummies <- .validate_count(nb_dummies, "nb_dummies")

  symmetric <- .resolve_symmetric(symmetric, C1, C2)

  G0 <- .validate_partial_init(G0, p, q, m, ns, nt)

  out <- .partial_fgw_exact_dispatch(
    M = .empty_feature_cost(),
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    m = m,
    alpha = 1,
    G0 = G0,
    max_iter = numItermax,
    tol = tol,
    symmetric = symmetric,
    lp_solver = match.arg(lp_solver),
    lp_scale = lp_scale,
    nb_dummies = nb_dummies,
    verbose = verbose
  )

  if (!isTRUE(log)) {
    return(out$plan)
  }

  objective_recomputed <- ot_gw_square(C1, C2, out$plan, symmetric = symmetric)
  .partial_log_with_nested_status(
    out, "partial_gw_dist", tol, numItermax,
    p, q, m, 1e-8, objective_recomputed,
    list(lin_loss = 0, quad_loss = objective_recomputed),
    mass_defaulted = mass_defaulted
  )
}

#' Partial Gromov-Wasserstein Objective Value
#'
#' @inheritParams partial_gromov_wasserstein
#' @return Partial GW value.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
partial_gromov_wasserstein2 <- function(...) {
  out <- partial_gromov_wasserstein(..., log = TRUE)
  out$partial_gw_dist
}

#' Partial Fused Gromov-Wasserstein (square loss)
#'
#' POT-compatible exact partial FGW solver using conditional gradient.
#'
#' @inheritParams partial_gromov_wasserstein
#' @param M Cross-domain feature cost matrix.
#' @param alpha FGW tradeoff in `[0, 1]`.
#' @return If `log = FALSE`, returns a coupling matrix. If `log = TRUE`, returns
#'   a list with `plan`, `partial_fgw_dist`, `lin_loss`, `quad_loss`,
#'   `iterations`, `error`, and `loss_trace`.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
partial_fused_gromov_wasserstein <- function(
    M,
    C1,
    C2,
    p = NULL,
    q = NULL,
    m = NULL,
    loss_fun = "square_loss",
    alpha = 0.5,
    nb_dummies = 1L,
    G0 = NULL,
    thres = 1,
    numItermax = 10000L,
    tol = 1e-8,
    symmetric = NULL,
    warn = TRUE,
    log = FALSE,
    verbose = FALSE,
    lp_solver = c("cpp_transport", "lp_matrix", "lp_transport"),
    lp_scale = 1e6,
    ...) {
  .check_square_loss(loss_fun)
  .assert_matrix(M, "M")
  .assert_matrix(C1, "C1")
  .assert_matrix(C2, "C2")

  ns <- nrow(C1)
  nt <- nrow(C2)
  if (ncol(C1) != ns || ncol(C2) != nt) {
    stop("`C1` and `C2` must be square.", call. = FALSE)
  }
  if (nrow(M) != ns || ncol(M) != nt) {
    stop("`M` must have shape nrow(C1) x nrow(C2).", call. = FALSE)
  }
  if (!is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be in [0, 1].", call. = FALSE)
  }

  if (is.null(p)) p <- rep(1 / ns, ns)
  if (is.null(q)) q <- rep(1 / nt, nt)
  p <- .assert_prob(p, ns, "p")
  q <- .assert_prob(q, nt, "q")
  mass_defaulted <- is.null(m)
  m <- .validate_mass(m, p, q)
  numItermax <- .validate_count(numItermax, "numItermax")
  nb_dummies <- .validate_count(nb_dummies, "nb_dummies")

  symmetric <- .resolve_symmetric(symmetric, C1, C2)

  G0 <- .validate_partial_init(G0, p, q, m, ns, nt)

  out <- .partial_fgw_exact_dispatch(
    M = if (alpha >= 1) .empty_feature_cost() else M,
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    m = m,
    alpha = alpha,
    G0 = G0,
    max_iter = numItermax,
    tol = tol,
    symmetric = symmetric,
    lp_solver = match.arg(lp_solver),
    lp_scale = lp_scale,
    nb_dummies = nb_dummies,
    verbose = verbose
  )

  if (!isTRUE(log)) {
    return(out$plan)
  }

  lin_recomputed <- (1 - alpha) * ot_linear_cost(M, out$plan)
  quad_recomputed <- alpha * ot_gw_square(C1, C2, out$plan, symmetric = symmetric)
  .partial_log_with_nested_status(
    out, "partial_fgw_dist", tol, numItermax,
    p, q, m, 1e-8, lin_recomputed + quad_recomputed,
    list(lin_loss = lin_recomputed, quad_loss = quad_recomputed),
    mass_defaulted = mass_defaulted
  )
}

#' Partial Fused Gromov-Wasserstein Objective Value
#'
#' @inheritParams partial_fused_gromov_wasserstein
#' @return Partial FGW value.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
partial_fused_gromov_wasserstein2 <- function(...) {
  out <- partial_fused_gromov_wasserstein(..., log = TRUE)
  out$partial_fgw_dist
}

#' Penalized Variable-Mass Partial Fused Gromov-Wasserstein
#'
#' Optimizes the square-loss FGW objective over nonnegative subcouplings while
#' choosing the transported mass. Unmatched source and target mass each incur
#' `discard_penalty` per unit. The outer problem is nonconvex: a converged
#' result certifies subcoupling feasibility, every exact linear direction
#' solve, objective and line-search consistency, and a Frank-Wolfe stationarity
#' gap, but not global optimality.
#'
#' The unrooted objective is
#'
#' ```text
#' (1 - alpha) * sum(M * G) + alpha * GW_square(C1, C2, G)
#' + discard_penalty * (sum(p) + sum(q) - 2 * sum(G)).
#' ```
#'
#' This differs from [partial_fused_gromov_wasserstein()], which fixes
#' `sum(G)`, and from [fugw_kl()], which uses generalized-KL marginal
#' relaxation.
#'
#' @param M Finite nonnegative source-by-target feature cost.
#' @param C1,C2 Finite square within-domain structure matrices.
#' @param discard_penalty Finite nonnegative cost per unmatched unit on each
#'   marginal. Increasing it weakly favors more transported mass in each exact
#'   linearized subproblem.
#' @param p,q Finite nonnegative source and target measures. Defaults are
#'   uniform probability measures.
#' @param alpha Feature/structure tradeoff in `[0, 1]`.
#' @param G0 Optional feasible nonnegative subcoupling. Its total mass is free;
#'   row and column sums must not exceed `p` and `q`.
#' @param numItermax Maximum Frank-Wolfe updates.
#' @param tol Relative Frank-Wolfe-gap tolerance.
#' @param inner_max_iter Maximum iterations for each certified penalized linear
#'   OT direction solve.
#' @param inner_tol Tolerance for each direction solve.
#' @param symmetric Whether to use the symmetric square-loss formula. By
#'   default it is inferred from `C1` and `C2`.
#' @param trace If `TRUE`, retain every accepted feasible iterate and its
#'   diagnostics in `debug_trace`.
#' @param verbose If `TRUE`, print iteration diagnostics.
#' @return An `rfugw_result` containing the plan, unrooted objective terms,
#'   transported/discarded masses, feasibility and objective certificates,
#'   Frank-Wolfe gap, line-search evidence, and nested exact-solver status.
#' @examples
#' C1 <- matrix(c(0, 1, 1, 0), 2, 2)
#' C2 <- matrix(c(0, 2, 2, 0), 2, 2)
#' M <- matrix(c(0.1, 2, 2, 0.2), 2, 2)
#' out <- penalized_partial_fused_gromov_wasserstein(
#'   M, C1, C2, discard_penalty = 1
#' )
#' out$transported_mass
#' out$frank_wolfe_gap
#' @export
penalized_partial_fused_gromov_wasserstein <- function(
    M,
    C1,
    C2,
    discard_penalty,
    p = NULL,
    q = NULL,
    alpha = 0.5,
    G0 = NULL,
    numItermax = 500L,
    tol = 1e-8,
    inner_max_iter = 20000L,
    inner_tol = 1e-12,
    symmetric = NULL,
    trace = FALSE,
    verbose = FALSE) {
  M <- .validate_finite_matrix(M, "M")
  C1 <- .validate_finite_matrix(C1, "C1", square = TRUE)
  C2 <- .validate_finite_matrix(C2, "C2", square = TRUE)
  ns <- nrow(C1)
  nt <- nrow(C2)
  if (nrow(M) != ns || ncol(M) != nt) {
    stop("`M` must have shape nrow(C1) x nrow(C2).", call. = FALSE)
  }
  if (any(M < 0)) {
    stop("`M` must be nonnegative.", call. = FALSE)
  }
  alpha <- .validate_alpha(alpha)
  if (!is.numeric(discard_penalty) || length(discard_penalty) != 1L ||
      !is.finite(discard_penalty) || discard_penalty < 0) {
    stop("`discard_penalty` must be one finite nonnegative number.", call. = FALSE)
  }
  numItermax <- .validate_count(numItermax, "numItermax")
  inner_max_iter <- .validate_count(inner_max_iter, "inner_max_iter")
  tol <- .validate_positive_scalar(tol, "tol")
  inner_tol <- .validate_positive_scalar(inner_tol, "inner_tol")
  if (!is.logical(trace) || length(trace) != 1L || is.na(trace)) {
    stop("`trace` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  p <- .validate_finite_measure(p, ns, "p", rep(1 / ns, ns))
  q <- .validate_finite_measure(q, nt, "q", rep(1 / nt, nt))
  symmetric <- .resolve_symmetric(symmetric, C1, C2)
  G <- .validate_variable_partial_init(G0, p, q, ns, nt)
  initialization <- if (is.null(G0)) "zero_subcoupling" else "user_subcoupling"

  terms <- .penalized_partial_fgw_terms(
    M, C1, C2, G, p, q, alpha, discard_penalty, symmetric
  )
  loss_trace <- terms$objective
  mass_trace <- terms$transported_mass
  debug_trace <- if (trace) {
    list(list(
      iteration = 0L,
      plan = G,
      objective = terms$objective,
      transported_mass = terms$transported_mass,
      feasibility = .variable_partial_feasibility(G, p, q)
    ))
  } else {
    NULL
  }
  direction_trace <- list()
  updates <- 0L
  frank_wolfe_gap <- Inf
  raw_frank_wolfe_gap <- Inf
  gap_tolerance <- max(tol, tol * (1 + abs(terms$objective)))
  max_line_search_residual <- 0
  max_line_search_tolerance <- 0
  line_search_consistent <- TRUE
  inner_converged <- TRUE
  stationarity_consistent <- TRUE
  stopped_by_gap <- FALSE

  evaluate_direction <- function(iteration, terms, G) {
    linear_cost <- (1 - alpha) * M + alpha * terms$gw_gradient
    inner <- .solve_penalized_partial_fgw_direction(
      linear_cost, p, q, discard_penalty, inner_max_iter, inner_tol
    )
    raw_gap <- if (isTRUE(inner$converged)) {
      sum(terms$gradient * (G - inner$plan))
    } else {
      Inf
    }
    list(
      inner = inner,
      raw_gap = raw_gap,
      summary = list(
        iteration = as.integer(iteration),
        status = inner$status,
        converged = isTRUE(inner$converged),
        residual = inner$residual,
        iterations = inner$iterations,
        transported_mass = inner$transported_mass,
        linear_cost_shift = inner$linear_cost_shift,
        adjusted_discard_penalty = inner$adjusted_discard_penalty,
        raw_frank_wolfe_gap = raw_gap
      )
    )
  }

  for (k in seq_len(numItermax)) {
    direction_eval <- evaluate_direction(k, terms, G)
    direction_trace[[length(direction_trace) + 1L]] <- direction_eval$summary
    if (!isTRUE(direction_eval$inner$converged)) {
      inner_converged <- FALSE
      break
    }

    raw_frank_wolfe_gap <- direction_eval$raw_gap
    gap_tolerance <- max(tol, tol * (1 + abs(terms$objective)))
    stationarity_consistent <- is.finite(raw_frank_wolfe_gap) &&
      raw_frank_wolfe_gap >= -gap_tolerance
    frank_wolfe_gap <- max(0, raw_frank_wolfe_gap)
    if (stationarity_consistent && frank_wolfe_gap <= gap_tolerance) {
      stopped_by_gap <- TRUE
      break
    }
    if (!stationarity_consistent) {
      break
    }

    line <- .penalized_partial_fgw_line_search(
      M, C1, C2, G, direction_eval$inner$plan, alpha,
      discard_penalty, symmetric, current_terms = terms
    )
    next_G <- G + line$step * line$delta
    next_terms <- .penalized_partial_fgw_terms(
      M, C1, C2, next_G, p, q, alpha, discard_penalty, symmetric
    )
    line_residual <- abs(next_terms$objective - line$predicted_objective)
    line_tolerance <- max(
      1e-10,
      20 * tol * (1 + abs(next_terms$objective) +
        abs(line$predicted_objective))
    )
    max_line_search_residual <- max(max_line_search_residual, line_residual)
    max_line_search_tolerance <- max(max_line_search_tolerance, line_tolerance)
    if (!is.finite(line_residual) || line_residual > line_tolerance) {
      line_search_consistent <- FALSE
      break
    }

    G <- next_G
    terms <- next_terms
    updates <- k
    loss_trace <- c(loss_trace, terms$objective)
    mass_trace <- c(mass_trace, terms$transported_mass)
    if (trace) {
      debug_trace[[length(debug_trace) + 1L]] <- list(
        iteration = as.integer(k),
        plan = G,
        objective = terms$objective,
        transported_mass = terms$transported_mass,
        step = line$step,
        line_search_quadratic = line$quadratic,
        line_search_linear = line$linear,
        line_search_residual = line_residual,
        feasibility = .variable_partial_feasibility(G, p, q)
      )
    }
    if (isTRUE(verbose) && (k == 1L || k %% 25L == 0L)) {
      cat(sprintf(
        "iter=%d objective=%.8e mass=%.8e fw_gap=%.3e step=%.3e\n",
        k, terms$objective, terms$transported_mass,
        frank_wolfe_gap, line$step
      ))
    }
  }

  if (inner_converged && line_search_consistent &&
      stationarity_consistent && !stopped_by_gap) {
    direction_eval <- evaluate_direction(updates + 1L, terms, G)
    direction_trace[[length(direction_trace) + 1L]] <- direction_eval$summary
    inner_converged <- isTRUE(direction_eval$inner$converged)
    if (inner_converged) {
      raw_frank_wolfe_gap <- direction_eval$raw_gap
      gap_tolerance <- max(tol, tol * (1 + abs(terms$objective)))
      stationarity_consistent <- is.finite(raw_frank_wolfe_gap) &&
        raw_frank_wolfe_gap >= -gap_tolerance
      frank_wolfe_gap <- max(0, raw_frank_wolfe_gap)
      stopped_by_gap <- stationarity_consistent &&
        frank_wolfe_gap <= gap_tolerance
    } else {
      frank_wolfe_gap <- Inf
    }
  }

  feasibility <- .variable_partial_feasibility(G, p, q)
  feature_recomputed <- (1 - alpha) * ot_linear_cost(M, G)
  structure_recomputed <- alpha * ot_gw_square(
    C1, C2, G, symmetric = symmetric
  )
  penalty_recomputed <- discard_penalty *
    (sum(p) + sum(q) - 2 * sum(G))
  objective_recomputed <- feature_recomputed + structure_recomputed +
    penalty_recomputed
  objective_residual <- abs(terms$objective - objective_recomputed)
  objective_tolerance <- max(
    1e-10,
    20 * tol * (1 + abs(terms$objective) + abs(objective_recomputed))
  )
  objective_consistent <- is.finite(objective_residual) &&
    objective_residual <= objective_tolerance
  components_consistent <- all(is.finite(c(
    feature_recomputed, structure_recomputed, penalty_recomputed
  ))) && abs(
    feature_recomputed + structure_recomputed + penalty_recomputed -
      objective_recomputed
  ) <= objective_tolerance
  trace_feasible <- !trace || all(vapply(
    debug_trace, function(x) isTRUE(x$feasibility$feasible), logical(1)
  ))
  certified <- isTRUE(inner_converged) && isTRUE(line_search_consistent) &&
    isTRUE(stationarity_consistent) && isTRUE(stopped_by_gap) &&
    isTRUE(feasibility$feasible) && isTRUE(objective_consistent) &&
    isTRUE(components_consistent) && isTRUE(trace_feasible)
  status <- if (!inner_converged) {
    "inner_failure"
  } else if (!line_search_consistent) {
    "objective_mismatch"
  } else if (!stationarity_consistent) {
    "stationarity_failure"
  } else if (!feasibility$feasible || !trace_feasible) {
    "infeasible"
  } else if (!objective_consistent || !components_consistent) {
    "objective_mismatch"
  } else if (certified) {
    "converged"
  } else {
    "max_iter"
  }
  termination_reason <- switch(
    status,
    converged = "frank_wolfe_gap",
    max_iter = "maximum_iterations",
    inner_failure = "penalized_linear_direction_failure",
    stationarity_failure = "negative_frank_wolfe_gap",
    infeasible = "subcoupling_feasibility_failure",
    objective_mismatch = if (!line_search_consistent) {
      "line_search_polynomial_mismatch"
    } else {
      "objective_recomputation_mismatch"
    }
  )
  inner_iterations <- sum(vapply(
    direction_trace, function(x) as.numeric(x$iterations), numeric(1)
  ), na.rm = TRUE)
  inner_residuals <- vapply(
    direction_trace, function(x) as.numeric(x$residual), numeric(1)
  )

  out <- list(
    plan = G,
    penalized_partial_objective = terms$objective,
    partial_fgw_dist = terms$objective,
    feature_term = terms$feature_term,
    structure_term = terms$structure_term,
    discard_penalty_term = terms$penalty_term,
    feature_cost_unweighted = terms$feature_raw,
    structure_cost_unweighted = terms$structure_raw,
    transported_mass = terms$transported_mass,
    mass = terms$transported_mass,
    discarded_source_mass = terms$discarded_source_mass,
    discarded_target_mass = terms$discarded_target_mass,
    original_source_mass = sum(p),
    original_target_mass = sum(q),
    source_measure_original = p,
    target_measure_original = q,
    measure_normalization = "finite_measure_preserved",
    discard_penalty = discard_penalty,
    penalty_direction =
      "larger_discard_penalty_weakly_favors_more_transport_in_linearized_subproblems",
    alpha = alpha,
    symmetric = symmetric,
    iterations = as.integer(updates),
    max_iter = as.integer(numItermax),
    residual = frank_wolfe_gap,
    error = frank_wolfe_gap,
    frank_wolfe_gap = frank_wolfe_gap,
    raw_frank_wolfe_gap = raw_frank_wolfe_gap,
    frank_wolfe_gap_tolerance = gap_tolerance,
    stationarity_consistent = stationarity_consistent,
    loss_trace = loss_trace,
    mass_trace = mass_trace,
    line_search_residual = max_line_search_residual,
    line_search_tolerance = max_line_search_tolerance,
    line_search_consistent = line_search_consistent,
    feasibility = "partial_variable_mass",
    feasibility_residual = feasibility$residual,
    feasibility_tolerance = feasibility$tolerance,
    feasible = feasibility$feasible,
    row_residual = feasibility$row_residual,
    col_residual = feasibility$col_residual,
    nonnegative_residual = feasibility$nonnegative_residual,
    mass_residual = feasibility$mass_residual,
    mass_target = NA_real_,
    mass_certified = feasibility$feasible,
    mass_certification = "optimized_subcoupling_mass",
    objective_recomputed = objective_recomputed,
    objective_residual = objective_residual,
    objective_tolerance = objective_tolerance,
    objective_consistent = objective_consistent,
    objective_components_consistent = components_consistent,
    feature_term_recomputed = feature_recomputed,
    structure_term_recomputed = structure_recomputed,
    discard_penalty_term_recomputed = penalty_recomputed,
    inner_converged = inner_converged,
    inner_status = if (inner_converged) "converged" else "failed",
    inner_iterations = as.integer(inner_iterations),
    inner_residual = if (length(inner_residuals)) tail(inner_residuals, 1L) else Inf,
    max_inner_residual = if (length(inner_residuals)) {
      max(inner_residuals, na.rm = TRUE)
    } else {
      Inf
    },
    direction_trace = direction_trace,
    debug_trace = debug_trace,
    trace_feasible = trace_feasible,
    initialization = initialization,
    formulation = "penalized_partial_fgw_square",
    backend = "frank_wolfe_certified_penalized_linear_ot",
    regularization = 0,
    status = status,
    converged = identical(status, "converged"),
    termination_reason = termination_reason,
    warning_payload = if (identical(status, "converged")) NULL else list(
      code = status,
      message = sprintf("Penalized partial FGW is not certified: %s.", status)
    )
  )
  out$runtime_provenance <- .solver_runtime_provenance(out)
  class(out) <- unique(c("rfugw_result", class(out)))
  out
}

#' Entropic Partial Gromov-Wasserstein (square loss)
#'
#' POT-compatible entropic partial GW solver via alternating gradient and
#' entropic partial OT projections.
#'
#' @inheritParams partial_gromov_wasserstein
#' @param reg Entropic regularization parameter (>0).
#' @param check_every Outer stopping check interval.
#' @param inner_max_iter Maximum iterations for inner entropic partial OT solve.
#' @param inner_tol Inner stopping tolerance for partial OT solve.
#' @param method Certified public partial-Sinkhorn backend: `"auto"`, bounded
#'   `"scaling"`, or genuine `"log"` Dykstra.
#' @return If `log = FALSE`, returns a coupling matrix. If `log = TRUE`, returns
#'   a list with `plan`, `partial_gw_dist`, `iterations`, `error`, and `err_trace`.
#' @export
entropic_partial_gromov_wasserstein <- function(
    C1,
    C2,
    p = NULL,
    q = NULL,
    reg = 1.0,
    m = NULL,
    loss_fun = "square_loss",
    G0 = NULL,
    numItermax = 1000L,
    tol = 1e-7,
    symmetric = NULL,
    log = FALSE,
    verbose = FALSE,
    check_every = 2L,
    inner_max_iter = 300L,
    inner_tol = 1e-12,
    method = c("auto", "scaling", "log")) {
  method <- match.arg(method)
  .check_square_loss(loss_fun)
  .assert_matrix(C1, "C1")
  .assert_matrix(C2, "C2")
  ns <- nrow(C1)
  nt <- nrow(C2)
  if (ncol(C1) != ns || ncol(C2) != nt) {
    stop("`C1` and `C2` must be square.", call. = FALSE)
  }
  if (!is.finite(reg) || reg <= 0) {
    stop("`reg` must be positive.", call. = FALSE)
  }

  if (is.null(p)) p <- rep(1 / ns, ns)
  if (is.null(q)) q <- rep(1 / nt, nt)
  p <- .assert_prob(p, ns, "p")
  q <- .assert_prob(q, nt, "q")
  mass_defaulted <- is.null(m)
  m <- .validate_mass(m, p, q)
  numItermax <- .validate_count(numItermax, "numItermax")
  inner_max_iter <- .validate_count(inner_max_iter, "inner_max_iter")
  check_every <- .validate_count(check_every, "check_every")

  symmetric <- .resolve_symmetric(symmetric, C1, C2)
  G0 <- .validate_partial_init(G0, p, q, m, ns, nt)

  out <- .partial_fgw_entropic_dispatch(
    M = .empty_feature_cost(),
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    m = m,
    reg = reg,
    alpha = 1,
    G0 = G0,
    max_iter = numItermax,
    tol = tol,
    symmetric = symmetric,
    inner_max_iter = inner_max_iter,
    inner_tol = inner_tol,
    verbose = verbose,
    check_every = check_every,
    method = method
  )

  if (!isTRUE(log)) {
    plan <- out$plan
    attr(plan, "sinkhorn_dispatch") <- out$sinkhorn_dispatch
    return(plan)
  }

  objective_recomputed <- ot_gw_square(C1, C2, out$plan, symmetric = symmetric)
  ans <- .partial_log_with_nested_status(
    out, "partial_gw_dist", tol, numItermax,
    p, q, m, inner_tol, objective_recomputed,
    list(lin_loss = 0, quad_loss = objective_recomputed),
    mass_defaulted = mass_defaulted
  )
  .attach_partial_sinkhorn_dispatch(ans, out$sinkhorn_dispatch)
}

#' Entropic Partial Gromov-Wasserstein Objective Value
#'
#' @inheritParams entropic_partial_gromov_wasserstein
#' @return Entropic partial GW value.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
entropic_partial_gromov_wasserstein2 <- function(...) {
  out <- entropic_partial_gromov_wasserstein(..., log = TRUE)
  out$partial_gw_dist
}

#' Entropic Partial Fused Gromov-Wasserstein (square loss)
#'
#' POT-compatible entropic partial FGW solver.
#'
#' @inheritParams entropic_partial_gromov_wasserstein
#' @param M Cross-domain feature cost matrix.
#' @param alpha FGW tradeoff in `[0, 1]`.
#' @return If `log = FALSE`, returns a coupling matrix. If `log = TRUE`, returns
#'   a list with `plan`, `partial_fgw_dist`, `lin_loss`, `quad_loss`,
#'   `iterations`, `error`, and `err_trace`.
#' @export
entropic_partial_fused_gromov_wasserstein <- function(
    M,
    C1,
    C2,
    p = NULL,
    q = NULL,
    reg = 1.0,
    m = NULL,
    loss_fun = "square_loss",
    alpha = 0.5,
    G0 = NULL,
    numItermax = 1000L,
    tol = 1e-7,
    symmetric = NULL,
    log = FALSE,
    verbose = FALSE,
    check_every = 2L,
    inner_max_iter = 300L,
    inner_tol = 1e-12,
    method = c("auto", "scaling", "log")) {
  method <- match.arg(method)
  .check_square_loss(loss_fun)
  .assert_matrix(M, "M")
  .assert_matrix(C1, "C1")
  .assert_matrix(C2, "C2")
  ns <- nrow(C1)
  nt <- nrow(C2)
  if (ncol(C1) != ns || ncol(C2) != nt) {
    stop("`C1` and `C2` must be square.", call. = FALSE)
  }
  if (nrow(M) != ns || ncol(M) != nt) {
    stop("`M` must have shape nrow(C1) x nrow(C2).", call. = FALSE)
  }
  if (!is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be in [0, 1].", call. = FALSE)
  }
  if (!is.finite(reg) || reg <= 0) {
    stop("`reg` must be positive.", call. = FALSE)
  }

  if (is.null(p)) p <- rep(1 / ns, ns)
  if (is.null(q)) q <- rep(1 / nt, nt)
  p <- .assert_prob(p, ns, "p")
  q <- .assert_prob(q, nt, "q")
  mass_defaulted <- is.null(m)
  m <- .validate_mass(m, p, q)
  numItermax <- .validate_count(numItermax, "numItermax")
  inner_max_iter <- .validate_count(inner_max_iter, "inner_max_iter")
  check_every <- .validate_count(check_every, "check_every")

  symmetric <- .resolve_symmetric(symmetric, C1, C2)
  G0 <- .validate_partial_init(G0, p, q, m, ns, nt)

  out <- .partial_fgw_entropic_dispatch(
    M = if (alpha >= 1) .empty_feature_cost() else M,
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    m = m,
    reg = reg,
    alpha = alpha,
    G0 = G0,
    max_iter = numItermax,
    tol = tol,
    symmetric = symmetric,
    inner_max_iter = inner_max_iter,
    inner_tol = inner_tol,
    verbose = verbose,
    check_every = check_every,
    method = method
  )

  if (!isTRUE(log)) {
    plan <- out$plan
    attr(plan, "sinkhorn_dispatch") <- out$sinkhorn_dispatch
    return(plan)
  }

  lin_recomputed <- (1 - alpha) * ot_linear_cost(M, out$plan)
  quad_recomputed <- alpha * ot_gw_square(C1, C2, out$plan, symmetric = symmetric)
  ans <- .partial_log_with_nested_status(
    out, "partial_fgw_dist", tol, numItermax,
    p, q, m, inner_tol, lin_recomputed + quad_recomputed,
    list(lin_loss = lin_recomputed, quad_loss = quad_recomputed),
    mass_defaulted = mass_defaulted
  )
  .attach_partial_sinkhorn_dispatch(ans, out$sinkhorn_dispatch)
}

#' Entropic Partial Fused Gromov-Wasserstein Objective Value
#'
#' @inheritParams entropic_partial_fused_gromov_wasserstein
#' @return Entropic partial FGW value.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
entropic_partial_fused_gromov_wasserstein2 <- function(...) {
  out <- entropic_partial_fused_gromov_wasserstein(..., log = TRUE)
  out$partial_fgw_dist
}

#' Semi-Relaxed Gromov-Wasserstein (non-entropic, square loss)
#'
#' POT-compatible non-entropic semirelaxed GW solver via conditional gradient.
#'
#' @inheritParams entropic_semirelaxed_gromov_wasserstein
#' @param tol_rel Relative stopping tolerance.
#' @param tol_abs Absolute stopping tolerance.
#' @param log If `TRUE`, include `loss_trace`.
#' @param random_state Ignored; retained for POT signature compatibility.
#' @return A list with `plan`, `q`, `srgw_dist`, `iterations`, `error`, and
#'   `abs_error`.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
semirelaxed_gromov_wasserstein <- function(
    C1,
    C2,
    p = NULL,
    loss_fun = "square_loss",
    symmetric = NULL,
    log = FALSE,
    G0 = NULL,
    max_iter = 10000L,
    tol_rel = 1e-9,
    tol_abs = 1e-9,
    random_state = 0,
    verbose = FALSE,
    ...) {
  .check_square_loss(loss_fun)
  max_iter <- .validate_count(max_iter, "max_iter")
  .assert_matrix(C1, "C1")
  .assert_matrix(C2, "C2")
  ns <- nrow(C1)
  nt <- nrow(C2)
  if (ncol(C1) != ns || ncol(C2) != nt) {
    stop("`C1` and `C2` must be square.", call. = FALSE)
  }

  if (is.null(p)) p <- rep(1 / ns, ns)
  p <- .assert_prob(p, ns, "p")

  if (is.null(symmetric)) {
    symmetric <- .is_symmetric_cost(C1) && .is_symmetric_cost(C2)
  } else {
    symmetric <- isTRUE(symmetric)
  }

  out <- .semirelaxed_fgw_exact_dispatch(
    M = .empty_feature_cost(),
    C1 = C1,
    C2 = C2,
    p = p,
    alpha = 1,
    symmetric = symmetric,
    G0 = G0,
    max_iter = max_iter,
    tol_rel = tol_rel,
    tol_abs = tol_abs,
    verbose = verbose
  )

  res <- list(
    plan = out$plan,
    q = out$q,
    srgw_dist = out$srgw_dist,
    iterations = out$iterations,
    error = out$error,
    abs_error = out$abs_error,
    symmetric = out$symmetric
  )

  if (isTRUE(log)) {
    res$loss_trace <- out$loss_trace
  }
  converged <- (is.finite(out$error) && out$error <= tol_rel) ||
    (is.finite(out$abs_error) && out$abs_error <= tol_abs)
  ans <- .attach_solver_diagnostics(
    res,
    residual = out$abs_error,
    converged = converged,
    iterations = out$iterations,
    max_iter = max_iter,
    p = p,
    plan = out$plan,
    feasibility = "semirelaxed",
    feasibility_tol = 1e-8,
    objective_recomputed = ot_gw_square(
      C1, C2, out$plan, symmetric = symmetric
    )
  )
  ans$termination_reason <- .termination_reason_from_result(ans, max_iter)
  ans
}

#' Semi-Relaxed Gromov-Wasserstein Objective Value
#'
#' @inheritParams semirelaxed_gromov_wasserstein
#' @return Semirelaxed GW value.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
semirelaxed_gromov_wasserstein2 <- function(...) {
  out <- semirelaxed_gromov_wasserstein(...)
  out$srgw_dist
}

#' Semi-Relaxed Fused Gromov-Wasserstein (non-entropic, square loss)
#'
#' POT-compatible non-entropic semirelaxed FGW solver via conditional gradient.
#'
#' @inheritParams entropic_semirelaxed_fused_gromov_wasserstein
#' @param tol_rel Relative stopping tolerance.
#' @param tol_abs Absolute stopping tolerance.
#' @param log If `TRUE`, include `loss_trace`.
#' @param random_state Ignored; retained for POT signature compatibility.
#' @return A list with `plan`, `q`, `srfgw_dist`, `lin_loss`, `quad_loss`,
#'   `iterations`, `error`, and `abs_error`.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
semirelaxed_fused_gromov_wasserstein <- function(
    M,
    C1,
    C2,
    p = NULL,
    loss_fun = "square_loss",
    symmetric = NULL,
    alpha = 0.5,
    G0 = NULL,
    log = FALSE,
    max_iter = 10000L,
    tol_rel = 1e-9,
    tol_abs = 1e-9,
    random_state = 0,
    verbose = FALSE,
    ...) {
  .check_square_loss(loss_fun)
  max_iter <- .validate_count(max_iter, "max_iter")
  .assert_matrix(M, "M")
  .assert_matrix(C1, "C1")
  .assert_matrix(C2, "C2")
  ns <- nrow(C1)
  nt <- nrow(C2)
  if (ncol(C1) != ns || ncol(C2) != nt) {
    stop("`C1` and `C2` must be square.", call. = FALSE)
  }
  if (nrow(M) != ns || ncol(M) != nt) {
    stop("`M` must have shape nrow(C1) x nrow(C2).", call. = FALSE)
  }
  if (!is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be in [0, 1].", call. = FALSE)
  }

  if (is.null(p)) p <- rep(1 / ns, ns)
  p <- .assert_prob(p, ns, "p")

  if (is.null(symmetric)) {
    symmetric <- .is_symmetric_cost(C1) && .is_symmetric_cost(C2)
  } else {
    symmetric <- isTRUE(symmetric)
  }

  out <- .semirelaxed_fgw_exact_dispatch(
    M = if (alpha >= 1) .empty_feature_cost() else M,
    C1 = C1,
    C2 = C2,
    p = p,
    alpha = alpha,
    symmetric = symmetric,
    G0 = G0,
    max_iter = max_iter,
    tol_rel = tol_rel,
    tol_abs = tol_abs,
    verbose = verbose
  )

  res <- list(
    plan = out$plan,
    q = out$q,
    srfgw_dist = out$srfgw_dist,
    lin_loss = out$lin_loss,
    quad_loss = out$quad_loss,
    iterations = out$iterations,
    error = out$error,
    abs_error = out$abs_error,
    symmetric = out$symmetric
  )

  if (isTRUE(log)) {
    res$loss_trace <- out$loss_trace
  }
  converged <- (is.finite(out$error) && out$error <= tol_rel) ||
    (is.finite(out$abs_error) && out$abs_error <= tol_abs)
  ans <- .attach_solver_diagnostics(
    res,
    residual = out$abs_error,
    converged = converged,
    iterations = out$iterations,
    max_iter = max_iter,
    p = p,
    plan = out$plan,
    feasibility = "semirelaxed",
    feasibility_tol = 1e-8,
    objective_recomputed = ot_fgw_square(
      M, C1, C2, out$plan, alpha = alpha, symmetric = symmetric
    ),
    objective_components = list(
      lin_loss = (1 - alpha) * ot_linear_cost(M, out$plan),
      quad_loss = alpha * ot_gw_square(C1, C2, out$plan, symmetric = symmetric)
    )
  )
  ans$termination_reason <- .termination_reason_from_result(ans, max_iter)
  ans
}

#' Semi-Relaxed Fused Gromov-Wasserstein Objective Value
#'
#' @inheritParams semirelaxed_fused_gromov_wasserstein
#' @return Semirelaxed FGW value.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
semirelaxed_fused_gromov_wasserstein2 <- function(...) {
  out <- semirelaxed_fused_gromov_wasserstein(...)
  out$srfgw_dist
}

#' Entropic Gromov-Wasserstein Barycenters
#'
#' POT-compatible GW barycenter wrapper using the optimized FGW barycenter core
#' with zero features and `alpha = 1`.
#'
#' @param N Number of barycenter nodes.
#' @param Cs List of structure matrices.
#' @param ps Optional list of source weights.
#' @param p Optional barycenter weights.
#' @param lambdas Optional barycenter sample weights.
#' @param loss_fun Currently only `"square_loss"`.
#' @param epsilon Entropic regularization.
#' @param symmetric Symmetry flag for inner GW solves.
#' @param max_iter Max outer barycenter iterations.
#' @param tol Outer stopping tolerance.
#' @param stop_criterion Currently only `"barycenter"` is supported.
#' @param warmstartT Warm-start inner GW solves.
#' @param verbose Print diagnostics.
#' @param log Return diagnostics/history.
#' @param init_C Optional initial barycenter structure.
#' @param random_state Optional seed for initialization.
#' @return If `log = FALSE`, returns barycenter structure matrix `C`. If
#'   `log = TRUE`, returns a list with `C`, `p`, `couplings`, `history`, and
#'   `objective`.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
entropic_gromov_barycenters <- function(
    N,
    Cs,
    ps = NULL,
    p = NULL,
    lambdas = NULL,
    loss_fun = "square_loss",
    epsilon = 0.1,
    symmetric = TRUE,
    max_iter = 1000L,
    tol = 1e-9,
    stop_criterion = c("barycenter", "loss"),
    warmstartT = FALSE,
    verbose = FALSE,
    log = FALSE,
    init_C = NULL,
    random_state = NULL,
    ...) {
  .check_square_loss(loss_fun)
  stop_criterion <- match.arg(stop_criterion)
  if (!identical(stop_criterion, "barycenter")) {
    stop("`stop_criterion = \"loss\"` is not currently supported.", call. = FALSE)
  }

  Ys <- lapply(Cs, function(C) matrix(0, nrow = nrow(C), ncol = 1L))
  init_X <- matrix(0, nrow = as.integer(N), ncol = 1L)

  out <- fgw_barycenters(
    N = N,
    Ys = Ys,
    Cs = Cs,
    ps = ps,
    p = p,
    lambdas = lambdas,
    alpha = 1,
    fixed_structure = FALSE,
    fixed_features = TRUE,
    loss_fun = loss_fun,
    epsilon = epsilon,
    symmetric = symmetric,
    max_iter = as.integer(max_iter),
    tol = tol,
    warmstartT = warmstartT,
    init_C = init_C,
    init_X = init_X,
    random_state = random_state,
    verbose = verbose,
    ...
  )

  if (!isTRUE(log)) {
    return(out$C)
  }

  list(
    C = out$C,
    p = out$p,
    couplings = out$couplings,
    history = out$history,
    objective = out$objective,
    iterations = out$iterations,
    error = out$error
  )
}

#' Entropic Fused Gromov-Wasserstein Barycenters
#'
#' POT-compatible alias for fixed-support entropic FGW barycenters.
#'
#' @inheritParams fgw_barycenters
#' @param stop_criterion Currently only `"barycenter"` is supported.
#' @param init_Y Optional initial barycenter features.
#' @param log If `TRUE`, return full diagnostics; otherwise `Y` and `C`.
#' @return If `log = FALSE`, returns a list with `Y` and `C`. If `log = TRUE`,
#'   returns full diagnostics.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
entropic_fused_gromov_barycenters <- function(
    N,
    Ys,
    Cs,
    ps = NULL,
    p = NULL,
    lambdas = NULL,
    loss_fun = "square_loss",
    epsilon = 0.1,
    symmetric = TRUE,
    alpha = 0.5,
    max_iter = 1000L,
    tol = 1e-9,
    stop_criterion = c("barycenter", "loss"),
    warmstartT = FALSE,
    verbose = FALSE,
    log = FALSE,
    init_C = NULL,
    init_Y = NULL,
    fixed_structure = FALSE,
    fixed_features = FALSE,
    random_state = NULL,
    ...) {
  .check_square_loss(loss_fun)
  stop_criterion <- match.arg(stop_criterion)
  if (!identical(stop_criterion, "barycenter")) {
    stop("`stop_criterion = \"loss\"` is not currently supported.", call. = FALSE)
  }

  out <- fgw_barycenters(
    N = N,
    Ys = Ys,
    Cs = Cs,
    ps = ps,
    p = p,
    lambdas = lambdas,
    alpha = alpha,
    fixed_structure = fixed_structure,
    fixed_features = fixed_features,
    loss_fun = loss_fun,
    epsilon = epsilon,
    symmetric = symmetric,
    max_iter = as.integer(max_iter),
    tol = tol,
    warmstartT = warmstartT,
    init_C = init_C,
    init_X = init_Y,
    random_state = random_state,
    verbose = verbose,
    ...
  )

  if (!isTRUE(log)) {
    return(list(Y = out$X, C = out$C))
  }

  list(
    Y = out$X,
    C = out$C,
    p = out$p,
    couplings = out$couplings,
    history = out$history,
    objective = out$objective,
    iterations = out$iterations,
    error = out$error
  )
}

#' Gromov-Wasserstein Barycenters
#'
#' Convenience alias to `entropic_gromov_barycenters`.
#'
#' @inheritParams entropic_gromov_barycenters
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
gromov_barycenters <- function(...) {
  entropic_gromov_barycenters(...)
}

#' Fused Gromov-Wasserstein Barycenters
#'
#' Convenience alias to `entropic_fused_gromov_barycenters`.
#'
#' @inheritParams entropic_fused_gromov_barycenters
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
fused_gromov_barycenters <- function(...) {
  entropic_fused_gromov_barycenters(...)
}

#' Sampled Gromov-Wasserstein (square loss)
#'
#' POT-style stochastic GW estimator using sampled gradients and balanced OT
#' projection steps.
#'
#' @param C1 Source structure matrix.
#' @param C2 Target structure matrix.
#' @param p Source weights (default uniform).
#' @param q Target weights (default uniform).
#' @param loss_fun Currently only `"square_loss"` is supported.
#' @param nb_samples_grad Number of sampled gradient points, or length-2 vector
#'   `(n_source_samples, n_target_samples)`. Values below 1 error. Source or
#'   target counts above `ns` / `nt` warn and clamp.
#' @param epsilon Entropic regularization for the OT projection step. If `<= 0`,
#'   exact LP projection is used.
#' @param max_iter Maximum stochastic iterations.
#' @param log If `TRUE`, return objective estimate and diagnostics.
#' @param verbose If `TRUE`, print iterative diagnostics.
#' @param random_state Optional seed.
#' @param sinkhorn_max_iter Sinkhorn iterations when `epsilon > 0`.
#' @param sinkhorn_tol Sinkhorn tolerance when `epsilon > 0`.
#' @param lp_solver LP backend used when `epsilon <= 0`.
#' @param lp_scale Integer scaling for LP marginals.
#' @return If `log = FALSE`, returns coupling matrix `T`. If `log = TRUE`,
#'   returns a list with `plan`, `gw_dist_estimated`, and `iterations`.
#'
#' @section Experimental:
#' Sampled GW is experimental. The certified 0.1 envelope is that a full
#' budget `(ns, nt)` is closer to dense [entropic_gromov_wasserstein()]
#' than a tiny budget such as `(2, 1)`, in square-loss GW and plan
#' Frobenius distance. Intermediate budgets are not certified as
#' monotone. See `inst/bench/sampled-budget-curves.md`.
#' @export
sampled_gromov_wasserstein <- function(
    C1,
    C2,
    p = NULL,
    q = NULL,
    loss_fun = "square_loss",
    nb_samples_grad = 100L,
    epsilon = 1,
    max_iter = 500L,
    log = FALSE,
    verbose = FALSE,
    random_state = NULL,
    sinkhorn_max_iter = 200L,
    sinkhorn_tol = 1e-9,
    lp_solver = c("lp_matrix", "lp_transport"),
    lp_scale = 1e6) {
  .check_square_loss(loss_fun)
  max_iter <- .validate_count(max_iter, "max_iter")
  sinkhorn_max_iter <- .validate_count(sinkhorn_max_iter, "sinkhorn_max_iter")
  .assert_matrix(C1, "C1")
  .assert_matrix(C2, "C2")
  ns <- nrow(C1)
  nt <- nrow(C2)
  if (ncol(C1) != ns || ncol(C2) != nt) {
    stop("`C1` and `C2` must be square.", call. = FALSE)
  }

  if (is.null(p)) p <- rep(1 / ns, ns)
  if (is.null(q)) q <- rep(1 / nt, nt)
  p <- .assert_prob(p, ns, "p")
  q <- .assert_prob(q, nt, "q")

  budget <- .parse_sampled_budget(nb_samples_grad, ns, nt)
  nb_p <- budget$nb_p
  nb_q <- budget$nb_q

  if (!is.null(random_state)) set.seed(as.integer(random_state))

  T <- p %o% q
  symmetric <- .is_symmetric_cost(C1) && .is_symmetric_cost(C2)
  it_last <- 0L

  if (epsilon > 0 && exists("cpp_sampled_gromov_wasserstein_entropic_square", mode = "function")) {
    use_mixed_precision <- .runtime_env_bool("RFUGW_SAMPLED_MIXED", FALSE)
    out_cpp <- cpp_sampled_gromov_wasserstein_entropic_square(
      C1 = C1,
      C2 = C2,
      p = p,
      q = q,
      nb_p = as.integer(nb_p),
      nb_q = as.integer(nb_q),
      epsilon = epsilon,
      max_iter = max_iter,
      sinkhorn_max_iter = sinkhorn_max_iter,
      sinkhorn_tol = sinkhorn_tol,
      symmetric = symmetric,
      init_plan = T,
      verbose = isTRUE(verbose),
      use_mixed_precision = use_mixed_precision
    )
    T <- out_cpp$plan
    it_last <- as.integer(out_cpp$iterations)
  } else {
    continue_small <- 0L
    lp_solver <- match.arg(lp_solver)
    lp_direction <- if (epsilon <= 0) .make_transport_lp_solver(p, q, scale = lp_scale, solver = lp_solver) else NULL

    for (it in seq_len(max_iter)) {
      it_last <- it
      idx0 <- sample.int(ns, size = min(nb_p, ns), prob = p, replace = FALSE)
      Lik <- matrix(0, nrow = ns, ncol = nt)

      for (i in idx0) {
        row_prob <- T[i, ]
        srow <- sum(row_prob)
        if (!is.finite(srow) || srow <= 0) {
          row_prob <- q
        } else {
          row_prob <- row_prob / srow
        }
        kq <- min(nb_q, nt)
        replace_q <- sum(row_prob > 0) < kq
        idx1 <- sample.int(nt, size = kq, prob = row_prob, replace = replace_q)

        C2_sub <- C2[idx1, , drop = FALSE]
        mu2 <- colMeans(C2_sub)
        mu2_sq <- colMeans(C2_sub^2)

        block <- tcrossprod(C1[i, ]^2, rep(1, nt)) +
          tcrossprod(rep(1, ns), mu2_sq) -
          2 * tcrossprod(C1[i, ], mu2)

        if (!symmetric && stats::runif(1) > 0.5) {
          C2_sub_t <- C2[, idx1, drop = FALSE]
          mu2_t <- rowMeans(C2_sub_t)
          mu2_t_sq <- rowMeans(C2_sub_t^2)
          block <- tcrossprod(C1[, i]^2, rep(1, nt)) +
            tcrossprod(rep(1, ns), mu2_t_sq) -
            2 * tcrossprod(C1[, i], mu2_t)
        }

        Lik <- Lik + block
      }

      max_lik <- suppressWarnings(max(Lik, na.rm = TRUE))
      if (is.finite(max_lik) && max_lik > 0) {
        Lik <- Lik / max_lik
      }

      if (epsilon > 0) {
        logT <- log(pmax(T, exp(-200)))
        logT[logT <= -200] <- -Inf
        Lik_eff <- Lik - epsilon * logT
        new_T <- .sinkhorn_balanced(
          a = p,
          b = q,
          M = Lik_eff,
          reg = epsilon,
          max_iter = sinkhorn_max_iter,
          tol = sinkhorn_tol
        )
      } else {
        new_T <- lp_direction(Lik)
      }

      change_T <- mean((T - new_T)^2)
      if (!is.finite(change_T)) {
        break
      }

      if (change_T <= 1e-19) {
        continue_small <- continue_small + 1L
        if (continue_small > 100L) {
          T <- new_T
          break
        }
      } else {
        continue_small <- 0L
      }

      if (isTRUE(verbose) && (it %% 10L == 0L || it == 1L)) {
        cat(sprintf("iter=%d change=%.8e\n", it, change_T))
      }

      T <- new_T
    }
  }

  if (!isTRUE(log)) {
    return(T)
  }

  gw_est <- .gw_square_value(C1, C2, T, p, q, symmetric = symmetric)
  list(
    plan = T,
    gw_dist_estimated = gw_est,
    iterations = as.integer(it_last),
    status = "experimental",
    converged = FALSE,
    termination_reason = "experimental_no_convergence_certificate",
    certification = "experimental_no_convergence_claim"
  )
}

#' Dense GW Plan Followed by Truncated SVD (Experimental)
#'
#' Computes a dense entropic GW plan from sample distances, then compresses
#' that completed plan with truncated SVD. The name is deliberately explicit:
#' this is not an end-to-end low-rank or low-memory GW solver.
#'
#' @param X_s Source samples (`ns x d`).
#' @param X_t Target samples (`nt x d`).
#' @param a Optional source weights.
#' @param b Optional target weights.
#' @param reg Entropic regularization used for the inner GW solve.
#' @param rank Target low-rank factorization rank. Values below 1 error.
#'   Values above `min(ns, nt)` warn and clamp.
#' @param rescale_cost Whether to normalize structure costs to `[0, 1]`.
#' @param stopThr Outer GW stopping tolerance.
#' @param numItermax Outer GW max iterations.
#' @param log If `TRUE`, return diagnostics.
#' @return A list with `Q`, `R`, `g`, `representation`,
#'   `dense_plan_materialized`, and byte counts for the dense structures and
#'   plan. If `log = TRUE`, it also includes `plan`, `value_quad`, and `value`.
#'
#' @section Experimental:
#' The solve materializes two dense square structure costs and a dense
#' rectangular plan before factorization. Peak memory is therefore not bounded
#' by the returned rank. Reconstruction error
#' `||T - T_r||_F / ||T||_F` decreases as rank
#' increases up to `min(ns, nt)`. See
#' `inst/bench/sampled-budget-curves.md`.
#' @export
dense_gromov_wasserstein_plan_svd <- function(
    X_s,
    X_t,
    a = NULL,
    b = NULL,
    reg = 0,
    rank = NULL,
    rescale_cost = TRUE,
    stopThr = 1e-4,
    numItermax = 1000L,
    log = FALSE) {
  numItermax <- .validate_count(numItermax, "numItermax")
  .assert_matrix(X_s, "X_s")
  .assert_matrix(X_t, "X_t")
  ns <- nrow(X_s)
  nt <- nrow(X_t)

  if (is.null(a)) a <- rep(1 / ns, ns)
  if (is.null(b)) b <- rep(1 / nt, nt)
  a <- .assert_prob(a, ns, "a")
  b <- .assert_prob(b, nt, "b")

  C1 <- as.matrix(stats::dist(X_s))
  C2 <- as.matrix(stats::dist(X_t))
  if (isTRUE(rescale_cost)) {
    if (max(C1) > 0) C1 <- C1 / max(C1)
    if (max(C2) > 0) C2 <- C2 / max(C2)
  }

  if (!is.finite(reg) || reg <= 0) {
    reg <- 0.05
  }

  out <- entropic_gromov_wasserstein(
    C1 = C1,
    C2 = C2,
    p = a,
    q = b,
    epsilon = reg,
    max_iter = numItermax,
    tol = stopThr,
    precision = "mixed",
    sinkhorn_max_iter = 500L,
    sinkhorn_tol = 1e-9,
    solver = "PGD"
  )

  T <- out$plan
  r <- .parse_lowrank_rank(rank, ns, nt)

  sv <- svd(T, nu = r, nv = r)
  s <- pmax(sv$d[seq_len(r)], 0)
  S <- diag(sqrt(s), nrow = r, ncol = r)
  Q <- sv$u[, seq_len(r), drop = FALSE] %*% S
  R <- sv$v[, seq_len(r), drop = FALSE] %*% S
  g <- rep(1, r)

  representation <- "posthoc_svd_of_dense_gw_plan"
  dense_structure_bytes <- as.numeric(object.size(C1) + object.size(C2))
  dense_plan_bytes <- as.numeric(object.size(T))
  common <- list(
    Q = Q,
    R = R,
    g = g,
    representation = representation,
    dense_plan_materialized = TRUE,
    dense_structure_bytes = dense_structure_bytes,
    dense_plan_bytes = dense_plan_bytes,
    solve_memory_order = "O(ns^2 + nt^2 + ns*nt)",
    factor_rank = r
  )

  if (!isTRUE(log)) return(common)

  value_quad <- .gw_square_value(C1, C2, T, a, b, symmetric = TRUE)
  value <- value_quad + reg * sum(T * log(pmax(T, 1e-300)))
  c(common, list(
    plan = T, value_quad = value_quad, value = value,
    status = "experimental",
    converged = FALSE,
    termination_reason = "experimental_no_convergence_certificate",
    certification = "experimental_no_convergence_claim",
    underlying_solver_status = out$status,
    underlying_solver_converged = out$converged
  ))
}

#' Deprecated pseudo-low-rank GW compatibility wrapper
#'
#' `lowrank_gromov_wasserstein_samples()` was misnamed: it always materializes
#' and solves a dense GW problem before applying SVD. Use
#' [dense_gromov_wasserstein_plan_svd()] for the same explicit operation.
#' POT-shaped factorized-cost, Dykstra, and initialization parameters were never
#' implemented and now error when changed from their legacy defaults.
#'
#' @inheritParams dense_gromov_wasserstein_plan_svd
#' @param alpha Deprecated unsupported POT-compatibility parameter.
#' @param gamma_init Deprecated unsupported POT-compatibility parameter.
#' @param cost_factorized_Xs Deprecated unsupported factorized source cost.
#' @param cost_factorized_Xt Deprecated unsupported factorized target cost.
#' @param stopThr_dykstra Deprecated unsupported Dykstra tolerance.
#' @param numItermax_dykstra Deprecated unsupported Dykstra iteration budget.
#' @param seed_init Deprecated unsupported initialization seed.
#' @param warn Deprecated unsupported warning control.
#' @param warn_dykstra Deprecated unsupported Dykstra warning control.
#' @return The result of [dense_gromov_wasserstein_plan_svd()].
#' @export
lowrank_gromov_wasserstein_samples <- function(
    X_s,
    X_t,
    a = NULL,
    b = NULL,
    reg = 0,
    rank = NULL,
    alpha = 1e-10,
    gamma_init = "rescale",
    rescale_cost = TRUE,
    cost_factorized_Xs = NULL,
    cost_factorized_Xt = NULL,
    stopThr = 1e-4,
    numItermax = 1000L,
    stopThr_dykstra = 1e-3,
    numItermax_dykstra = 10000L,
    seed_init = 49,
    warn = TRUE,
    warn_dykstra = FALSE,
    log = FALSE) {
  unsupported <- character()
  if (!missing(alpha)) unsupported <- c(unsupported, "alpha")
  if (!missing(gamma_init)) unsupported <- c(unsupported, "gamma_init")
  if (!missing(cost_factorized_Xs)) unsupported <- c(unsupported, "cost_factorized_Xs")
  if (!missing(cost_factorized_Xt)) unsupported <- c(unsupported, "cost_factorized_Xt")
  if (!missing(stopThr_dykstra)) unsupported <- c(unsupported, "stopThr_dykstra")
  if (!missing(numItermax_dykstra)) unsupported <- c(unsupported, "numItermax_dykstra")
  if (!missing(seed_init)) unsupported <- c(unsupported, "seed_init")
  if (!missing(warn)) unsupported <- c(unsupported, "warn")
  if (!missing(warn_dykstra)) unsupported <- c(unsupported, "warn_dykstra")
  if (length(unsupported)) {
    stop(
      sprintf(
        paste0(
          "Unsupported pseudo-low-rank compatibility parameter(s): %s. ",
          "They were never implemented; migrate to ",
          "`dense_gromov_wasserstein_plan_svd()`."
        ),
        paste(unsupported, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  .Deprecated("dense_gromov_wasserstein_plan_svd")
  dense_gromov_wasserstein_plan_svd(
    X_s = X_s,
    X_t = X_t,
    a = a,
    b = b,
    reg = reg,
    rank = rank,
    rescale_cost = rescale_cost,
    stopThr = stopThr,
    numItermax = numItermax,
    log = log
  )
}

#' Fused Unbalanced Across-Spaces Divergence (KL Sinkhorn)
#'
#' POT-compatible across-spaces unbalanced divergence solver (sample/feature
#' joint alignment) for `divergence = "kl"` with Sinkhorn updates.
#'
#' @param X Source matrix (`n_sample_x x n_feature_x`).
#' @param Y Target matrix (`n_sample_y x n_feature_y`).
#' @param wx_samp Source sample weights.
#' @param wx_feat Source feature weights.
#' @param wy_samp Target sample weights.
#' @param wy_feat Target feature weights.
#' @param reg_marginals Marginal relaxation(s), length 1 or 2.
#' @param epsilon Entropic regularization(s), scalar or length 2. Must be
#'   positive; default `1e-2`.
#' @param reg_type Either `"joint"` (FUGW-style) or `"independent"` (UCOOT).
#' @param divergence Only `"kl"` is supported. `"l2"` is rejected.
#' @param unbalanced_solver `"sinkhorn"` is the supported scaling-domain
#'   implementation. The former `"sinkhorn_log"` scaling alias is deprecated
#'   and errors because it was not a genuine log-domain solver. `"mm"` and
#'   `"lbfgsb"` are also rejected.
#' @param alpha Linear-term coefficient(s), scalar or length 2.
#' @param M_samp Optional sample linear cost matrix.
#' @param M_feat Optional feature linear cost matrix.
#' @param rescale_plan Rescale sample/feature plans to equal mass each BCD step.
#' @param init_pi Optional list with `pi_samp` and `pi_feat` initial couplings.
#' @param init_duals Accepted for POT-shaped signatures and ignored.
#' @param max_iter Max BCD iterations.
#' @param tol BCD stopping tolerance on sample coupling change.
#' @param max_iter_ot Max iterations for inner unbalanced Sinkhorn solves.
#' @param tol_ot Inner unbalanced Sinkhorn tolerance.
#' @param log If `TRUE`, return diagnostics and objective decomposition.
#' @param verbose If `TRUE`, print BCD diagnostics.
#' @return A list with `pi_samp`, `pi_feat`, `status`, and `converged`. If
#'   `log = TRUE`, also includes `error`, `linear_cost`, and `ucoot_cost`.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
fused_unbalanced_across_spaces_divergence <- function(
    X,
    Y,
    wx_samp = NULL,
    wx_feat = NULL,
    wy_samp = NULL,
    wy_feat = NULL,
    reg_marginals = 10,
    epsilon = 1e-2,
    reg_type = c("joint", "independent"),
    divergence = c("kl"),
    unbalanced_solver = c("sinkhorn", "sinkhorn_log"),
    alpha = 0,
    M_samp = NULL,
    M_feat = NULL,
    rescale_plan = TRUE,
    init_pi = NULL,
    init_duals = NULL,
    max_iter = 100L,
    tol = 1e-7,
    max_iter_ot = 500L,
    tol_ot = 1e-7,
    log = FALSE,
    verbose = FALSE,
    ...) {
  .reject_unused_dots(...)
  X <- .validate_finite_matrix(X, "X")
  Y <- .validate_finite_matrix(Y, "Y")
  max_iter <- .validate_count(max_iter, "max_iter")
  max_iter_ot <- .validate_count(max_iter_ot, "max_iter_ot")
  tol <- .validate_nonneg_scalar(tol, "tol")
  tol_ot <- .validate_positive_scalar(tol_ot, "tol_ot")
  reg_type <- match.arg(reg_type)
  divergence <- as.character(divergence)[1]
  unbalanced_solver <- as.character(unbalanced_solver)[1]
  if (!identical(divergence, "kl")) {
    stop(
      "`divergence = \"l2\"` is unsupported. Use `divergence = \"kl\"`.",
      call. = FALSE
    )
  }
  if (!unbalanced_solver %in% c("sinkhorn", "sinkhorn_log")) {
    stop(
      paste(
        "Unsupported `unbalanced_solver`.",
        "Supported choices: \"sinkhorn\", \"sinkhorn_log\".",
        "\"mm\" and \"lbfgsb\" are not implemented."
      ),
      call. = FALSE
    )
  }
  if (identical(unbalanced_solver, "sinkhorn_log")) {
    stop(
      paste(
        "`unbalanced_solver = \"sinkhorn_log\"` is deprecated and unsupported:",
        "the previous implementation was a scaling-domain alias.",
        "Use `unbalanced_solver = \"sinkhorn\"` within its documented regime."
      ),
      call. = FALSE
    )
  }

  nx_samp <- nrow(X)
  nx_feat <- ncol(X)
  ny_samp <- nrow(Y)
  ny_feat <- ncol(Y)

  if (is.null(wx_samp)) wx_samp <- rep(1 / nx_samp, nx_samp)
  if (is.null(wx_feat)) wx_feat <- rep(1 / nx_feat, nx_feat)
  if (is.null(wy_samp)) wy_samp <- rep(1 / ny_samp, ny_samp)
  if (is.null(wy_feat)) wy_feat <- rep(1 / ny_feat, ny_feat)
  wx_samp <- .assert_prob(wx_samp, nx_samp, "wx_samp")
  wx_feat <- .assert_prob(wx_feat, nx_feat, "wx_feat")
  wy_samp <- .assert_prob(wy_samp, ny_samp, "wy_samp")
  wy_feat <- .assert_prob(wy_feat, ny_feat, "wy_feat")

  reg_marginals <- .parse_pair(reg_marginals, "reg_marginals")
  epsilon <- .parse_pair(epsilon, "epsilon")
  alpha <- .parse_pair(alpha, "alpha")
  if (any(reg_marginals <= 0)) {
    stop("`reg_marginals` must contain positive values.", call. = FALSE)
  }

  rho_x <- reg_marginals[[1]]
  rho_y <- reg_marginals[[2]]
  eps_samp <- epsilon[[1]]
  eps_feat <- epsilon[[2]]

  if (identical(reg_type, "joint")) {
    eps_feat <- eps_samp
  }
  if (eps_samp <= 0 || eps_feat <= 0) {
    stop("Current KL Sinkhorn implementation requires positive `epsilon` values.", call. = FALSE)
  }

  if (!is.null(M_samp)) {
    M_samp <- .validate_finite_matrix(M_samp, "M_samp")
    if (nrow(M_samp) != nx_samp || ncol(M_samp) != ny_samp) {
      stop("`M_samp` has incompatible shape.", call. = FALSE)
    }
    M_samp <- alpha[[1]] * M_samp
  }

  if (!is.null(M_feat)) {
    M_feat <- .validate_finite_matrix(M_feat, "M_feat")
    if (nrow(M_feat) != nx_feat || ncol(M_feat) != ny_feat) {
      stop("`M_feat` has incompatible shape.", call. = FALSE)
    }
    M_feat <- alpha[[2]] * M_feat
  }

  wxy_samp <- wx_samp %o% wy_samp
  wxy_feat <- wx_feat %o% wy_feat

  if (is.null(init_pi)) {
    pi_samp <- wxy_samp
    pi_feat <- wxy_feat
  } else {
    if (!is.list(init_pi) || !all(c("pi_samp", "pi_feat") %in% names(init_pi))) {
      stop("`init_pi` must be NULL or a list containing `pi_samp` and `pi_feat`.", call. = FALSE)
    }
    pi_samp <- init_pi$pi_samp
    pi_feat <- init_pi$pi_feat
    .assert_matrix(pi_samp, "init_pi$pi_samp")
    .assert_matrix(pi_feat, "init_pi$pi_feat")
    if (nrow(pi_samp) != nx_samp || ncol(pi_samp) != ny_samp) {
      stop("`init_pi$pi_samp` has incompatible shape.", call. = FALSE)
    }
    if (nrow(pi_feat) != nx_feat || ncol(pi_feat) != ny_feat) {
      stop("`init_pi$pi_feat` has incompatible shape.", call. = FALSE)
    }
  }

  if (exists("cpp_ucoot_kl", mode = "function")) {
    out <- cpp_ucoot_kl(
      X = X,
      Y = Y,
      wx_samp = wx_samp,
      wx_feat = wx_feat,
      wy_samp = wy_samp,
      wy_feat = wy_feat,
      reg_marginals = c(rho_x, rho_y),
      epsilon = c(eps_samp, eps_feat),
      M_samp = if (is.null(M_samp)) .empty_feature_cost() else M_samp,
      M_feat = if (is.null(M_feat)) .empty_feature_cost() else M_feat,
      init_pi_samp = pi_samp,
      init_pi_feat = pi_feat,
      joint = identical(reg_type, "joint"),
      rescale_plan = isTRUE(rescale_plan),
      max_iter = as.integer(max_iter),
      tol = tol,
      max_iter_ot = as.integer(max_iter_ot),
      tol_ot = tol_ot,
      use_warm_start = TRUE
    )
    residual <- if (length(out$err_trace)) out$err_trace[length(out$err_trace)] else out$error
    ucoot_recomputed <- .fused_unbalanced_across_spaces_cost_kl(
      M_linear = list(M_samp, M_feat),
      data = list(X^2, Y^2, X, Y),
      tuple_pxy_samp = list(wx_samp, wy_samp, wxy_samp),
      tuple_pxy_feat = list(wx_feat, wy_feat, wxy_feat),
      pi_samp = out$pi_samp,
      pi_feat = out$pi_feat,
      hyperparams = c(rho_x, rho_y, eps_samp, eps_feat),
      reg_type = reg_type
    )
    if (isTRUE(log)) {
      out$error <- out$err_trace
    }
    ans <- .attach_solver_diagnostics(
      out,
      residual = residual,
      converged = is.finite(residual) && residual < tol,
      iterations = out$iterations,
      max_iter = as.integer(max_iter),
      plan = out$pi_samp,
      inner_residual = out$inner_residual,
      max_inner_residual = out$max_inner_residual,
      inner_iterations = if (!is.null(out$inner_iters_total)) {
        as.integer(out$inner_iters_total)
      } else {
        NA_integer_
      },
      inner_converged = out$inner_converged,
      inner_status = out$inner_status,
      feasibility = "unbalanced",
      feasibility_tol = tol_ot,
      objective_recomputed = ucoot_recomputed$ucoot_cost,
      objective_components = list(linear_cost = ucoot_recomputed$linear_cost)
    )
    ans$termination_reason <- .termination_reason_from_result(ans, as.integer(max_iter))
    if (!isTRUE(log)) {
      ans$linear_cost <- NULL
      ans$ucoot_cost <- NULL
      ans$error <- NULL
      ans$err_trace <- NULL
    }
    return(ans)
  }

  out <- .ucoot_kl_r_core(
    X = X,
    Y = Y,
    wx_samp = wx_samp,
    wx_feat = wx_feat,
    wy_samp = wy_samp,
    wy_feat = wy_feat,
    wxy_samp = wxy_samp,
    wxy_feat = wxy_feat,
    rho_x = rho_x,
    rho_y = rho_y,
    eps_samp = eps_samp,
    eps_feat = eps_feat,
    M_samp = M_samp,
    M_feat = M_feat,
    pi_samp = pi_samp,
    pi_feat = pi_feat,
    reg_type = reg_type,
    rescale_plan = rescale_plan,
    max_iter = max_iter,
    tol = tol,
    max_iter_ot = max_iter_ot,
    tol_ot = tol_ot,
    log = log,
    verbose = verbose
  )
  residual <- if (length(out$err_trace)) out$err_trace[length(out$err_trace)] else Inf
  if (!isTRUE(log)) {
    out$error <- NULL
    out$err_trace <- NULL
    out$linear_cost <- NULL
    out$ucoot_cost <- NULL
    out$duals_sample <- NULL
    out$duals_feature <- NULL
  } else {
    out$error <- out$err_trace
  }
  return(.attach_solver_diagnostics(
    out,
    residual = residual,
    converged = is.finite(residual) && residual < tol,
    iterations = length(out$err_trace),
    max_iter = as.integer(max_iter),
    plan = out$pi_samp
  ))
}

.ucoot_kl_r_core <- function(
    X,
    Y,
    wx_samp,
    wx_feat,
    wy_samp,
    wy_feat,
    wxy_samp,
    wxy_feat,
    rho_x,
    rho_y,
    eps_samp,
    eps_feat,
    M_samp,
    M_feat,
    pi_samp,
    pi_feat,
    reg_type,
    rescale_plan,
    max_iter,
    tol,
    max_iter_ot,
    tol_ot,
    log = TRUE,
    verbose = FALSE) {
  X_sqr <- X^2
  Y_sqr <- Y^2

  data_samp <- list(
    X_sqr = X_sqr,
    Y_sqr = Y_sqr,
    X = X,
    Y = Y,
    Y_t = t(Y),
    M = M_samp
  )
  data_feat <- list(
    X_sqr = t(X_sqr),
    Y_sqr = t(Y_sqr),
    X = t(X),
    Y = t(Y),
    Y_t = Y,
    M = M_feat
  )

  err_trace <- numeric()
  duals_samp <- NULL
  duals_feat <- NULL

  for (it in seq_len(as.integer(max_iter))) {
    pi_samp_prev <- pi_samp

    mass <- sum(pi_samp)
    uot_cost_feat <- .uot_cost_matrix_kl(
      data = data_feat,
      pi = pi_samp,
      tuple_p = list(wx_samp, wy_samp),
      hyperparams = c(rho_x, rho_y, eps_samp),
      reg_type = reg_type
    )

    new_rho <- c(rho_x * mass, rho_y * mass)
    new_eps <- if (identical(reg_type, "joint")) mass * eps_feat else eps_feat

    res_feat <- .sinkhorn_unbalanced_kl(
      M = uot_cost_feat,
      a = wx_feat,
      b = wy_feat,
      reg = new_eps,
      rho = new_rho,
      c = wxy_feat,
      max_iter = as.integer(max_iter_ot),
      tol = tol_ot,
      plan_init = pi_feat
    )
    pi_feat <- res_feat$plan
    duals_feat <- res_feat$potentials

    if (isTRUE(rescale_plan)) {
      pi_feat <- pi_feat * sqrt(mass / pmax(sum(pi_feat), 1e-300))
    }

    mass <- sum(pi_feat)
    uot_cost_samp <- .uot_cost_matrix_kl(
      data = data_samp,
      pi = pi_feat,
      tuple_p = list(wx_feat, wy_feat),
      hyperparams = c(rho_x, rho_y, eps_feat),
      reg_type = reg_type
    )

    new_rho <- c(rho_x * mass, rho_y * mass)
    new_eps <- if (identical(reg_type, "joint")) mass * eps_feat else eps_feat

    res_samp <- .sinkhorn_unbalanced_kl(
      M = uot_cost_samp,
      a = wx_samp,
      b = wy_samp,
      reg = new_eps,
      rho = new_rho,
      c = wxy_samp,
      max_iter = as.integer(max_iter_ot),
      tol = tol_ot,
      plan_init = pi_samp
    )
    pi_samp <- res_samp$plan
    duals_samp <- res_samp$potentials

    if (isTRUE(rescale_plan)) {
      pi_samp <- pi_samp * sqrt(mass / pmax(sum(pi_samp), 1e-300))
    }

    err <- sum(abs(pi_samp - pi_samp_prev))
    err_trace <- c(err_trace, err)

    if (isTRUE(verbose) && (it %% 10L == 0L || it == 1L)) {
      cat(sprintf("iter=%d err=%.8e\n", it, err))
    }

    if (err < tol) {
      break
    }
  }

  if (any(!is.finite(pi_samp)) || any(!is.finite(pi_feat))) {
    stop("Encountered non-finite values in coupling matrices. Adjust hyperparameters.", call. = FALSE)
  }

  out <- list(pi_samp = pi_samp, pi_feat = pi_feat, err_trace = err_trace)
  if (isTRUE(log)) {
    costs <- .fused_unbalanced_across_spaces_cost_kl(
      M_linear = list(M_samp, M_feat),
      data = list(X_sqr, Y_sqr, X, Y),
      tuple_pxy_samp = list(wx_samp, wy_samp, wxy_samp),
      tuple_pxy_feat = list(wx_feat, wy_feat, wxy_feat),
      pi_samp = pi_samp,
      pi_feat = pi_feat,
      hyperparams = c(rho_x, rho_y, eps_samp, eps_feat),
      reg_type = reg_type
    )
    out$duals_sample <- duals_samp
    out$duals_feature <- duals_feat
    out$linear_cost <- costs$linear_cost
    out$ucoot_cost <- costs$ucoot_cost
  }
  out
}

#' Unbalanced Co-Optimal Transport
#'
#' POT-compatible UCOOT wrapper (`reg_type = "independent"`).
#'
#' @inheritParams fused_unbalanced_across_spaces_divergence
#' @return A list with sample and feature couplings; with `log = TRUE`, includes
#'   objective diagnostics.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
unbalanced_co_optimal_transport <- function(
    X,
    Y,
    wx_samp = NULL,
    wx_feat = NULL,
    wy_samp = NULL,
    wy_feat = NULL,
    reg_marginals = 10,
    epsilon = 1e-2,
    divergence = c("kl"),
    unbalanced_solver = c("sinkhorn", "sinkhorn_log"),
    alpha = 0,
    M_samp = NULL,
    M_feat = NULL,
    rescale_plan = TRUE,
    init_pi = NULL,
    init_duals = NULL,
    max_iter = 100L,
    tol = 1e-7,
    max_iter_ot = 500L,
    tol_ot = 1e-7,
    log = FALSE,
    verbose = FALSE,
    ...) {
  fused_unbalanced_across_spaces_divergence(
    X = X,
    Y = Y,
    wx_samp = wx_samp,
    wx_feat = wx_feat,
    wy_samp = wy_samp,
    wy_feat = wy_feat,
    reg_marginals = reg_marginals,
    epsilon = epsilon,
    reg_type = "independent",
    divergence = divergence,
    unbalanced_solver = unbalanced_solver,
    alpha = alpha,
    M_samp = M_samp,
    M_feat = M_feat,
    rescale_plan = rescale_plan,
    init_pi = init_pi,
    init_duals = init_duals,
    max_iter = max_iter,
    tol = tol,
    max_iter_ot = max_iter_ot,
    tol_ot = tol_ot,
    log = log,
    verbose = verbose,
    ...
  )
}

#' Unbalanced Co-Optimal Transport Objective Value
#'
#' @inheritParams unbalanced_co_optimal_transport
#' @return Numeric UCOOT objective value. If `log = TRUE`, returns a list with
#'   `ucoot` and detailed diagnostics.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
unbalanced_co_optimal_transport2 <- function(..., log = FALSE) {
  out <- unbalanced_co_optimal_transport(..., log = TRUE)
  if (isTRUE(log)) {
    return(list(ucoot = out$ucoot_cost, log = out))
  }
  out$ucoot_cost
}
