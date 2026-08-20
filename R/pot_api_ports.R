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

.entropic_partial_wasserstein <- function(
    a,
    b,
    M,
    reg,
    m,
    numItermax = 1000L,
    stopThr = 1e-100,
    verbose = FALSE,
    log = FALSE) {
  if (exists("cpp_entropic_partial_wasserstein", mode = "function")) {
    return(cpp_entropic_partial_wasserstein(
      a = as.numeric(a),
      b = as.numeric(b),
      M = M,
      reg = as.numeric(reg),
      m = as.numeric(m),
      numItermax = as.integer(numItermax),
      stopThr = as.numeric(stopThr),
      verbose = isTRUE(verbose),
      log = isTRUE(log)
    ))
  }
  .assert_matrix(M, "M")
  if (!is.finite(reg) || reg <= 0) {
    stop("`reg` must be positive.", call. = FALSE)
  }

  ns <- nrow(M)
  nt <- ncol(M)
  if (length(a) != ns || length(b) != nt) {
    stop("`a` and `b` must match the dimensions of `M`.", call. = FALSE)
  }

  dx <- rep(1, ns)
  dy <- rep(1, nt)

  K <- exp(-M / reg)
  K <- K * (m / (sum(K) + 1e-300))

  q1 <- matrix(1, nrow = ns, ncol = nt)
  q2 <- matrix(1, nrow = ns, ncol = nt)
  q3 <- matrix(1, nrow = ns, ncol = nt)

  err <- Inf
  cpt <- 0L
  log_e <- list(err = numeric())

  while (err > stopThr && cpt < as.integer(numItermax)) {
    Kprev <- K

    K <- K * q1
    row_scale <- pmin(a / pmax(rowSums(K), 1e-300), dx)
    K1 <- sweep(K, 1L, row_scale, `*`)
    q1 <- q1 * (Kprev / pmax(K1, 1e-300))

    K1prev <- K1
    K1 <- K1 * q2
    col_scale <- pmin(b / pmax(colSums(K1), 1e-300), dy)
    K2 <- sweep(K1, 2L, col_scale, `*`)
    q2 <- q2 * (K1prev / pmax(K2, 1e-300))

    K2prev <- K2
    K2 <- K2 * q3
    K <- K2 * (m / pmax(sum(K2), 1e-300))
    q3 <- q3 * (K2prev / pmax(K, 1e-300))

    if (any(!is.finite(K))) {
      warning(sprintf("Numerical instability in entropic partial OT at iteration %d.", cpt), call. = FALSE)
      break
    }

    if ((cpt %% 10L) == 0L) {
      err <- sqrt(sum((Kprev - K)^2))
      if (isTRUE(log)) {
        log_e$err <- c(log_e$err, err)
      }
      if (isTRUE(verbose) && ((cpt %% 200L) == 0L || cpt == 0L)) {
        cat(sprintf("it=%d err=%.8e\n", cpt, err))
      }
    }

    cpt <- cpt + 1L
  }

  if (isTRUE(log)) {
    log_e$partial_w_dist <- sum(M * K)
    return(list(plan = K, log = log_e))
  }
  K
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
    check_every = 10L) {
  lin_w <- 1 - alpha
  quad_w <- alpha

  G <- G0
  err <- Inf
  it <- 0L
  err_trace <- numeric()

  for (k in seq_len(max_iter)) {
    Gprev <- G
    gw <- .gw_square_terms(C1, C2, G, symmetric = symmetric)
    M_entr <- quad_w * gw$grad + lin_w * M

    G <- .entropic_partial_wasserstein(
      a = p,
      b = q,
      M = M_entr,
      reg = reg,
      m = m,
      numItermax = inner_max_iter,
      stopThr = inner_tol,
      verbose = FALSE,
      log = FALSE
    )

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
    err_trace = err_trace
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
    check_every = 10L) {
  if (exists("cpp_partial_fgw_entropic_square", mode = "function")) {
    return(cpp_partial_fgw_entropic_square(
      M = M,
      C1 = C1,
      C2 = C2,
      p = p,
      q = q,
      m = m,
      reg = reg,
      alpha = alpha,
      symmetric = isTRUE(symmetric),
      init_plan = if (is.null(G0)) .empty_feature_cost() else G0,
      max_iter = as.integer(max_iter),
      tol = tol,
      inner_max_iter = as.integer(inner_max_iter),
      inner_tol = inner_tol,
      check_every = as.integer(check_every)
    ))
  }
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
    check_every = check_every
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

.partial_sinkhorn_dispatch <- function(
    method,
    reg,
    C1,
    C2,
    M = NULL,
    alpha = 1) {
  structure_bound <- 2 * alpha * (max(abs(C1)) + max(abs(C2)))^2
  feature_bound <- if (is.null(M)) 0 else (1 - alpha) * max(abs(M))
  proxy <- c(0, structure_bound + feature_bound, -(structure_bound + feature_bound))
  dispatch <- .select_sinkhorn_method(
    method,
    proxy,
    reg,
    precision = "double",
    context = "Entropic partial Sinkhorn"
  )
  if (identical(dispatch$effective, "log")) {
    detail <- if (identical(dispatch$requested, "log")) {
      "`method = \"log\"` was explicitly requested"
    } else {
      sprintf(
        "auto found metric %.6g above the scaling threshold %.6g",
        dispatch$metric, dispatch$threshold
      )
    }
    stop(
      paste0(
        "Entropic partial Sinkhorn requires a genuine log-domain Dykstra ",
        "backend because ", detail, "; that backend is not implemented. ",
        "Rescale the costs or increase `reg`."
      ),
      call. = FALSE
    )
  }
  dispatch
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
#' @param method Scaling-domain inner Dykstra method or `"auto"`. The public
#'   partial primitive is certified only when its conservative dynamic-range
#'   proxy is at most 500. `"log"` and unsafe auto requests error until a
#'   genuine log-domain Dykstra backend exists.
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
    method = c("scaling", "auto", "log")) {
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
  dispatch <- .partial_sinkhorn_dispatch(method, reg, C1, C2)

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
    check_every = check_every
  )

  if (!isTRUE(log)) {
    plan <- out$plan
    attr(plan, "sinkhorn_dispatch") <- dispatch
    return(plan)
  }

  objective_recomputed <- ot_gw_square(C1, C2, out$plan, symmetric = symmetric)
  ans <- .partial_log_with_nested_status(
    out, "partial_gw_dist", tol, numItermax,
    p, q, m, inner_tol, objective_recomputed,
    list(lin_loss = 0, quad_loss = objective_recomputed),
    mass_defaulted = mass_defaulted
  )
  .attach_partial_sinkhorn_dispatch(ans, dispatch)
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
    method = c("scaling", "auto", "log")) {
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
  dispatch <- .partial_sinkhorn_dispatch(
    method, reg, C1, C2, M = M, alpha = alpha
  )

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
    check_every = check_every
  )

  if (!isTRUE(log)) {
    plan <- out$plan
    attr(plan, "sinkhorn_dispatch") <- dispatch
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
  .attach_partial_sinkhorn_dispatch(ans, dispatch)
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

#' Low-rank GW from Samples (Experimental)
#'
#' Approximate low-rank sampled GW interface for large problems. This function
#' computes an entropic GW plan then derives a rank-`r` factorization via
#' truncated SVD. It is provided as an experimental POT-compatibility tier.
#'
#' @param X_s Source samples (`ns x d`).
#' @param X_t Target samples (`nt x d`).
#' @param a Optional source weights.
#' @param b Optional target weights.
#' @param reg Entropic regularization used for the inner GW solve.
#' @param rank Target low-rank factorization rank. Values below 1 error.
#'   Values above `min(ns, nt)` warn and clamp.
#' @param alpha Lower bound parameter retained for POT compatibility.
#' @param gamma_init Retained for POT compatibility.
#' @param rescale_cost Whether to normalize structure costs to `[0, 1]`.
#' @param cost_factorized_Xs Retained for POT compatibility.
#' @param cost_factorized_Xt Retained for POT compatibility.
#' @param stopThr Outer GW stopping tolerance.
#' @param numItermax Outer GW max iterations.
#' @param stopThr_dykstra Retained for POT compatibility.
#' @param numItermax_dykstra Retained for POT compatibility.
#' @param seed_init Optional random seed.
#' @param warn Retained for POT compatibility.
#' @param warn_dykstra Retained for POT compatibility.
#' @param log If `TRUE`, return diagnostics.
#' @return A list with `Q`, `R`, `g`. If `log = TRUE`, includes `plan`,
#'   `value_quad`, and `value`.
#'
#' @section Experimental:
#' This is a post-hoc SVD of a dense GW plan, not a factorized solver.
#' Reconstruction error `||T - T_r||_F / ||T||_F` decreases as rank
#' increases up to `min(ns, nt)`. It is experimental and may be renamed
#' (`bd-01M05QY9GKDDTB3CXTKAJE0E8G`). See `inst/bench/sampled-budget-curves.md`.
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
  numItermax <- .validate_count(numItermax, "numItermax")
  numItermax_dykstra <- .validate_count(
    numItermax_dykstra, "numItermax_dykstra"
  )
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

  if (!isTRUE(log)) {
    return(list(Q = Q, R = R, g = g))
  }

  value_quad <- .gw_square_value(C1, C2, T, a, b, symmetric = TRUE)
  value <- value_quad + reg * sum(T * log(pmax(T, 1e-300)))
  list(
    Q = Q, R = R, g = g, plan = T, value_quad = value_quad, value = value,
    status = "experimental",
    converged = FALSE,
    termination_reason = "experimental_no_convergence_certificate",
    certification = "experimental_no_convergence_claim",
    underlying_solver_status = out$status,
    underlying_solver_converged = out$converged
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
