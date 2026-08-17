.assert_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix.", name), call. = FALSE)
  }
}

.assert_prob <- function(x, n, name) {
  if (!is.numeric(x) || length(x) != n) {
    stop(sprintf("`%s` must be numeric with length %d.", name, n), call. = FALSE)
  }
  if (any(!is.finite(x)) || any(x < 0)) {
    stop(sprintf("`%s` must be finite and nonnegative.", name), call. = FALSE)
  }
  s <- sum(x)
  if (s <= 0) {
    stop(sprintf("`%s` must have positive total mass.", name), call. = FALSE)
  }
  x / s
}

#' Entropic Fused Gromov-Wasserstein (square loss)
#'
#' POT-style projected gradient iterations with Sinkhorn projections.
#' Reported `fgw_dist` is the unregularized FGW objective at the returned plan;
#' it does not add the entropic term.
#'
#' @param M Cross-domain feature cost matrix (`ns x nt`).
#' @param C1 Source structure cost matrix (`ns x ns`).
#' @param C2 Target structure cost matrix (`nt x nt`).
#' @param p Source weights (default uniform). Renormalized to sum 1.
#' @param q Target weights (default uniform). Renormalized to sum 1.
#' @param alpha Trade-off between feature and structure terms, in `[0, 1]`.
#' @param epsilon Entropic regularization (must be positive).
#' @param max_iter Max outer FGW iterations.
#' @param tol Outer stopping tolerance on Frobenius norm of plan updates.
#' @param sinkhorn_max_iter Max Sinkhorn iterations per outer step.
#' @param sinkhorn_tol Sinkhorn stopping tolerance on marginal residual.
#' @param init_plan Optional initial coupling (`ns x nt`) used as a warm start.
#' @param structure_rank Optional low-rank approximation rank for structure
#'   matrices in the FGW tensor product (`0` disables approximation).
#' @param sinkhorn_method Sinkhorn variant: `"scaling"` or `"log"`.
#' @param precision Numeric precision mode: `"mixed"` (default) or `"double"`
#'   (`"mixed"` uses float inner iterations with double final objective
#'   evaluation). For larger problems, `"double"` uses a mixed-precision
#'   accelerated inner Sinkhorn path while returning double outputs.
#' @param symmetric `NULL` auto-detects symmetry of `C1` and `C2`. `TRUE`
#'   requires both costs to be symmetric within `1e-10`. `FALSE` uses the
#'   two-sided tensor.
#' @param solver Either `"PGD"` or `"PPA"`.
#' @param check_every Evaluate outer stopping error every `check_every` iterations.
#' @return A list with `plan`, `fgw_dist`, `iterations`, `error`, `residual`,
#'   `status`, `converged`, and marginal residuals. See
#'   `inst/solver-contract.md`.
#' @examples
#' set.seed(1)
#' C1 <- as.matrix(dist(matrix(rnorm(8), 4, 2)))
#' C2 <- as.matrix(dist(matrix(rnorm(10), 5, 2)))
#' M <- matrix(runif(20), 4, 5)
#' out <- fgw_entropic(M, C1, C2, epsilon = 0.1, max_iter = 40L)
#' out$status
#' out$residual
#' @export
fgw_entropic <- function(
    M,
    C1,
    C2,
    p = NULL,
    q = NULL,
    alpha = 0.5,
    epsilon = 0.1,
    max_iter = 1000L,
    tol = 1e-9,
    sinkhorn_max_iter = 1000L,
    sinkhorn_tol = 1e-9,
    init_plan = NULL,
    structure_rank = 0L,
    sinkhorn_method = c("scaling", "log"),
    precision = c("mixed", "double"),
    symmetric = NULL,
    solver = c("PGD", "PPA"),
    check_every = 1L) {
  M <- .validate_finite_matrix(M, "M")
  C1 <- .validate_finite_matrix(C1, "C1", square = TRUE)
  C2 <- .validate_finite_matrix(C2, "C2", square = TRUE)

  ns <- nrow(C1)
  nt <- nrow(C2)
  if (nrow(M) != ns || ncol(M) != nt) {
    stop("`M` must have shape nrow(C1) x nrow(C2).", call. = FALSE)
  }

  alpha <- .validate_alpha(alpha)
  epsilon <- .validate_positive_scalar(epsilon, "epsilon")
  tol <- .validate_positive_scalar(tol, "tol")
  sinkhorn_tol <- .validate_positive_scalar(sinkhorn_tol, "sinkhorn_tol")
  max_iter <- .validate_count(max_iter, "max_iter")
  sinkhorn_max_iter <- .validate_count(sinkhorn_max_iter, "sinkhorn_max_iter")
  check_every <- .validate_count(check_every, "check_every")
  structure_rank <- .validate_count(structure_rank, "structure_rank", min = 0L)
  symmetric <- .resolve_symmetric(symmetric, C1, C2)

  if (is.null(p)) p <- rep(1 / ns, ns)
  if (is.null(q)) q <- rep(1 / nt, nt)
  p <- .assert_prob(p, ns, "p")
  q <- .assert_prob(q, nt, "q")
  init_plan <- .validate_optional_init_plan(init_plan, ns, nt, "init_plan")

  solver <- match.arg(solver)
  sinkhorn_method <- match.arg(sinkhorn_method)
  precision <- match.arg(precision)
  use_ppa <- identical(solver, "PPA")
  use_log_sinkhorn <- identical(sinkhorn_method, "log")
  use_mixed_precision <- identical(precision, "mixed")

  out <- cpp_fgw_entropic_square(
    M = unname(M),
    C1 = unname(C1),
    C2 = unname(C2),
    p = unname(p),
    q = unname(q),
    alpha = alpha,
    epsilon = epsilon,
    max_iter = max_iter,
    tol = tol,
    sinkhorn_max_iter = sinkhorn_max_iter,
    sinkhorn_tol = sinkhorn_tol,
    init_plan = unname(init_plan),
    symmetric = symmetric,
    use_ppa = use_ppa,
    use_log_sinkhorn = use_log_sinkhorn,
    use_mixed_precision = use_mixed_precision,
    check_every = check_every,
    approx_rank = structure_rank
  )
  residual <- if (!is.null(out$error)) out$error else Inf
  .attach_solver_diagnostics(
    out,
    residual = residual,
    converged = is.finite(residual) && residual <= tol,
    iterations = out$iterations,
    max_iter = max_iter,
    p = p,
    q = q,
    plan = out$plan
  )
}

#' Entropic Fused Gromov-Wasserstein Objective Value
#'
#' @inheritParams fgw_entropic
#' @return Numeric scalar FGW value.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
fgw_entropic2 <- function(...) {
  out <- fgw_entropic(...)
  out$fgw_dist
}

#' Fused Unbalanced Gromov-Wasserstein (KL divergence, Sinkhorn inner solver)
#'
#' POT-style BCD on two couplings with KL-unbalanced OT subproblems.
#'
#' @param Cx Source structure matrix.
#' @param Cy Target structure matrix.
#' @param wx Source weights (default uniform). Renormalized to sum 1.
#' @param wy Target weights (default uniform). Renormalized to sum 1.
#' @param reg_marginals Length-1 or length-2 marginal relaxation parameters.
#' @param epsilon Joint KL regularization parameter (entropic term).
#' @param alpha Linear FGW term coefficient, in `[0, 1]`.
#' @param M Linear sample cost matrix (default `NULL`).
#' @param init_pi Optional initial sample coupling (`nx x ny`).
#' @param max_iter Max BCD iterations.
#' @param tol BCD stopping tolerance (`L1` delta on sample coupling).
#' @param max_iter_ot Max unbalanced Sinkhorn iterations per inner solve.
#' @param tol_ot Inner Sinkhorn tolerance.
#' @param rescale_plan Whether to rescale successive plans to equal mass.
#' @param check_every Evaluate BCD stopping criterion every `check_every` iterations.
#' @param precision Numeric precision mode for the C++ solver (`"double"` or
#'   `"mixed"`).
#' @return A list with `pi_samp`, `pi_feat`, `fugw_cost`, `linear_cost`,
#'   `iterations`, `error`, `status`, and `converged`.
#' @examples
#' set.seed(1)
#' Cx <- as.matrix(dist(matrix(rnorm(8), 4, 2)))
#' Cy <- as.matrix(dist(matrix(rnorm(10), 5, 2)))
#' out <- fugw_kl(Cx, Cy, epsilon = 0.05, max_iter = 20L)
#' out$status
#' @export
fugw_kl <- function(
    Cx,
    Cy,
    wx = NULL,
    wy = NULL,
    reg_marginals = c(10, 10),
    epsilon = 1e-2,
    alpha = 0.5,
    M = NULL,
    init_pi = NULL,
    max_iter = 100L,
    tol = 1e-7,
    max_iter_ot = 500L,
    tol_ot = 1e-7,
    rescale_plan = TRUE,
    check_every = 1L,
    precision = c("double", "mixed")) {
  precision <- match.arg(precision)
  Cx <- .validate_finite_matrix(Cx, "Cx", square = TRUE)
  Cy <- .validate_finite_matrix(Cy, "Cy", square = TRUE)
  nx <- nrow(Cx)
  ny <- nrow(Cy)
  alpha <- .validate_alpha(alpha)
  epsilon <- .validate_positive_scalar(epsilon, "epsilon")
  tol <- .validate_positive_scalar(tol, "tol")
  tol_ot <- .validate_positive_scalar(tol_ot, "tol_ot")
  max_iter <- .validate_count(max_iter, "max_iter")
  max_iter_ot <- .validate_count(max_iter_ot, "max_iter_ot")
  check_every <- .validate_count(check_every, "check_every")

  if (is.null(wx)) wx <- rep(1 / nx, nx)
  if (is.null(wy)) wy <- rep(1 / ny, ny)
  wx <- .assert_prob(wx, nx, "wx")
  wy <- .assert_prob(wy, ny, "wy")

  if (is.null(M)) {
    M <- matrix(0, nrow = nx, ncol = ny)
  } else {
    M <- .validate_finite_matrix(M, "M")
    if (nrow(M) != nx || ncol(M) != ny) {
      stop("`M` must have shape nrow(Cx) x nrow(Cy).", call. = FALSE)
    }
  }

  if (length(reg_marginals) == 1L) {
    reg_marginals <- c(reg_marginals, reg_marginals)
  }
  if (length(reg_marginals) != 2L || any(!is.finite(reg_marginals)) || any(reg_marginals <= 0)) {
    stop("`reg_marginals` must be one or two finite positive numbers.", call. = FALSE)
  }

  if (is.null(init_pi)) {
    init_pi <- tcrossprod(wx, wy)
  } else {
    init_pi <- .validate_optional_init_plan(init_pi, nx, ny, "init_pi")
    if (length(init_pi) == 0L) {
      stop("`init_pi` must have shape nrow(Cx) x nrow(Cy).", call. = FALSE)
    }
  }

  out <- cpp_fugw_kl_square(
    Cx = unname(Cx),
    Cy = unname(Cy),
    wx = unname(wx),
    wy = unname(wy),
    reg_marginals = as.numeric(reg_marginals),
    epsilon = epsilon,
    alpha = alpha,
    M = unname(M),
    init_pi = unname(init_pi),
    max_iter = max_iter,
    tol = tol,
    max_iter_ot = max_iter_ot,
    tol_ot = tol_ot,
    rescale_plan = isTRUE(rescale_plan),
    check_every = check_every,
    use_mixed_precision = identical(precision, "mixed")
  )
  residual <- if (!is.null(out$error)) out$error else Inf
  inner_iters <- if (!is.null(out$inner_iters_total)) {
    as.integer(out$inner_iters_total)
  } else {
    NA_integer_
  }
  .attach_solver_diagnostics(
    out,
    residual = residual,
    converged = is.finite(residual) && residual <= tol,
    iterations = out$iterations,
    max_iter = max_iter,
    plan = out$pi_samp,
    inner_iterations = inner_iters
  )
}

#' Fused Unbalanced Gromov-Wasserstein Objective Value
#'
#' @inheritParams fugw_kl
#' @return Numeric scalar FUGW objective value.
#' @param ... Additional arguments. Unused extras are rejected when the solver uses `.reject_unused_dots()`; otherwise they are forwarded to the primary solver.
#' @export
fugw_kl2 <- function(...) {
  out <- fugw_kl(...)
  out$fugw_cost
}

.solve_1d_linesearch_quad <- function(a, b) {
  if (!is.finite(a) || !is.finite(b)) return(0)
  if (a > 0) {
    return(-b / (2 * a))
  }
  if ((a + b) < 0) return(1)
  0
}

.prepare_integer_marginals <- function(p, q, scale = 1e6) {
  n <- length(p)
  m <- length(q)
  mass <- as.integer(round(scale))
  if (mass <= 0L) {
    stop("`lp_scale` must be positive.", call. = FALSE)
  }

  a <- as.integer(round(p * mass))
  b <- as.integer(round(q * mass))

  da <- mass - sum(a)
  if (da != 0L) {
    ia <- which.max(a)
    a[ia] <- a[ia] + da
  }
  db <- mass - sum(b)
  if (db != 0L) {
    ib <- which.max(b)
    b[ib] <- b[ib] + db
  }

  dmass <- sum(a) - sum(b)
  if (dmass != 0L) {
    ib <- which.max(b)
    b[ib] <- b[ib] + dmass
  }
  if (any(a < 0L) || any(b < 0L)) {
    stop("Failed to build nonnegative integer marginals for LP solve.", call. = FALSE)
  }

  list(a = a, b = b, n = n, m = m, mass = sum(a))
}

.make_transport_lp_solver <- function(
    p,
    q,
    scale = 1e6,
    solver = c("lp_matrix", "lp_transport")) {
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop(
      paste(
        "`lp_solver = \"lp_transport\"` and `\"lp_matrix\"` require package `lpSolve`.",
        "Install it with install.packages(\"lpSolve\"), or use the default",
        "`lp_solver = \"cpp_transport\"` which has no extra dependency."
      ),
      call. = FALSE
    )
  }
  solver <- match.arg(solver)
  prep <- .prepare_integer_marginals(p, q, scale = scale)
  a <- prep$a
  b <- prep$b
  n <- prep$n
  m <- prep$m
  mass <- prep$mass

  # lp.transport rebuilds constraints each call; use lp() with cached constraints for speed.
  if (identical(solver, "lp_matrix")) {
    const.mat <- matrix(0, nrow = n + m, ncol = n * m)
    for (i in seq_len(n)) {
      idx <- ((i - 1L) * m + 1L):(i * m)
      const.mat[i, idx] <- 1
    }
    for (j in seq_len(m)) {
      idx <- seq.int(j, by = m, length.out = n)
      const.mat[n + j, idx] <- 1
    }
    const.dir <- rep("=", n + m)
    const.rhs <- c(a, b)

    function(cost) {
      sol <- lpSolve::lp(
        direction = "min",
        objective.in = as.vector(t(cost)),
        const.mat = const.mat,
        const.dir = const.dir,
        const.rhs = const.rhs
      )
      if (!identical(sol$status, 0L)) {
        stop(sprintf("lp() failed with status %s.", sol$status), call. = FALSE)
      }
      matrix(sol$solution, nrow = n, byrow = TRUE) / mass
    }
  } else {
    function(cost) {
      sol <- lpSolve::lp.transport(
        cost.mat = cost,
        direction = "min",
        row.signs = rep("=", n),
        row.rhs = a,
        col.signs = rep("=", m),
        col.rhs = b
      )
      if (!identical(sol$status, 0L)) {
        stop(sprintf("lp.transport failed with status %s.", sol$status), call. = FALSE)
      }
      sol$solution / mass
    }
  }
}

#' Unregularized FGW via Conditional Gradient + LP Direction (square loss)
#'
#' Conditional-gradient iterations with an exact linear-OT direction step.
#' The outer GW/FGW problem is non-convex; the returned plan is a stationary
#' point of this procedure, not a certified global minimizer. Reported
#' `fgw_dist` is the unregularized FGW objective.
#'
#' @inheritParams fgw_entropic
#' @param G0 Optional feasible initial coupling (`ns x nt`). Honored by every
#'   `lp_solver` backend.
#' @param tol_rel Relative stopping tolerance on objective change.
#' @param tol_abs Absolute stopping tolerance on objective change.
#' @param lp_scale Integer mass scaling used for LP marginals (lpSolve backends).
#' @param lp_solver LP backend for direction step: `"cpp_transport"` (default,
#'   C++ transport-simplex backend), `"lp_transport"` (lpSolve), or `"lp_matrix"`
#'   (cached dense lpSolve constraints).
#' @param lp_max_iter Maximum iterations for the C++ transport-simplex LP backend.
#' @param lp_tol Optimality tolerance for the C++ transport-simplex LP backend.
#' @return A list with `plan`, `fgw_dist`, `iterations`, `error`, `rel_error`,
#'   `loss_trace`, `status`, and `converged`.
#' @examples
#' set.seed(1)
#' C1 <- as.matrix(dist(matrix(rnorm(8), 4, 2)))
#' C2 <- as.matrix(dist(matrix(rnorm(8), 4, 2)))
#' M <- matrix(runif(16), 4, 4)
#' out <- fgw_exact_cg(M, C1, C2, max_iter = 20L)
#' out$status
#' @export
fgw_exact_cg <- function(
    M,
    C1,
    C2,
    p = NULL,
    q = NULL,
    alpha = 0.5,
    symmetric = NULL,
    G0 = NULL,
    max_iter = 500L,
    tol_rel = 1e-9,
    tol_abs = 1e-9,
    lp_scale = 1e6,
    lp_solver = c("cpp_transport", "lp_transport", "lp_matrix"),
    lp_max_iter = 20000L,
    lp_tol = 1e-12) {
  M <- .validate_finite_matrix(M, "M")
  C1 <- .validate_finite_matrix(C1, "C1", square = TRUE)
  C2 <- .validate_finite_matrix(C2, "C2", square = TRUE)

  ns <- nrow(C1)
  nt <- nrow(C2)
  if (nrow(M) != ns || ncol(M) != nt) {
    stop("`M` must have shape nrow(C1) x nrow(C2).", call. = FALSE)
  }

  alpha <- .validate_alpha(alpha)
  tol_rel <- .validate_positive_scalar(tol_rel, "tol_rel")
  tol_abs <- .validate_positive_scalar(tol_abs, "tol_abs")
  lp_tol <- .validate_positive_scalar(lp_tol, "lp_tol")
  max_iter <- .validate_count(max_iter, "max_iter", min = 0L)
  lp_max_iter <- .validate_count(lp_max_iter, "lp_max_iter")
  lp_scale <- .validate_positive_scalar(lp_scale, "lp_scale")
  symmetric <- .resolve_symmetric(symmetric, C1, C2)

  if (is.null(p)) p <- rep(1 / ns, ns)
  if (is.null(q)) q <- rep(1 / nt, nt)
  p <- .assert_prob(p, ns, "p")
  q <- .assert_prob(q, nt, "q")
  lp_solver <- match.arg(lp_solver)
  if (!is.null(G0)) {
    G0 <- .validate_balanced_plan(G0, p, q, "G0")
  }

  if (identical(lp_solver, "cpp_transport")) {
    init_plan <- if (is.null(G0)) {
      matrix(numeric(0), nrow = 0, ncol = 0)
    } else {
      unname(G0)
    }
    out <- cpp_fgw_exact_cg_square(
      M = unname(M),
      C1 = unname(C1),
      C2 = unname(C2),
      p = unname(p),
      q = unname(q),
      alpha = alpha,
      symmetric = symmetric,
      max_iter = max_iter,
      tol_rel = tol_rel,
      tol_abs = tol_abs,
      lp_max_iter = lp_max_iter,
      lp_tol = lp_tol,
      init_plan = init_plan
    )
    residual <- if (!is.null(out$error)) out$error else Inf
    rel_error <- if (!is.null(out$rel_error)) out$rel_error else Inf
    converged <- (is.finite(rel_error) && rel_error < tol_rel) ||
      (is.finite(residual) && residual < tol_abs)
    return(.attach_solver_diagnostics(
      out,
      residual = residual,
      converged = converged,
      iterations = out$iterations,
      max_iter = max_iter,
      p = p,
      q = q,
      plan = out$plan,
      inner_iterations = if (!is.null(out$inner_iterations)) out$inner_iterations else NA_integer_,
      lp_ok = isTRUE(out$lp_ok)
    ))
  }

  if (is.null(G0)) {
    G <- tcrossprod(p, q)
  } else {
    G <- G0
  }

  fC1 <- C1^2
  fC2 <- C2^2
  hC1 <- C1
  hC2 <- 2 * C2
  constC <- tcrossprod(as.vector(fC1 %*% p), rep(1, nt)) +
    tcrossprod(rep(1, ns), as.vector(q %*% t(fC2)))

  if (!symmetric) {
    fC1t <- t(C1)^2
    fC2t <- t(C2)^2
    hC1t <- t(C1)
    hC2t <- 2 * t(C2)
    constCt <- tcrossprod(as.vector(fC1t %*% p), rep(1, nt)) +
      tcrossprod(rep(1, ns), as.vector(q %*% t(fC2t)))
  }

  M_lin <- (1 - alpha) * M
  lp_direction <- .make_transport_lp_solver(p, q, scale = lp_scale, solver = lp_solver)

  Acur <- hC1 %*% G %*% t(hC2)
  tens <- constC - Acur
  if (!symmetric) {
    Acurt <- hC1t %*% G %*% t(hC2t)
    tenst <- constCt - Acurt
  }

  gw_term <- if (symmetric) sum(tens * G) else 0.5 * (sum(tens * G) + sum(tenst * G))
  cost_G <- sum(M_lin * G) + alpha * gw_term
  loss_trace <- numeric(max_iter + 1L)
  loss_trace[1] <- cost_G
  rel_delta <- NaN
  abs_delta <- NaN
  it <- 0L

  for (k in seq_len(max_iter)) {
    gG <- if (symmetric) {
      2 * tens
    } else {
      tens + tenst
    }
    Mi <- M_lin + alpha * gG
    Gc <- lp_direction(Mi)
    deltaG <- Gc - G

    dot <- hC1 %*% deltaG %*% t(hC2)
    a_ls <- -alpha * sum(dot * deltaG)
    if (symmetric) {
      b_ls <- sum(M_lin * deltaG) - 2 * alpha * sum(dot * G)
    } else {
      b_ls <- sum(M_lin * deltaG) - alpha * (
        sum(dot * G) + sum(Acur * deltaG)
      )
    }

    step <- .solve_1d_linesearch_quad(a_ls, b_ls)
    step <- min(1, max(0, step))
    new_cost <- cost_G + a_ls * step^2 + b_ls * step
    G <- G + step * deltaG
    Acur <- Acur + step * dot
    tens <- constC - Acur
    if (!symmetric) {
      dot_t <- hC1t %*% deltaG %*% t(hC2t)
      Acurt <- Acurt + step * dot_t
      tenst <- constCt - Acurt
    }

    abs_delta <- abs(new_cost - cost_G)
    rel_delta <- if (new_cost != 0) abs_delta / abs(new_cost) else NaN
    cost_G <- new_cost
    loss_trace[k + 1L] <- cost_G
    it <- k

    if ((is.finite(rel_delta) && rel_delta < tol_rel) || (is.finite(abs_delta) && abs_delta < tol_abs)) {
      break
    }
  }

  out <- list(
    plan = G,
    fgw_dist = cost_G,
    iterations = it,
    error = abs_delta,
    rel_error = rel_delta,
    loss_trace = loss_trace[seq_len(it + 1L)]
  )
  converged <- (is.finite(rel_delta) && rel_delta < tol_rel) ||
    (is.finite(abs_delta) && abs_delta < tol_abs)
  .attach_solver_diagnostics(
    out,
    residual = abs_delta,
    converged = converged,
    iterations = it,
    max_iter = max_iter,
    p = p,
    q = q,
    plan = G
  )
}
