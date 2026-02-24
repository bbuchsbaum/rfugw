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
#'
#' @param M Cross-domain feature cost matrix (`ns x nt`).
#' @param C1 Source structure cost matrix (`ns x ns`).
#' @param C2 Target structure cost matrix (`nt x nt`).
#' @param p Source weights (default uniform).
#' @param q Target weights (default uniform).
#' @param alpha Trade-off between feature and structure terms.
#' @param epsilon Entropic regularization.
#' @param max_iter Max outer FGW iterations.
#' @param tol Outer stopping tolerance on Frobenius norm of plan updates.
#' @param sinkhorn_max_iter Max Sinkhorn iterations per outer step.
#' @param sinkhorn_tol Sinkhorn stopping tolerance on marginal residual.
#' @param init_plan Optional initial coupling (`ns x nt`) used as a warm start.
#' @param structure_rank Optional low-rank approximation rank for structure
#'   matrices in the FGW tensor product (`0` disables approximation).
#' @param precision Numeric precision mode: `"mixed"` (default) or `"double"`
#'   (`"mixed"` uses float inner iterations with double final objective
#'   evaluation). For larger problems, `"double"` uses a mixed-precision
#'   accelerated inner Sinkhorn path while returning double outputs.
#' @param symmetric Assume `C1` and `C2` are symmetric if `TRUE`.
#' @param solver Either `"PGD"` or `"PPA"`.
#' @param check_every Evaluate outer stopping error every `check_every` iterations.
#' @return A list with `plan`, `fgw_dist`, `iterations`, `error`.
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
    symmetric = TRUE,
    solver = c("PGD", "PPA"),
    check_every = 1L) {
  .assert_matrix(M, "M")
  .assert_matrix(C1, "C1")
  .assert_matrix(C2, "C2")

  ns <- nrow(C1)
  nt <- nrow(C2)
  if (ncol(C1) != ns) stop("`C1` must be square.", call. = FALSE)
  if (ncol(C2) != nt) stop("`C2` must be square.", call. = FALSE)
  if (nrow(M) != ns || ncol(M) != nt) {
    stop("`M` must have shape nrow(C1) x nrow(C2).", call. = FALSE)
  }

  if (is.null(p)) p <- rep(1 / ns, ns)
  if (is.null(q)) q <- rep(1 / nt, nt)
  p <- .assert_prob(p, ns, "p")
  q <- .assert_prob(q, nt, "q")
  if (is.null(init_plan)) {
    init_plan <- matrix(numeric(0), nrow = 0, ncol = 0)
  } else {
    .assert_matrix(init_plan, "init_plan")
    if (nrow(init_plan) != ns || ncol(init_plan) != nt) {
      stop("`init_plan` must have shape nrow(C1) x nrow(C2).", call. = FALSE)
    }
  }
  structure_rank <- as.integer(structure_rank)[1]
  if (!is.finite(structure_rank) || structure_rank < 0) {
    stop("`structure_rank` must be a nonnegative integer.", call. = FALSE)
  }

  solver <- match.arg(solver)
  sinkhorn_method <- match.arg(sinkhorn_method)
  precision <- match.arg(precision)
  use_ppa <- identical(solver, "PPA")
  use_log_sinkhorn <- identical(sinkhorn_method, "log")
  use_mixed_precision <- identical(precision, "mixed")

  cpp_fgw_entropic_square(
    M = unname(M),
    C1 = unname(C1),
    C2 = unname(C2),
    p = unname(p),
    q = unname(q),
    alpha = alpha,
    epsilon = epsilon,
    max_iter = as.integer(max_iter),
    tol = tol,
    sinkhorn_max_iter = as.integer(sinkhorn_max_iter),
    sinkhorn_tol = sinkhorn_tol,
    init_plan = unname(init_plan),
    symmetric = isTRUE(symmetric),
    use_ppa = use_ppa,
    use_log_sinkhorn = use_log_sinkhorn,
    use_mixed_precision = use_mixed_precision,
    check_every = as.integer(check_every),
    approx_rank = structure_rank
  )
}

#' Entropic Fused Gromov-Wasserstein Objective Value
#'
#' @inheritParams fgw_entropic
#' @return Numeric scalar FGW value.
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
#' @param wx Source weights (default uniform).
#' @param wy Target weights (default uniform).
#' @param reg_marginals Length-1 or length-2 marginal relaxation parameters.
#' @param epsilon Joint KL regularization parameter (entropic term).
#' @param alpha Linear FGW term coefficient.
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
#'   `iterations`, and `error`.
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
  .assert_matrix(Cx, "Cx")
  .assert_matrix(Cy, "Cy")
  nx <- nrow(Cx)
  ny <- nrow(Cy)
  if (ncol(Cx) != nx) stop("`Cx` must be square.", call. = FALSE)
  if (ncol(Cy) != ny) stop("`Cy` must be square.", call. = FALSE)

  if (is.null(wx)) wx <- rep(1 / nx, nx)
  if (is.null(wy)) wy <- rep(1 / ny, ny)
  wx <- .assert_prob(wx, nx, "wx")
  wy <- .assert_prob(wy, ny, "wy")

  if (is.null(M)) {
    M <- matrix(0, nrow = nx, ncol = ny)
  } else {
    .assert_matrix(M, "M")
    if (nrow(M) != nx || ncol(M) != ny) {
      stop("`M` must have shape nrow(Cx) x nrow(Cy).", call. = FALSE)
    }
  }

  if (length(reg_marginals) == 1L) {
    reg_marginals <- c(reg_marginals, reg_marginals)
  }
  if (length(reg_marginals) != 2L) {
    stop("`reg_marginals` must have length 1 or 2.", call. = FALSE)
  }

  if (is.null(init_pi)) {
    init_pi <- tcrossprod(wx, wy)
  } else {
    .assert_matrix(init_pi, "init_pi")
    if (nrow(init_pi) != nx || ncol(init_pi) != ny) {
      stop("`init_pi` must have shape nrow(Cx) x nrow(Cy).", call. = FALSE)
    }
  }

  cpp_fugw_kl_square(
    Cx = unname(Cx),
    Cy = unname(Cy),
    wx = unname(wx),
    wy = unname(wy),
    reg_marginals = as.numeric(reg_marginals),
    epsilon = epsilon,
    alpha = alpha,
    M = unname(M),
    init_pi = unname(init_pi),
    max_iter = as.integer(max_iter),
    tol = tol,
    max_iter_ot = as.integer(max_iter_ot),
    tol_ot = tol_ot,
    rescale_plan = isTRUE(rescale_plan),
    check_every = as.integer(check_every),
    use_mixed_precision = identical(precision, "mixed")
  )
}

#' Fused Unbalanced Gromov-Wasserstein Objective Value
#'
#' @inheritParams fugw_kl
#' @return Numeric scalar FUGW objective value.
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
      "`fgw_exact_cg()` requires package `lpSolve` for the EMD/LP direction step.",
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

#' Exact FGW via Conditional Gradient + LP Direction (square loss)
#'
#' POT-style CG iterations with a linear-program transport direction step.
#'
#' @inheritParams fgw_entropic
#' @param lp_scale Integer mass scaling used for LP marginals.
#' @param lp_solver LP backend for direction step: `"cpp_transport"` (default,
#'   C++ transport-simplex backend), `"lp_transport"` (lpSolve), or `"lp_matrix"`
#'   (cached dense lpSolve constraints).
#' @param lp_max_iter Maximum iterations for the C++ transport-simplex LP backend.
#' @param lp_tol Optimality tolerance for the C++ transport-simplex LP backend.
#' @return A list with `plan`, `fgw_dist`, `iterations`, `error`, `loss_trace`.
#' @export
fgw_exact_cg <- function(
    M,
    C1,
    C2,
    p = NULL,
    q = NULL,
    alpha = 0.5,
    symmetric = TRUE,
    G0 = NULL,
    max_iter = 500L,
    tol_rel = 1e-9,
    tol_abs = 1e-9,
    lp_scale = 1e6,
    lp_solver = c("cpp_transport", "lp_transport", "lp_matrix"),
    lp_max_iter = 20000L,
    lp_tol = 1e-12) {
  .assert_matrix(M, "M")
  .assert_matrix(C1, "C1")
  .assert_matrix(C2, "C2")

  ns <- nrow(C1)
  nt <- nrow(C2)
  if (ncol(C1) != ns) stop("`C1` must be square.", call. = FALSE)
  if (ncol(C2) != nt) stop("`C2` must be square.", call. = FALSE)
  if (nrow(M) != ns || ncol(M) != nt) {
    stop("`M` must have shape nrow(C1) x nrow(C2).", call. = FALSE)
  }

  if (is.null(p)) p <- rep(1 / ns, ns)
  if (is.null(q)) q <- rep(1 / nt, nt)
  p <- .assert_prob(p, ns, "p")
  q <- .assert_prob(q, nt, "q")
  lp_solver <- match.arg(lp_solver)

  if (identical(lp_solver, "cpp_transport")) {
    return(cpp_fgw_exact_cg_square(
      M = unname(M),
      C1 = unname(C1),
      C2 = unname(C2),
      p = unname(p),
      q = unname(q),
      alpha = alpha,
      symmetric = isTRUE(symmetric),
      max_iter = as.integer(max_iter),
      tol_rel = tol_rel,
      tol_abs = tol_abs,
      lp_max_iter = as.integer(lp_max_iter),
      lp_tol = lp_tol
    ))
  }

  if (is.null(G0)) {
    G <- tcrossprod(p, q)
  } else {
    .assert_matrix(G0, "G0")
    if (nrow(G0) != ns || ncol(G0) != nt) {
      stop("`G0` must have shape nrow(C1) x nrow(C2).", call. = FALSE)
    }
    if (max(abs(rowSums(G0) - p)) > 1e-8 || max(abs(colSums(G0) - q)) > 1e-8) {
      stop("`G0` must satisfy the marginal constraints.", call. = FALSE)
    }
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

  list(
    plan = G,
    fgw_dist = cost_G,
    iterations = it,
    error = abs_delta,
    rel_error = rel_delta,
    loss_trace = loss_trace[seq_len(it + 1L)]
  )
}
