.symmetry_tol <- 1e-10

.resolve_symmetric <- function(symmetric, C1, C2, tol = .symmetry_tol) {
  detected <- .is_symmetric_cost(C1, tol) && .is_symmetric_cost(C2, tol)
  if (is.null(symmetric)) {
    return(detected)
  }
  if (length(symmetric) != 1L || is.na(symmetric)) {
    stop("`symmetric` must be TRUE, FALSE, or NULL.", call. = FALSE)
  }
  symmetric <- isTRUE(symmetric)
  if (symmetric && !detected) {
    stop(
      "`symmetric = TRUE` requires C1 and C2 to be symmetric within 1e-10.",
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
  x <- suppressWarnings(as.integer(x)[1])
  if (!is.finite(x) || is.na(x) || x < min) {
    stop(sprintf("`%s` must be an integer >= %d.", name, min), call. = FALSE)
  }
  x
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

.attach_solver_diagnostics <- function(out,
                                       residual,
                                       converged,
                                       iterations,
                                       max_iter,
                                       p = NULL,
                                       q = NULL,
                                       plan = NULL,
                                       inner_residual = NA_real_,
                                       inner_iterations = NA_integer_,
                                       numerical_ok = TRUE,
                                       lp_ok = TRUE) {
  if (is.null(plan)) {
    plan <- out$plan
    if (is.null(plan) && !is.null(out$pi_samp)) {
      plan <- out$pi_samp
    }
  }
  residuals <- .plan_marginal_residuals(plan, p, q)
  plan_ok <- is.null(plan) || (is.matrix(plan) && all(is.finite(plan)))
  value_fields <- c("fgw_dist", "gw_dist", "fugw_cost", "ucoot_cost",
                    "srgw_dist", "srfgw_dist", "partial_gw_dist", "partial_fgw_dist")
  value_ok <- TRUE
  for (nm in value_fields) {
    if (!is.null(out[[nm]]) && !all(is.finite(out[[nm]]))) {
      value_ok <- FALSE
    }
  }
  if (!isTRUE(numerical_ok) || !isTRUE(plan_ok) || !isTRUE(value_ok)) {
    status <- "numerical_failure"
  } else if (!isTRUE(lp_ok)) {
    status <- "lp_failure"
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
  out$inner_residual <- as.numeric(inner_residual)
  out$inner_iterations <- as.integer(inner_iterations)
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
  ns <- as.integer(ns)[1]
  nt <- as.integer(nt)[1]
  if (!is.finite(ns) || is.na(ns) || ns < 1L ||
      !is.finite(nt) || is.na(nt) || nt < 1L) {
    stop("`ns` and `nt` must be integers >= 1.", call. = FALSE)
  }

  if (length(nb_samples_grad) == 1L) {
    nb <- suppressWarnings(as.integer(nb_samples_grad)[1])
    if (!is.finite(nb) || is.na(nb) || nb < 1L) {
      stop("`nb_samples_grad` must be >= 1.", call. = FALSE)
    }
    if (nb > ns) {
      nb_p <- ns
      nb_q <- max(1L, as.integer(nb %/% ns))
    } else {
      nb_p <- nb
      nb_q <- 1L
    }
  } else {
    nb <- suppressWarnings(as.integer(nb_samples_grad))
    if (length(nb) != 2L || any(!is.finite(nb)) || any(is.na(nb)) || any(nb < 1L)) {
      stop("`nb_samples_grad` must be an integer >=1 or a length-2 integer vector.", call. = FALSE)
    }
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
  rmax <- min(as.integer(ns)[1], as.integer(nt)[1])
  if (!is.finite(rmax) || is.na(rmax) || rmax < 1L) {
    stop("`ns` and `nt` must be integers >= 1.", call. = FALSE)
  }
  if (is.null(rank)) {
    return(as.integer(rmax))
  }
  r <- suppressWarnings(as.integer(rank)[1])
  if (!is.finite(r) || is.na(r) || r < 1L) {
    stop("`rank` must be an integer >= 1.", call. = FALSE)
  }
  if (r > rmax) {
    warning(
      sprintf("`rank` = %d exceeds min(ns, nt) = %d; clamping.", r, rmax),
      call. = FALSE
    )
    r <- rmax
  }
  as.integer(r)
}
