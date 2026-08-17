make_pair <- function(seed = 1, ns = 6, nt = 7) {
  set.seed(seed)
  C1 <- as.matrix(dist(matrix(rnorm(ns * 2), ns, 2)))
  C2 <- as.matrix(dist(matrix(rnorm(nt * 2), nt, 2)))
  M <- as.matrix(dist(rbind(
    matrix(rnorm(ns * 2), ns, 2),
    matrix(rnorm(nt * 2), nt, 2)
  )))[seq_len(ns), ns + seq_len(nt)]
  list(
    C1 = C1 / max(C1),
    C2 = C2 / max(C2),
    M = M / max(M),
    p = rep(1 / ns, ns),
    q = rep(1 / nt, nt)
  )
}

test_that("package metadata is real and citation validates", {
  desc <- read.dcf(system.file("DESCRIPTION", package = "rfugw"))
  expect_false(grepl("Codex|example.com", desc[, "Authors@R"]))
  expect_equal(unname(desc[, "URL"]), "https://github.com/bbuchsbaum/rfugw")
  expect_equal(unname(desc[, "BugReports"]), "https://github.com/bbuchsbaum/rfugw/issues")
  cit <- citation("rfugw")
  expect_s3_class(cit, "citation")
  expect_true(any(grepl("Buchsbaum", format(cit, "text"))))
})

test_that("UCOOT defaults succeed and unsupported choices fail early", {
  set.seed(3)
  X <- matrix(rnorm(12), 4, 3)
  Y <- matrix(rnorm(15), 5, 3)

  out <- unbalanced_co_optimal_transport(X, Y)
  expect_true(is.matrix(out$pi_samp))
  expect_true(is.matrix(out$pi_feat))
  expect_true(out$status %in% c("converged", "max_iter"))
  expect_type(out$converged, "logical")

  out2 <- fused_unbalanced_across_spaces_divergence(X, Y)
  expect_true(is.matrix(out2$pi_samp))

  expect_error(
    unbalanced_co_optimal_transport(X, Y, divergence = "l2"),
    "kl"
  )
  expect_error(
    unbalanced_co_optimal_transport(X, Y, unbalanced_solver = "mm"),
    "sinkhorn"
  )
  expect_error(
    unbalanced_co_optimal_transport(X, Y, unbalanced_solver = "lbfgsb"),
    "sinkhorn"
  )
  expect_error(
    unbalanced_co_optimal_transport(X, Y, epsilon = 0),
    "positive"
  )
  expect_error(
    unbalanced_co_optimal_transport(X, Y, not_a_real_arg = 1),
    "Unsupported argument"
  )
})

test_that("exact CG honors G0 on cpp and lpSolve backends", {
  d <- make_pair(11, ns = 5, nt = 5)
  G_prod <- tcrossprod(d$p, d$q)
  G_perm <- diag(d$p)

  eval_start <- function(G0, solver) {
    fgw_exact_cg(
      d$M, d$C1, d$C2,
      p = d$p,
      q = d$q,
      alpha = 0.4,
      G0 = G0,
      max_iter = 1L,
      lp_solver = solver
    )
  }

  cpp_prod <- eval_start(G_prod, "cpp_transport")
  cpp_perm <- eval_start(G_perm, "cpp_transport")
  cpp_default <- fgw_exact_cg(
    d$M, d$C1, d$C2,
    p = d$p,
    q = d$q,
    alpha = 0.4,
    max_iter = 1L,
    lp_solver = "cpp_transport"
  )
  expect_gt(abs(cpp_prod$loss_trace[[1]] - cpp_perm$loss_trace[[1]]), 1e-10)
  expect_equal(cpp_default$loss_trace[[1]], cpp_prod$loss_trace[[1]], tolerance = 1e-10)

  skip_if_not_installed("lpSolve")
  lp_prod <- eval_start(G_prod, "lp_transport")
  lp_perm <- eval_start(G_perm, "lp_transport")
  expect_equal(cpp_prod$loss_trace[[1]], lp_prod$loss_trace[[1]], tolerance = 1e-8)
  expect_equal(cpp_perm$loss_trace[[1]], lp_perm$loss_trace[[1]], tolerance = 1e-8)
  expect_equal(cpp_prod$fgw_dist, lp_prod$fgw_dist, tolerance = 1e-6)
  expect_equal(cpp_perm$fgw_dist, lp_perm$fgw_dist, tolerance = 1e-6)

  expect_error(
    fgw_exact_cg(d$M, d$C1, d$C2, G0 = matrix(-1, 5, 5)),
    "nonnegative"
  )
  expect_error(
    fgw_exact_cg(d$M, d$C1, d$C2, G0 = matrix(1, 4, 5)),
    "shape"
  )
  bad <- G_prod
  bad[1, 1] <- bad[1, 1] + 0.2
  expect_error(fgw_exact_cg(d$M, d$C1, d$C2, G0 = bad), "marginal")
})

test_that("symmetry is auto-detected and symmetric=TRUE is validated", {
  d <- make_pair(4)
  C1_asymm <- d$C1
  C1_asymm[1, 2] <- C1_asymm[1, 2] + 0.3

  expect_error(
    fgw_entropic(d$M, C1_asymm, d$C2, symmetric = TRUE, max_iter = 2L),
    "symmetric"
  )
  expect_error(
    fgw_exact_cg(d$M, C1_asymm, d$C2, symmetric = TRUE, max_iter = 2L),
    "symmetric"
  )

  out_auto <- fgw_entropic(d$M, C1_asymm, d$C2, epsilon = 0.1, max_iter = 40L)
  out_false <- fgw_entropic(
    d$M, C1_asymm, d$C2,
    epsilon = 0.1,
    max_iter = 40L,
    symmetric = FALSE
  )
  expect_equal(out_auto$fgw_dist, out_false$fgw_dist, tolerance = 1e-8)
  expect_equal(out_auto$plan, out_false$plan, tolerance = 1e-6)

  out_fast <- fgw_entropic(d$M, d$C1, d$C2, epsilon = 0.1, max_iter = 40L, symmetric = TRUE)
  out_gen <- fgw_entropic(d$M, d$C1, d$C2, epsilon = 0.1, max_iter = 40L, symmetric = FALSE)
  expect_equal(out_fast$fgw_dist, out_gen$fgw_dist, tolerance = 1e-6)
})

test_that("flagship solvers expose status and distinguish max_iter from convergence", {
  d <- make_pair(8)
  ok <- fgw_entropic(d$M, d$C1, d$C2, epsilon = 0.05, max_iter = 200L, tol = 1e-6)
  expect_true(ok$status %in% c("converged", "max_iter"))
  expect_true(is.logical(ok$converged))
  expect_true(is.finite(ok$residual))
  expect_true(is.finite(ok$row_residual))
  expect_true(is.finite(ok$col_residual))

  limited <- fgw_entropic(d$M, d$C1, d$C2, epsilon = 0.05, max_iter = 1L, tol = 1e-15)
  expect_equal(limited$status, "max_iter")
  expect_false(limited$converged)

  exact <- fgw_exact_cg(d$M, d$C1, d$C2, max_iter = 80L)
  expect_true(exact$status %in% c("converged", "max_iter", "lp_failure"))
  expect_true(is.finite(exact$residual))

  fugw <- fugw_kl(d$C1, d$C2, M = d$M, max_iter = 1L, tol = 1e-15)
  expect_equal(fugw$status, "max_iter")
  expect_false(fugw$converged)
})

test_that("parameter validation rejects invalid alpha, epsilon, and empty mass", {
  d <- make_pair(2)
  expect_error(fgw_entropic(d$M, d$C1, d$C2, alpha = 1.5), "alpha")
  expect_error(fgw_entropic(d$M, d$C1, d$C2, epsilon = 0), "epsilon")
  expect_error(fgw_entropic(d$M, d$C1, d$C2, max_iter = 0), "max_iter")
  expect_error(fgw_entropic(d$M, d$C1, d$C2, p = rep(0, nrow(d$C1))), "positive total mass")
  expect_error(fgw_exact_cg(d$M, d$C1, d$C2, alpha = -0.1), "alpha")
  expect_error(fugw_kl(d$C1, d$C2, epsilon = -1), "epsilon")
  expect_error(fgw_entropic(d$M + NA, d$C1, d$C2), "finite")
  expect_error(
    fgw_entropic(matrix(0, 0, 0), matrix(0, 0, 0), matrix(0, 0, 0)),
    "nonempty"
  )
})

test_that("singleton problems stay finite and weights are renormalized", {
  C1 <- matrix(0, 1, 1)
  C2 <- matrix(0, 1, 1)
  M <- matrix(2, 1, 1)
  out <- fgw_entropic(M, C1, C2, p = 4, q = 4, epsilon = 0.1, max_iter = 20L)
  expect_true(is.finite(out$fgw_dist))
  expect_equal(sum(out$plan), 1, tolerance = 1e-8)
  exact <- fgw_exact_cg(M, C1, C2, p = 3, q = 3, max_iter = 5L)
  expect_true(is.finite(exact$fgw_dist))
  expect_equal(sum(exact$plan), 1, tolerance = 1e-8)
})

test_that("numerical breakdown cannot be reported as success", {
  broken <- list(plan = matrix(NaN, 2, 2), fgw_dist = NaN, iterations = 3L)
  out <- rfugw:::.attach_solver_diagnostics(
    broken,
    residual = 0,
    converged = TRUE,
    iterations = 3L,
    max_iter = 10L,
    plan = broken$plan
  )
  expect_equal(out$status, "numerical_failure")
  expect_false(out$converged)
})

independent_fgw <- function(M, C1, C2, G, alpha) {
  ns <- nrow(C1)
  nt <- nrow(C2)
  quad <- 0
  for (i in seq_len(ns)) {
    for (j in seq_len(ns)) {
      for (k in seq_len(nt)) {
        for (l in seq_len(nt)) {
          quad <- quad + (C1[i, j] - C2[k, l])^2 * G[i, k] * G[j, l]
        }
      }
    }
  }
  (1 - alpha) * sum(M * G) + alpha * quad
}

test_that("reported unregularized objectives match an independent recomputation", {
  d <- make_pair(5, ns = 4, nt = 4)
  alpha <- 0.35
  exact <- fgw_exact_cg(d$M, d$C1, d$C2, p = d$p, q = d$q, alpha = alpha, max_iter = 40L)
  expect_equal(
    exact$fgw_dist,
    independent_fgw(d$M, d$C1, d$C2, exact$plan, alpha),
    tolerance = 1e-7
  )
  expect_equal(fused_gromov_wasserstein2(d$M, d$C1, d$C2, alpha = alpha, max_iter = 40L),
    exact$fgw_dist,
    tolerance = 1e-8
  )
  ent <- fgw_entropic(d$M, d$C1, d$C2, p = d$p, q = d$q, alpha = alpha, epsilon = 0.08, max_iter = 80L)
  expect_equal(
    ent$fgw_dist,
    independent_fgw(d$M, d$C1, d$C2, ent$plan, alpha),
    tolerance = 1e-6
  )
})
