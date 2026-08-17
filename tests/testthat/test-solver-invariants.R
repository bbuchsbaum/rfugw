# Invariant and metamorphic suite for flagship solvers.
# Seed: 20260816. Tolerances are absolute-plus-relative: 1e-6 abs on small
# objectives, 1e-5 on plans, 1e-4 on adversarial / mixed-precision paths.

.inv_data <- function(seed = 20260816, ns = 5, nt = 6) {
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

.expect_balanced_plan <- function(out, p, q, tol = 1e-5) {
  plan <- out$plan
  expect_true(is.matrix(plan))
  expect_true(all(is.finite(plan)))
  expect_true(all(plan >= -1e-12))
  expect_equal(rowSums(plan), p, tolerance = tol)
  expect_equal(colSums(plan), q, tolerance = tol)
  expect_true(is.finite(out$fgw_dist %||% out$gw_dist %||% out$fugw_cost))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

test_that("flagship solvers keep nonnegative finite plans and marginals", {
  d <- .inv_data()
  ent <- fgw_entropic(d$M, d$C1, d$C2, p = d$p, q = d$q, epsilon = 0.08, max_iter = 80L)
  .expect_balanced_plan(ent, d$p, d$q)
  exact <- fgw_exact_cg(d$M, d$C1, d$C2, p = d$p, q = d$q, max_iter = 40L)
  .expect_balanced_plan(exact, d$p, d$q)
  gw <- gromov_wasserstein(d$C1, d$C2, p = d$p, q = d$q, max_iter = 40L)
  .expect_balanced_plan(list(plan = gw$plan, fgw_dist = gw$gw_dist), d$p, d$q)
  fugw <- fugw_kl(d$C1, d$C2, M = d$M, epsilon = 0.02, max_iter = 30L)
  expect_true(all(is.finite(fugw$pi_samp)) && all(fugw$pi_samp >= -1e-12))
  expect_true(is.finite(fugw$fugw_cost))
})

test_that("precision, Sinkhorn, and solver paths stay consistent", {
  d <- .inv_data(seed = 7, ns = 5, nt = 5)
  args <- list(M = d$M, C1 = d$C1, C2 = d$C2, p = d$p, q = d$q, epsilon = 0.1, max_iter = 60L)
  mixed <- do.call(fgw_entropic, c(args, list(precision = "mixed")))
  dbl <- do.call(fgw_entropic, c(args, list(precision = "double")))
  logp <- do.call(fgw_entropic, c(args, list(sinkhorn_method = "log", precision = "double")))
  ppa <- do.call(fgw_entropic, c(args, list(solver = "PPA", precision = "double")))
  expect_equal(mixed$fgw_dist, dbl$fgw_dist, tolerance = 1e-4)
  expect_equal(logp$fgw_dist, dbl$fgw_dist, tolerance = 5e-4)
  expect_true(is.finite(ppa$fgw_dist))
  .expect_balanced_plan(mixed, d$p, d$q, tol = 5e-5)
  .expect_balanced_plan(logp, d$p, d$q, tol = 5e-5)
})

test_that("source permutation is a metamorphic invariant", {
  d <- .inv_data(seed = 11, ns = 5, nt = 5)
  perm <- c(3L, 1L, 5L, 2L, 4L)
  ent <- fgw_entropic(d$M, d$C1, d$C2, p = d$p, q = d$q, epsilon = 0.1, max_iter = 80L)
  ent_p <- fgw_entropic(
    d$M[perm, ],
    d$C1[perm, perm],
    d$C2,
    p = d$p[perm],
    q = d$q,
    epsilon = 0.1,
    max_iter = 80L
  )
  expect_equal(ent_p$fgw_dist, ent$fgw_dist, tolerance = 1e-5)
  expect_equal(ent_p$plan, ent$plan[perm, , drop = FALSE], tolerance = 1e-4)
})

test_that("structure scaling scales the unregularized GW term by the square", {
  d <- .inv_data(seed = 13, ns = 4, nt = 4)
  scale <- 2.5
  gw <- gromov_wasserstein(d$C1, d$C2, p = d$p, q = d$q, max_iter = 40L)
  gw_s <- gromov_wasserstein(scale * d$C1, scale * d$C2, p = d$p, q = d$q, max_iter = 40L)
  expect_equal(gw_s$gw_dist, (scale^2) * gw$gw_dist, tolerance = 1e-5)
})

test_that("common additive offset on both structure costs is invariant", {
  d <- .inv_data(seed = 17, ns = 4, nt = 4)
  off <- 3
  gw <- gromov_wasserstein(d$C1, d$C2, p = d$p, q = d$q, max_iter = 40L)
  gw_o <- gromov_wasserstein(d$C1 + off, d$C2 + off, p = d$p, q = d$q, max_iter = 40L)
  expect_equal(gw_o$gw_dist, gw$gw_dist, tolerance = 1e-6)
  expect_equal(gw_o$plan, gw$plan, tolerance = 1e-5)
})

test_that("role swap stays feasible and finite", {
  d <- .inv_data(seed = 19, ns = 5, nt = 5)
  fwd <- fgw_entropic(d$M, d$C1, d$C2, p = d$p, q = d$q, epsilon = 0.1, max_iter = 60L)
  rev <- fgw_entropic(t(d$M), d$C2, d$C1, p = d$q, q = d$p, epsilon = 0.1, max_iter = 60L)
  .expect_balanced_plan(rev, d$q, d$p)
  expect_equal(rev$fgw_dist, fwd$fgw_dist, tolerance = 5e-4)
})

test_that("adversarial weights, duplicates, and tiny epsilon stay finite", {
  d <- .inv_data(seed = 23, ns = 5, nt = 5)
  p_conc <- c(0.96, 0.01, 0.01, 0.01, 0.01)
  q_conc <- c(0.01, 0.01, 0.01, 0.01, 0.96)
  tiny <- fgw_entropic(
    d$M, d$C1, d$C2,
    p = p_conc,
    q = q_conc,
    epsilon = 1e-4,
    max_iter = 40L,
    sinkhorn_method = "log",
    precision = "double"
  )
  .expect_balanced_plan(tiny, p_conc / sum(p_conc), q_conc / sum(q_conc), tol = 5e-4)

  C1_dup <- d$C1
  C1_dup[2, ] <- C1_dup[1, ]
  C1_dup[, 2] <- C1_dup[, 1]
  C1_dup[2, 2] <- 0
  dup <- fgw_entropic(d$M, C1_dup, d$C2, epsilon = 0.1, max_iter = 40L)
  expect_true(is.finite(dup$fgw_dist))
  expect_true(all(is.finite(dup$plan)))

  ill <- fgw_entropic(1e6 * d$M, 1e3 * d$C1, 1e3 * d$C2, epsilon = 5, max_iter = 20L)
  expect_true(is.finite(ill$fgw_dist))

  fugw <- fugw_kl(
    d$C1, d$C2,
    M = d$M,
    epsilon = 1e-3,
    reg_marginals = c(1e-2, 1e3),
    max_iter = 20L
  )
  expect_true(is.finite(fugw$fugw_cost))
  expect_true(all(is.finite(fugw$pi_samp)))
})

test_that("semirelaxed and sampled experimental paths stay finite", {
  d <- .inv_data(seed = 29, ns = 5, nt = 6)
  sr <- entropic_semirelaxed_gromov_wasserstein(d$C1, d$C2, p = d$p, epsilon = 0.1, max_iter = 40L)
  expect_true(is.finite(sr$srgw_dist))
  expect_equal(rowSums(sr$plan), d$p, tolerance = 1e-5)
  expect_true(all(sr$plan >= -1e-12))

  samp <- sampled_gromov_wasserstein(
    d$C1, d$C2,
    p = d$p,
    q = d$q,
    epsilon = 0.1,
    nb_samples_grad = 8L,
    max_iter = 20L,
    log = TRUE
  )
  expect_true(is.finite(samp$gw_dist_estimated) || is.matrix(samp$plan))
})
