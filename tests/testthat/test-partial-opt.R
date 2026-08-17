make_partial_problem <- function(ns, nt, seed) {
  set.seed(seed)
  Xs <- matrix(stats::rnorm(ns * 3L), ns, 3L)
  Xt <- matrix(stats::rnorm(nt * 3L), nt, 3L)
  Fs <- matrix(stats::rnorm(ns * 2L), ns, 2L)
  Ft <- matrix(stats::rnorm(nt * 2L), nt, 2L)
  C1 <- as.matrix(stats::dist(Xs)); C1 <- C1 / max(C1)
  C2 <- as.matrix(stats::dist(Xt)); C2 <- C2 / max(C2)
  M <- as.matrix(stats::dist(rbind(Fs, Ft)))[seq_len(ns), ns + seq_len(nt)]
  M <- M / max(M)
  list(C1 = C1, C2 = C2, M = M, p = rep(1 / ns, ns), q = rep(1 / nt, nt))
}

test_that("native exact partial matches the lpSolve CG core on mass and objective", {
  skip_if_not_installed("lpSolve")
  d <- make_partial_problem(6L, 5L, seed = 20260816L)
  m <- 0.6
  cpp <- rfugw::partial_gromov_wasserstein(
    d$C1, d$C2, p = d$p, q = d$q, m = m,
    numItermax = 40L, tol = 1e-8, log = TRUE, lp_solver = "cpp_transport"
  )
  lp <- rfugw::partial_gromov_wasserstein(
    d$C1, d$C2, p = d$p, q = d$q, m = m,
    numItermax = 40L, tol = 1e-8, log = TRUE, lp_solver = "lp_matrix"
  )
  expect_equal(sum(cpp$plan), m, tolerance = 1e-5)
  expect_equal(sum(lp$plan), m, tolerance = 1e-5)
  expect_true(all(rowSums(cpp$plan) <= d$p + 5e-4))
  expect_true(all(colSums(cpp$plan) <= d$q + 5e-4))
  expect_equal(cpp$partial_gw_dist, lp$partial_gw_dist, tolerance = 5e-3)
})

test_that("C++ exact partial FGW matches the R CG core", {
  skip_if_not_installed("lpSolve")
  d <- make_partial_problem(5L, 5L, seed = 9L)
  m <- 0.55
  G0 <- (d$p %o% d$q) * (m / (sum(d$p) * sum(d$q)))
  cpp <- rfugw:::cpp_partial_fgw_exact_square(
    d$M, d$C1, d$C2, d$p, d$q, m, 0.5, TRUE,
    G0, 30L, 1e-8, 1L, 20000L, 1e-12
  )
  ref <- rfugw:::.partial_fgw_cg_core(
    d$M, d$C1, d$C2, d$p, d$q, m, 0.5, G0, 30L, 1e-8, TRUE,
    "lp_matrix", 1e6, 1L
  )
  expect_equal(cpp$objective, ref$objective, tolerance = 5e-3)
  expect_equal(sum(cpp$plan), m, tolerance = 1e-5)
})

test_that("C++ entropic partial FGW matches the R outer loop", {
  d <- make_partial_problem(6L, 5L, seed = 4L)
  m <- 0.7
  G0 <- (d$p %o% d$q) * (m / (sum(d$p) * sum(d$q)))
  cpp <- rfugw:::cpp_partial_fgw_entropic_square(
    d$M, d$C1, d$C2, d$p, d$q, m, 0.2, 0.6, TRUE,
    G0, 40L, 1e-7, 200L, 1e-12, 2L
  )
  ref <- rfugw:::.partial_fgw_entropic_core(
    d$M, d$C1, d$C2, d$p, d$q, m, 0.2, 0.6, G0, 40L, 1e-7, TRUE,
    200L, 1e-12, FALSE, 2L
  )
  expect_equal(cpp$objective, ref$objective, tolerance = 1e-8)
  expect_equal(unname(cpp$plan), unname(ref$plan), tolerance = 1e-8)
  expect_equal(sum(cpp$plan), m, tolerance = 1e-4)
})

test_that("empty M at alpha=1 matches an explicit zero feature cost", {
  d <- make_partial_problem(5L, 4L, seed = 2L)
  z <- matrix(0, nrow(d$C1), nrow(d$C2))
  empty <- rfugw:::cpp_partial_fgw_entropic_square(
    matrix(numeric(0), 0, 0), d$C1, d$C2, d$p, d$q, 0.6, 0.25, 1, TRUE,
    matrix(numeric(0), 0, 0), 20L, 1e-7, 100L, 1e-12, 2L
  )
  zero <- rfugw:::cpp_partial_fgw_entropic_square(
    z, d$C1, d$C2, d$p, d$q, 0.6, 0.25, 1, TRUE,
    matrix(numeric(0), 0, 0), 20L, 1e-7, 100L, 1e-12, 2L
  )
  expect_equal(empty$objective, zero$objective, tolerance = 1e-12)
  expect_equal(unname(empty$plan), unname(zero$plan), tolerance = 1e-12)
})

test_that("POT multi-size entropic partial objectives stay in range", {
  fx <- read_fixture("partial_pot_multisize.json")
  cases <- fx$cases
  if (is.data.frame(cases)) {
    # jsonlite may simplify a list of objects to a data.frame
    n_cases <- nrow(cases)
    get_case <- function(i) {
      list(
        C1 = cases$C1[[i]], C2 = cases$C2[[i]], M = cases$M[[i]],
        p = cases$p[[i]], q = cases$q[[i]],
        m = cases$m[[i]], reg = cases$reg[[i]], alpha = cases$alpha[[i]],
        numItermax = cases$numItermax[[i]], tol = cases$tol[[i]],
        epgw_dist = cases$epgw_dist[[i]], epfgw_dist = cases$epfgw_dist[[i]]
      )
    }
  } else {
    n_cases <- length(cases)
    get_case <- function(i) cases[[i]]
  }
  expect_gte(n_cases, 4L)
  for (i in seq_len(n_cases)) {
    cs <- get_case(i)
    gw <- rfugw::entropic_partial_gromov_wasserstein(
      C1 = cs$C1, C2 = cs$C2, p = cs$p, q = cs$q, m = cs$m, reg = cs$reg,
      numItermax = as.integer(cs$numItermax), tol = cs$tol, log = TRUE
    )
    fgw <- rfugw::entropic_partial_fused_gromov_wasserstein(
      M = cs$M, C1 = cs$C1, C2 = cs$C2, p = cs$p, q = cs$q, m = cs$m,
      reg = cs$reg, alpha = cs$alpha,
      numItermax = as.integer(cs$numItermax), tol = cs$tol, log = TRUE
    )
    # GW matches POT tightly. POT's entropic partial FGW gradient adds the
    # scalar (1-alpha)*sum(G*M) instead of (1-alpha)*M, so FGW objectives
    # are checked against the square-loss definition, not POT's FGW log.
    expect_lt(abs(gw$partial_gw_dist - cs$epgw_dist), 1e-4)
    expect_equal(
      fgw$partial_fgw_dist,
      rfugw::ot_fgw_square(cs$M, cs$C1, cs$C2, fgw, alpha = cs$alpha),
      tolerance = 1e-8
    )
    expect_equal(sum(gw$plan), cs$m, tolerance = 1e-4)
    expect_equal(sum(fgw$plan), cs$m, tolerance = 1e-4)
    expect_true(all(rowSums(gw$plan) <= cs$p + 1e-6))
    expect_true(all(colSums(gw$plan) <= cs$q + 1e-6))
    expect_true(all(rowSums(fgw$plan) <= cs$p + 1e-6))
    expect_true(all(colSums(fgw$plan) <= cs$q + 1e-6))
  }
})

test_that("supplied partial G0 is honored", {
  d <- make_partial_problem(4L, 4L, seed = 1L)
  m <- 0.5
  G0 <- diag(4) * (m / 4)
  out <- rfugw::entropic_partial_gromov_wasserstein(
    d$C1, d$C2, p = d$p, q = d$q, m = m, reg = 0.3, G0 = G0,
    numItermax = 1L, tol = 0, log = TRUE
  )
  expect_equal(sum(out$plan), m, tolerance = 1e-4)
  expect_false(isTRUE(all.equal(
    unname(out$plan),
    unname((d$p %o% d$q) * m),
    tolerance = 1e-6
  )))
})
