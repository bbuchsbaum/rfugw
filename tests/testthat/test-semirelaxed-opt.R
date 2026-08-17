make_sr_problem <- function(ns, nt, seed, asymmetric = FALSE) {
  set.seed(seed)
  Xs <- matrix(stats::rnorm(ns * 3L), ns, 3L)
  Xt <- matrix(stats::rnorm(nt * 3L), nt, 3L)
  Fs <- matrix(stats::rnorm(ns * 2L), ns, 2L)
  Ft <- matrix(stats::rnorm(nt * 2L), nt, 2L)
  C1 <- as.matrix(stats::dist(Xs)); C1 <- C1 / max(C1)
  C2 <- as.matrix(stats::dist(Xt)); C2 <- C2 / max(C2)
  if (isTRUE(asymmetric)) {
    C1[1L, 2L] <- C1[1L, 2L] + 0.2
    C2[2L, 1L] <- C2[2L, 1L] + 0.15
  }
  M <- as.matrix(stats::dist(rbind(Fs, Ft)))[seq_len(ns), ns + seq_len(nt)]
  M <- M / max(M)
  list(C1 = C1, C2 = C2, M = M, p = rep(1 / ns, ns))
}

empty_M <- function() matrix(numeric(0), nrow = 0, ncol = 0)

test_that("asymmetric unregularized SR C++ matches the R CG core", {
  d <- make_sr_problem(8L, 6L, seed = 20260816L, asymmetric = TRUE)
  ctrl <- list(max_iter = 80L, tol_rel = 1e-9, tol_abs = 1e-9)
  cpp <- rfugw:::cpp_semirelaxed_fgw_exact_square(
    M = d$M, C1 = d$C1, C2 = d$C2, p = d$p, alpha = 0.55,
    symmetric = FALSE, init_plan = empty_M(),
    max_iter = ctrl$max_iter, tol_rel = ctrl$tol_rel, tol_abs = ctrl$tol_abs
  )
  ref <- rfugw:::.semirelaxed_fgw_cg_core(
    M = d$M, C1 = d$C1, C2 = d$C2, p = d$p, alpha = 0.55,
    symmetric = FALSE, G0 = NULL,
    max_iter = ctrl$max_iter, tol_rel = ctrl$tol_rel, tol_abs = ctrl$tol_abs
  )
  expect_equal(unname(cpp$plan), unname(ref$plan), tolerance = 1e-8)
  expect_equal(cpp$srfgw_dist, ref$srfgw_dist, tolerance = 1e-8)
  expect_equal(unname(rowSums(cpp$plan)), d$p, tolerance = 1e-10)
})

test_that("public asymmetric SR GW uses the C++ exact path", {
  d <- make_sr_problem(7L, 5L, seed = 11L, asymmetric = TRUE)
  out <- rfugw::semirelaxed_gromov_wasserstein(
    C1 = d$C1, C2 = d$C2, p = d$p, symmetric = FALSE,
    max_iter = 60L, tol_rel = 1e-9, tol_abs = 1e-9
  )
  ref <- rfugw:::.semirelaxed_fgw_cg_core(
    M = matrix(0, nrow(d$C1), nrow(d$C2)),
    C1 = d$C1, C2 = d$C2, p = d$p, alpha = 1, symmetric = FALSE,
    G0 = NULL, max_iter = 60L, tol_rel = 1e-9, tol_abs = 1e-9
  )
  expect_equal(unname(out$plan), unname(ref$plan), tolerance = 1e-8)
  expect_equal(out$srgw_dist, ref$srgw_dist, tolerance = 1e-8)
})

test_that("supplied G0 is honored on exact and fast paths", {
  d <- make_sr_problem(6L, 6L, seed = 3L)
  G0 <- d$p %o% c(0.7, rep(0.3 / 5, 5))
  exact0 <- rfugw:::cpp_semirelaxed_fgw_exact_square(
    M = d$M, C1 = d$C1, C2 = d$C2, p = d$p, alpha = 0.4,
    symmetric = TRUE, init_plan = G0,
    max_iter = 1L, tol_rel = 0, tol_abs = 0
  )
  fast0 <- rfugw:::cpp_semirelaxed_fgw_cg_square_fast(
    M = d$M, C1 = d$C1, C2 = d$C2, p = d$p, alpha = 0.4,
    init_plan = G0, max_iter = 1L, tol_rel = 0, tol_abs = 0,
    verbose = FALSE, use_mixed_precision = FALSE
  )
  expect_false(isTRUE(all.equal(unname(exact0$plan), unname(d$p %o% rep(1 / 6, 6)),
                                tolerance = 1e-8)))
  expect_equal(unname(exact0$plan), unname(fast0$plan), tolerance = 5e-8)
})

test_that("symmetric exact and fast double paths agree", {
  d <- make_sr_problem(10L, 8L, seed = 21L)
  exact <- rfugw:::cpp_semirelaxed_fgw_exact_square(
    M = d$M, C1 = d$C1, C2 = d$C2, p = d$p, alpha = 0.5,
    symmetric = TRUE, init_plan = empty_M(),
    max_iter = 80L, tol_rel = 1e-9, tol_abs = 1e-9
  )
  fast <- rfugw:::cpp_semirelaxed_fgw_cg_square_fast(
    M = d$M, C1 = d$C1, C2 = d$C2, p = d$p, alpha = 0.5,
    init_plan = empty_M(), max_iter = 80L, tol_rel = 1e-9, tol_abs = 1e-9,
    verbose = FALSE, use_mixed_precision = FALSE
  )
  expect_equal(unname(exact$plan), unname(fast$plan), tolerance = 1e-7)
  expect_equal(exact$srfgw_dist, fast$srfgw_dist, tolerance = 1e-8)
})

test_that("mixed precision stays close to double on a larger symmetric problem", {
  d <- make_sr_problem(64L, 64L, seed = 42L)
  dbl <- rfugw:::cpp_semirelaxed_fgw_cg_square_fast(
    M = d$M, C1 = d$C1, C2 = d$C2, p = d$p, alpha = 0.5,
    init_plan = empty_M(), max_iter = 40L, tol_rel = 1e-6, tol_abs = 1e-6,
    verbose = FALSE, use_mixed_precision = FALSE
  )
  mix <- rfugw:::cpp_semirelaxed_fgw_cg_square_fast(
    M = d$M, C1 = d$C1, C2 = d$C2, p = d$p, alpha = 0.5,
    init_plan = empty_M(), max_iter = 40L, tol_rel = 1e-6, tol_abs = 1e-6,
    verbose = FALSE, use_mixed_precision = TRUE
  )
  expect_equal(mix$srfgw_dist, dbl$srfgw_dist, tolerance = 5e-4)
  expect_lt(sqrt(mean((mix$plan - dbl$plan)^2)), 5e-3)
  expect_equal(unname(rowSums(mix$plan)), d$p, tolerance = 1e-6)
})

test_that("empty M at alpha=1 matches an explicit zero feature cost", {
  d <- make_sr_problem(8L, 7L, seed = 9L)
  z <- matrix(0, nrow(d$C1), nrow(d$C2))
  empty_exact <- rfugw:::cpp_semirelaxed_fgw_exact_square(
    M = empty_M(), C1 = d$C1, C2 = d$C2, p = d$p, alpha = 1,
    symmetric = FALSE, init_plan = empty_M(),
    max_iter = 50L, tol_rel = 1e-9, tol_abs = 1e-9
  )
  zero_exact <- rfugw:::cpp_semirelaxed_fgw_exact_square(
    M = z, C1 = d$C1, C2 = d$C2, p = d$p, alpha = 1,
    symmetric = FALSE, init_plan = empty_M(),
    max_iter = 50L, tol_rel = 1e-9, tol_abs = 1e-9
  )
  empty_fast <- rfugw:::cpp_semirelaxed_fgw_cg_square_fast(
    M = empty_M(), C1 = d$C1, C2 = d$C2, p = d$p, alpha = 1,
    init_plan = empty_M(), max_iter = 50L, tol_rel = 1e-9, tol_abs = 1e-9,
    verbose = FALSE, use_mixed_precision = FALSE
  )
  zero_fast <- rfugw:::cpp_semirelaxed_fgw_cg_square_fast(
    M = z, C1 = d$C1, C2 = d$C2, p = d$p, alpha = 1,
    init_plan = empty_M(), max_iter = 50L, tol_rel = 1e-9, tol_abs = 1e-9,
    verbose = FALSE, use_mixed_precision = FALSE
  )
  expect_equal(empty_exact$srgw_dist, zero_exact$srgw_dist, tolerance = 1e-12)
  expect_equal(unname(empty_exact$plan), unname(zero_exact$plan), tolerance = 1e-12)
  expect_equal(empty_fast$srgw_dist, zero_fast$srgw_dist, tolerance = 1e-12)
  expect_equal(empty_exact$lin_loss, 0)
  expect_equal(empty_fast$lin_loss, 0)
})

test_that("entropic SR GW accepts an empty feature cost", {
  d <- make_sr_problem(8L, 6L, seed = 5L)
  empty <- rfugw:::cpp_entropic_semirelaxed_fgw_square(
    M = empty_M(), C1 = d$C1, C2 = d$C2, p = d$p, epsilon = 0.08,
    alpha = 1, symmetric = TRUE, use_mixed_precision = FALSE,
    max_iter = 80L, tol = 1e-8, check_every = 5L, init_plan = empty_M()
  )
  zero <- rfugw:::cpp_entropic_semirelaxed_fgw_square(
    M = matrix(0, nrow(d$C1), nrow(d$C2)), C1 = d$C1, C2 = d$C2, p = d$p,
    epsilon = 0.08, alpha = 1, symmetric = TRUE, use_mixed_precision = FALSE,
    max_iter = 80L, tol = 1e-8, check_every = 5L, init_plan = empty_M()
  )
  expect_equal(empty$srgw_dist, zero$srgw_dist, tolerance = 1e-10)
  expect_equal(unname(rowSums(empty$plan)), d$p, tolerance = 1e-10)
})
