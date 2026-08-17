make_ucoot_problem <- function(ns, nf, nt, mf, seed) {
  set.seed(seed)
  list(
    X = matrix(stats::rnorm(ns * nf), ns, nf),
    Y = matrix(stats::rnorm(nt * mf), nt, mf)
  )
}

test_that("C++ UCOOT matches the R BCD core on independent and joint modes", {
  d <- make_ucoot_problem(6L, 4L, 5L, 3L, seed = 11L)
  wx_s <- rep(1 / 6, 6)
  wx_f <- rep(1 / 4, 4)
  wy_s <- rep(1 / 5, 5)
  wy_f <- rep(1 / 3, 3)
  wxy_s <- wx_s %o% wy_s
  wxy_f <- wx_f %o% wy_f
  for (joint in c(FALSE, TRUE)) {
    cpp <- rfugw:::cpp_ucoot_kl(
      d$X, d$Y, wx_s, wx_f, wy_s, wy_f, c(10, 8), c(0.05, 0.03),
      matrix(numeric(0), 0, 0), matrix(numeric(0), 0, 0),
      wxy_s, wxy_f, joint, TRUE, 25L, 1e-7, 150L, 1e-7, FALSE
    )
    ref <- rfugw:::.ucoot_kl_r_core(
      d$X, d$Y, wx_s, wx_f, wy_s, wy_f, wxy_s, wxy_f,
      10, 8, 0.05, if (joint) 0.05 else 0.03, NULL, NULL,
      wxy_s, wxy_f, if (joint) "joint" else "independent",
      TRUE, 25L, 1e-7, 150L, 1e-7, TRUE, FALSE
    )
    expect_equal(cpp$ucoot_cost, ref$ucoot_cost, tolerance = 5e-2)
    expect_equal(unname(cpp$pi_samp), unname(ref$pi_samp), tolerance = 5e-2)
    expect_true(all(is.finite(cpp$pi_samp)))
    expect_true(all(is.finite(cpp$pi_feat)))
  }
})

test_that("warm-start telemetry records use or a safe fallback", {
  d <- make_ucoot_problem(5L, 4L, 5L, 4L, seed = 3L)
  out <- rfugw::unbalanced_co_optimal_transport(
    d$X, d$Y, max_iter = 12L, tol = 1e-8, max_iter_ot = 80L, tol_ot = 1e-7,
    log = TRUE
  )
  expect_true(is.logical(out$inner_warm_feat))
  expect_true(is.logical(out$inner_warm_samp))
  expect_true(is.logical(out$inner_warm_fallback_feat))
  expect_true(is.logical(out$inner_warm_fallback_samp))
  expect_equal(length(out$inner_warm_feat), out$iterations)
  expect_equal(
    out$inner_iters_total,
    sum(out$inner_iters_feat) + sum(out$inner_iters_samp)
  )
  if (out$iterations >= 2L) {
    later <- seq.int(2L, out$iterations)
    used_or_safe <- out$inner_warm_feat[later] | out$inner_warm_fallback_feat[later]
    expect_true(all(used_or_safe))
  }
})

test_that("supplied UCOOT init_pi is consumed", {
  d <- make_ucoot_problem(4L, 3L, 4L, 3L, seed = 2L)
  wx_s <- rep(1 / 4, 4)
  wy_s <- rep(1 / 4, 4)
  wx_f <- rep(1 / 3, 3)
  wy_f <- rep(1 / 3, 3)
  init <- list(
    pi_samp = (wx_s %o% wy_s) * 0.4 + diag(4) * 0.05,
    pi_feat = wx_f %o% wy_f
  )
  warm <- rfugw::unbalanced_co_optimal_transport(
    d$X, d$Y, init_pi = init, max_iter = 1L, tol = 0, log = TRUE
  )
  cold <- rfugw::unbalanced_co_optimal_transport(
    d$X, d$Y, max_iter = 1L, tol = 0, log = TRUE
  )
  expect_false(isTRUE(all.equal(unname(warm$pi_samp), unname(cold$pi_samp),
                                tolerance = 1e-8)))
})

test_that("POT joint and independent UCOOT differentials stay in range", {
  fx <- read_fixture("ucoot_pot_modes.json")
  ind <- rfugw::unbalanced_co_optimal_transport(
    X = fx$inputs$X, Y = fx$inputs$Y,
    reg_marginals = unlist(fx$params$reg_marginals),
    epsilon = unlist(fx$params$epsilon),
    max_iter = as.integer(fx$params$max_iter),
    tol = fx$params$tol,
    max_iter_ot = as.integer(fx$params$max_iter_ot),
    tol_ot = fx$params$tol_ot,
    log = TRUE
  )
  joint <- rfugw::fused_unbalanced_across_spaces_divergence(
    X = fx$inputs$X, Y = fx$inputs$Y,
    reg_marginals = unlist(fx$params$reg_marginals),
    epsilon = unlist(fx$params$epsilon),
    reg_type = "joint",
    max_iter = as.integer(fx$params$max_iter),
    tol = fx$params$tol,
    max_iter_ot = as.integer(fx$params$max_iter_ot),
    tol_ot = fx$params$tol_ot,
    log = TRUE
  )
  expect_equal(ind$ucoot_cost, fx$independent$ucoot_cost, tolerance = 2e-2)
  expect_equal(joint$ucoot_cost, fx$joint$ucoot_cost, tolerance = 2e-2)
  expect_equal(unname(ind$pi_samp), unname(fx$independent$pi_samp), tolerance = 2e-2)
  expect_equal(unname(joint$pi_samp), unname(fx$joint$pi_samp), tolerance = 2e-2)
})
