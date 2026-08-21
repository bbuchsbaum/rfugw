test_that("fractional public count parameters are rejected instead of truncated", {
  M <- C <- matrix(c(0, 1, 1, 0), 2, 2)
  X <- matrix(c(0, 0, 1, 1, 2, 2), 3, 2, byrow = TRUE)
  calls <- list(
    fgw_max_iter = function() fgw_entropic(M, C, C, max_iter = 1.5),
    sinkhorn_max_iter = function() ot_sinkhorn(M, max_iter = 2.5),
    partial_outer = function() partial_gromov_wasserstein(
      C, C, numItermax = 2.5
    ),
    partial_inner = function() entropic_partial_gromov_wasserstein(
      C, C, inner_max_iter = 2.5
    ),
    semirelaxed = function() semirelaxed_gromov_wasserstein(
      C, C, max_iter = 2.5
    ),
    sampled_budget = function() sampled_gromov_wasserstein(
      C, C, nb_samples_grad = 1.5
    ),
    sampled_iterations = function() sampled_gromov_wasserstein(
      C, C, max_iter = 2.5
    ),
    lowrank_rank = function() lowrank_gromov_wasserstein_samples(
      X, X, rank = 1.5
    ),
    ucoot_inner = function() unbalanced_co_optimal_transport(
      X, X, max_iter_ot = 2.5
    )
  )
  for (label in names(calls)) {
    expect_error(calls[[label]](), "integer", info = label)
  }
})

test_that("measure validation names finite, sign, and mass violations", {
  M <- matrix(c(0, 1, 1, 0), 2, 2)
  invalid <- list(
    na = c(NA_real_, 1),
    nan = c(NaN, 1),
    inf = c(Inf, 1),
    negative = c(-0.1, 1.1),
    zero_mass = c(0, 0)
  )
  for (label in names(invalid)) {
    pattern <- if (identical(label, "zero_mass")) "positive total mass" else
      "finite and nonnegative"
    expect_error(ot_emd(M, p = invalid[[label]]), pattern, info = label)
    expect_error(
      partial_gromov_wasserstein(M, M, p = invalid[[label]]),
      pattern,
      info = label
    )
  }
  expect_error(ot_sinkhorn(matrix(numeric(), 0, 0)), "nonempty")
  expect_error(ot_sinkhorn(matrix(c(0, Inf, 1, 0), 2, 2)), "`M` must be finite")
})

test_that("partial mass defaults and zero transport are explicit", {
  C <- matrix(c(0, 1, 1, 0), 2, 2)
  defaulted <- partial_gromov_wasserstein(
    C, C, numItermax = 3L, log = TRUE
  )
  zero <- partial_gromov_wasserstein(
    C, C, m = 0, numItermax = 3L, log = TRUE
  )
  expect_true(defaulted$transported_mass_defaulted)
  expect_identical(defaulted$transported_mass_target, 1)
  expect_equal(sum(defaulted$plan), 1, tolerance = 1e-10)
  expect_false(zero$transported_mass_defaulted)
  expect_identical(zero$transported_mass_target, 0)
  expect_equal(zero$plan, matrix(0, 2, 2), tolerance = 1e-12)
  expect_equal(zero$mass_residual, 0, tolerance = 0)
  expect_true(zero$mass_certified)
  expect_identical(zero$mass_certification, "certified")
  expect_error(
    partial_gromov_wasserstein(C, C, m = c(0.2, 0.3)),
    "`m` must be one finite number"
  )
  expect_error(
    partial_gromov_wasserstein(C, C, m = 1.1),
    "min\\(sum\\(p\\), sum\\(q\\)\\)"
  )
})

test_that("partial diagnostics refuse a mass certificate without a target", {
  out <- rfugw:::.attach_solver_diagnostics(
    list(plan = matrix(0.25, 2, 2), objective = 0),
    residual = 0,
    converged = TRUE,
    iterations = 1L,
    max_iter = 2L,
    p = c(0.5, 0.5),
    q = c(0.5, 0.5),
    feasibility = "partial"
  )
  expect_false(out$mass_certified)
  expect_identical(out$mass_certification, "not_certified_missing_target")
  expect_false(out$converged)
})

test_that("singleton, sparse-support, and extreme-imbalance transport is certified", {
  singleton <- ot_emd(matrix(7, 1, 1), 1, 1)
  singleton_sink <- ot_sinkhorn(
    matrix(7, 1, 1), 1, 1, epsilon = 0.1, max_iter = 20L, tol = 1e-10
  )
  expect_equal(singleton$plan, matrix(1, 1, 1), tolerance = 0)
  expect_true(singleton$converged)
  expect_equal(singleton_sink$plan, matrix(1, 1, 1), tolerance = 1e-12)

  M <- matrix(c(0, 2, 4, 1, 3, 0, 5, 2, 1), 3, 3)
  sparse <- ot_emd(M, c(0, 1e-12, 1), c(1, 0, 1e-12))
  expect_true(sparse$converged)
  expect_true(sparse$feasible)
  expect_equal(rowSums(sparse$plan), c(0, 1e-12, 1) / (1 + 1e-12), tolerance = 1e-10)
  expect_equal(colSums(sparse$plan), c(1, 0, 1e-12) / (1 + 1e-12), tolerance = 1e-10)
})

test_that("generalized KL handles empty transported support and reference zeros", {
  p <- c(1, 0)
  q <- c(0.25, 0.75)
  expect_equal(ot_kl(matrix(0, 2, 2), p, q), 1, tolerance = 0)
  supported <- matrix(c(0.25, 0.75, 0, 0), 2, 2, byrow = TRUE)
  expect_equal(ot_kl(supported, p, q), 0, tolerance = 1e-15)
  outside <- supported
  outside[2, 1] <- 1e-12
  expect_identical(ot_kl(outside, p, q), Inf)
})

test_that("symmetry detection uses combined absolute and relative tolerance", {
  large <- matrix(c(0, 1e8 + 5e-5, 1e8, 0), 2, 2, byrow = TRUE)
  small <- matrix(c(0, 1 + 5e-8, 1, 0), 2, 2, byrow = TRUE)
  exact <- matrix(c(0, 1, 1, 0), 2, 2)
  expect_true(rfugw:::.is_symmetric_cost(large))
  expect_false(rfugw:::.is_symmetric_cost(small))
  expect_true(rfugw:::.resolve_symmetric(TRUE, large, large))
  expect_error(
    rfugw:::.resolve_symmetric(TRUE, small, exact),
    "atol.*rtol.*scale"
  )
})

test_that("native FGW weight aliases are normalized and alpha-equivalent", {
  C1 <- matrix(c(0, 1, 1, 0), 2, 2)
  C2 <- matrix(c(0, 2, 2, 0), 2, 2)
  M <- matrix(c(0, 1, 2, 0), 2, 2)
  via_alpha <- fgw_entropic(
    M, C1, C2, alpha = 0.75, epsilon = 0.2, max_iter = 20L,
    precision = "strict_double"
  )
  via_weights <- fgw_entropic(
    M, C1, C2, feature_weight = 1, structure_weight = 3,
    epsilon = 0.2, max_iter = 20L, precision = "strict_double"
  )

  expect_equal(via_weights$plan, via_alpha$plan, tolerance = 1e-12)
  expect_equal(via_weights$fgw_dist, via_alpha$fgw_dist, tolerance = 1e-12)
  expect_identical(via_weights$alpha_convention, "structure_share")
  expect_equal(via_weights$feature_weight, 0.25)
  expect_equal(via_weights$structure_weight, 0.75)
  expect_true(via_weights$weights_normalized)
  expect_equal(
    via_weights$objective_decomposition$total,
    via_weights$fgw_dist,
    tolerance = 1e-10
  )

  expect_error(
    fgw_entropic(M, C1, C2, feature_weight = 1),
    "must be supplied together"
  )
  expect_error(
    fgw_entropic(
      M, C1, C2, alpha = 0.5, feature_weight = 1, structure_weight = 1
    ),
    "either `alpha`"
  )
  expect_error(
    fgw_exact_cg(M, C1, C2, feature_weight = 0, structure_weight = 0),
    "cannot both be zero"
  )
  exact_alpha <- fgw_exact_cg(M, C1, C2, alpha = 0.75, max_iter = 4L)
  exact_weights <- fused_gromov_wasserstein(
    M, C1, C2, feature_weight = 1, structure_weight = 3,
    max_iter = 4L
  )
  expect_equal(exact_weights$plan, exact_alpha$plan, tolerance = 1e-12)
  expect_identical(exact_weights$alpha_convention, "structure_share")
})

test_that("FUGW reports its distinct coefficient convention and decomposition", {
  C <- matrix(c(0, 1, 1, 0), 2, 2)
  M <- matrix(c(0, 1, 2, 0), 2, 2)
  legacy <- fugw_kl(C, C, M = M, alpha = 0.25, max_iter = 4L)
  explicit <- fugw_kl(
    C, C, M = M, feature_weight = 0.25, structure_weight = 1,
    max_iter = 4L
  )
  scaled <- fugw_kl(
    C, C, M = M, feature_weight = 0.25, structure_weight = 2,
    max_iter = 4L
  )

  expect_equal(explicit$pi_samp, legacy$pi_samp, tolerance = 1e-12)
  expect_equal(explicit$fugw_cost, legacy$fugw_cost, tolerance = 1e-12)
  expect_identical(
    legacy$alpha_convention,
    "feature_coefficient_with_unit_structure"
  )
  expect_identical(
    explicit$alpha_convention,
    "feature_coefficient_with_explicit_structure_coefficient"
  )
  expect_equal(scaled$structure_weight, 2)
  expect_equal(
    scaled$objective_decomposition$total,
    scaled$fugw_cost,
    tolerance = 1e-8
  )
})
