moment_fugw_oracle_fixture <- function(seed = 4102L) {
  set.seed(seed)
  ns <- 4L
  nt <- 5L
  source_coordinates <- matrix(rnorm(ns * 2L), ns, 2L)
  target_coordinates <- matrix(rnorm(nt * 2L), nt, 2L)
  source_features <- matrix(rnorm(ns * 3L), ns, 3L)
  target_features <- matrix(rnorm(nt * 3L), nt, 3L)
  wx <- runif(ns)
  wy <- runif(nt)
  wx <- wx / sum(wx)
  wy <- wy / sum(wy)
  list(
    Cx = sqeuclidean_cost(source_coordinates),
    Cy = sqeuclidean_cost(target_coordinates),
    M = sqeuclidean_cost(source_features, target_features),
    wx = wx,
    wy = wy
  )
}

test_that("the direct Moment-FUGW oracle has no solver or moment dependency", {
  oracle_functions <- c(
    "moment_fugw_oracle_scalar_gkl",
    "moment_fugw_oracle_product_gkl",
    "moment_fugw_oracle_structure",
    "moment_fugw_oracle_validate",
    "moment_fugw_oracle"
  )
  used_names <- unique(unlist(lapply(oracle_functions, function(name) {
    all.names(body(get(name, mode = "function")), functions = TRUE)
  })))
  expect_length(
    intersect(
      used_names,
      c(
        "fugw_kl", "fugw_factorized", ".fugw_factorized_objective",
        ".fused_unbalanced_across_spaces_cost_kl",
        ".factorized_structure_moments", ".divergence_between_products"
      )
    ),
    0L
  )
})

test_that("the four-index oracle certifies both exact implementations", {
  fixture <- moment_fugw_oracle_fixture()
  controls <- list(
    wx = fixture$wx,
    wy = fixture$wy,
    reg_marginals = c(3, 7),
    epsilon = 0.3,
    feature_weight = 0.45,
    structure_weight = 1.2,
    max_iter = 80L,
    tol = 1e-9,
    max_iter_ot = 10000L,
    tol_ot = 1e-11,
    rescale_plan = TRUE
  )
  dense <- do.call(fugw_kl, c(list(
    Cx = cost_block(fixture$Cx),
    Cy = cost_block(fixture$Cy),
    M = cost_block(fixture$M),
    precision = "strict_double"
  ), controls))
  implicit <- do.call(fugw_factorized, c(list(
    Cx = fixture$Cx,
    Cy = fixture$Cy,
    M = fixture$M,
    block_size = 2L
  ), controls))
  sample <- transport_plan_materialize(implicit$plans$sample)
  feature <- transport_plan_materialize(implicit$plans$feature)
  oracle <- moment_fugw_oracle(
    cost_block(fixture$Cx),
    cost_block(fixture$Cy),
    cost_block(fixture$M),
    fixture$wx,
    fixture$wy,
    sample,
    feature,
    controls$reg_marginals,
    controls$epsilon,
    controls$feature_weight,
    controls$structure_weight
  )

  expect_true(dense$converged)
  expect_true(implicit$converged)
  expect_true(implicit$inner_uot_certified)
  expect_equal(oracle$fugw_cost, implicit$fugw_cost, tolerance = 2e-10)
  expect_equal(oracle$structure_unweighted,
               implicit$objective_decomposition$structure_unweighted,
               tolerance = 2e-10)
  expect_equal(oracle$feature_unweighted,
               implicit$objective_decomposition$feature_unweighted,
               tolerance = 2e-10)
  expect_equal(oracle$regularization,
               implicit$objective_decomposition$regularization,
               tolerance = 2e-10)
  expect_equal(oracle$fugw_cost, dense$fugw_cost, tolerance = 2e-9)
  expect_equal(sample, dense$pi_samp, tolerance = 2e-8)
  expect_equal(feature, dense$pi_feat, tolerance = 2e-8)
})

test_that("oracle fixtures kill the four declared objective mutants", {
  Cx <- matrix(c(0, .4, 1.1, .2), 2, 2, byrow = TRUE)
  Cy <- matrix(c(.1, .8, .3, 0), 2, 2, byrow = TRUE)
  M <- matrix(c(.2, 1.3, .7, .1), 2, 2, byrow = TRUE)
  wx <- c(.65, .35)
  wy <- c(.4, .6)
  sample <- matrix(c(.17, .08, .04, .12), 2, 2, byrow = TRUE)
  feature <- matrix(c(.11, .03, .09, .07), 2, 2, byrow = TRUE)
  args <- list(
    Cx = Cx, Cy = Cy, M = M, wx = wx, wy = wy,
    sample_plan = sample, feature_plan = feature,
    reg_marginals = c(2.5, 6), epsilon = .35,
    feature_weight = .4, structure_weight = 1.3
  )
  truth <- do.call(moment_fugw_oracle, args)$fugw_cost
  mutants <- do.call(moment_fugw_oracle_mutants, args)

  expect_setequal(
    names(mutants),
    c(
      "wrong_structure_cross_sign", "omitted_source_product_kl",
      "omitted_plan_kl_mass_terms", "normalized_coefficient_convention"
    )
  )
  expect_true(all(is.finite(mutants)))
  expect_true(all(abs(mutants - truth) > 1e-4))
})

test_that("randomized exact-factor smoke cases retain strict dense parity", {
  cases <- data.frame(
    seed = c(5101L, 5102L, 5103L),
    ns = c(3L, 5L, 6L),
    nt = c(4L, 4L, 7L),
    block_size = c(1L, 3L, 4L),
    zero_weight = c(FALSE, TRUE, FALSE),
    near_zero_weight = c(FALSE, FALSE, TRUE)
  )
  for (case in seq_len(nrow(cases))) {
    spec <- cases[case, ]
    set.seed(spec$seed)
    source_left <- matrix(runif(spec$ns * 2L, -.15, .15), spec$ns, 2L)
    source_right <- matrix(runif(spec$ns * 2L, -.15, .15), spec$ns, 2L)
    target_left <- matrix(runif(spec$nt * 3L, -.15, .15), spec$nt, 3L)
    target_right <- matrix(runif(spec$nt * 3L, -.15, .15), spec$nt, 3L)
    Cx <- factorized_cost(
      source_left, source_right,
      row = runif(spec$ns, .2, .5), column = runif(spec$ns, .2, .5)
    )
    Cy <- factorized_cost(
      target_left, target_right,
      row = runif(spec$nt, .2, .5), column = runif(spec$nt, .2, .5)
    )
    M <- factorized_cost(
      matrix(runif(spec$ns * 2L, -.2, .2), spec$ns, 2L),
      matrix(runif(spec$nt * 2L, -.2, .2), spec$nt, 2L),
      row = runif(spec$ns, .1, .3), column = runif(spec$nt, .1, .3)
    )
    wx <- runif(spec$ns)
    wy <- runif(spec$nt)
    if (spec$zero_weight) {
      wx[[1L]] <- 0
      wy[[spec$nt]] <- 0
    }
    if (spec$near_zero_weight) {
      wx[[1L]] <- 1e-12
      wy[[spec$nt]] <- 1e-13
    }
    wx <- wx / sum(wx)
    wy <- wy / sum(wy)
    controls <- list(
      wx = wx, wy = wy, reg_marginals = c(3, 6), epsilon = .35,
      feature_weight = .3, structure_weight = 1.1,
      max_iter = 80L, tol = 2e-9,
      max_iter_ot = 10000L, tol_ot = 1e-11,
      rescale_plan = TRUE
    )
    dense <- do.call(fugw_kl, c(list(
      Cx = cost_block(Cx), Cy = cost_block(Cy), M = cost_block(M),
      precision = "strict_double"
    ), controls))
    implicit <- do.call(fugw_factorized, c(list(
      Cx = Cx, Cy = Cy, M = M, block_size = spec$block_size
    ), controls))
    sample <- transport_plan_materialize(implicit$plans$sample)
    feature <- transport_plan_materialize(implicit$plans$feature)
    label <- paste("random case", case)

    expect_true(dense$converged, info = label)
    expect_true(implicit$converged, info = label)
    expect_true(implicit$inner_uot_certified, info = label)
    expect_equal(implicit$fugw_cost, dense$fugw_cost,
                 tolerance = 1e-8, info = label)
    expect_equal(sample, dense$pi_samp, tolerance = 1e-6, info = label)
    expect_equal(feature, dense$pi_feat, tolerance = 1e-6, info = label)
  }
})

test_that("the retained exact-factor matrix satisfies its declared gates", {
  evidence_path <- system.file(
    "bench", "moment-fugw-accuracy-evidence.csv", package = "rfugw"
  )
  if (!nzchar(evidence_path) || !file.exists(evidence_path)) {
    evidence_path <- testthat::test_path(
      "..", "..", "inst", "bench", "moment-fugw-accuracy-evidence.csv"
    )
  }
  evidence <- utils::read.csv(evidence_path, stringsAsFactors = FALSE)

  expect_equal(nrow(evidence), 12L)
  expect_setequal(unique(c(evidence$ns, evidence$nt)),
                  c(4L, 5L, 8L, 11L, 16L, 19L, 32L, 37L,
                    64L, 71L, 127L, 128L))
  expect_true(all(evidence$pass))
  expect_true(all(evidence$dense_converged))
  expect_true(all(evidence$implicit_converged))
  expect_true(all(evidence$inner_certified))
  expect_true(all(evidence$outer_certified))
  expect_true(any(evidence$zero_weight))
  expect_true(any(evidence$near_zero_weight))
  expect_true(all(!evidence$zero_weight | !evidence$near_zero_weight))
  expect_lte(max(evidence$objective_relative_error), 1e-8)
  expect_lte(max(evidence$component_relative_error), 1e-8)
  expect_lte(max(evidence$action_relative_error), 1e-6)
  expect_lte(max(evidence$adjoint_relative_error), 1e-6)
  expect_lte(max(evidence$mass_relative_error), 1e-7)
  expect_lte(max(evidence$moment_relative_error), 1e-7)
  expect_true(all(is.na(evidence$reject_reason) | evidence$reject_reason == ""))
})
