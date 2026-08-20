source(trust_test_resource("numerical-trust", "mutation-proof-lib.R"), local = TRUE)

test_that("the bounded mutation harness kills every reviewed defect class", {
  report <- run_rfugw_mutation_proof()
  expect_true(report$all_killed)
  expect_length(report$cases, 6L)
  expect_true(all(vapply(report$cases, `[[`, logical(1), "killed")))
  expect_true(all(vapply(report$cases, function(x) x$runtime_ms >= 0, logical(1))))
  expect_setequal(
    vapply(report$cases, `[[`, character(1), "mutant"),
    c(
      "wrong_C2_transpose", "false_early_simplex_convergence",
      "mixed_tolerance_floor_1e-6", "ignored_inner_residual",
      "omitted_generalized_KL_mass_terms",
      "finite_positive_mass_outside_zero_support"
    )
  )
})

test_that("objective-only and symmetric-only checks miss the transpose mutant", {
  evidence <- run_rfugw_mutation_proof()$insufficient_test_demonstration
  expect_false(evidence$scalar_objective_detects_transpose_mutant)
  expect_false(evidence$symmetric_gradient_detects_transpose_mutant)
  expect_true(evidence$asymmetric_gradient_detects_transpose_mutant)
})
