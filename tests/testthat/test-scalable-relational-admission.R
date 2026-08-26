test_that("scalable relational OT remains fail-closed", {
  decision <- paste(readLines(
    trust_test_resource("scalable-relational-ot-decision.md"),
    warn = FALSE
  ), collapse = "\n")
  expect_match(decision, "No relational-OT approximation is promoted")
  expect_match(decision, "five double matrices")
  expect_match(decision, "1.0 GB at `n = 5,000`")
  expect_match(decision, "operator structure cost")
  expect_match(decision, "remain experimental")
})

test_that("admission evidence exposes the quadratic coupling boundary", {
  evidence <- utils::read.csv(
    bench_test_resource("scalable-relational-admission.csv"),
    stringsAsFactors = FALSE
  )
  required <- c(
    "candidate", "evidence_type", "n_source", "n_target",
    "budget_source", "budget_target", "median_seconds",
    "r_visible_allocated_bytes", "relative_objective_gap_to_dense_sampled",
    "returned_plan_logical_bytes", "minimum_simultaneous_dense_workspaces",
    "dense_workspace_lower_bound_bytes", "runtime_doubling_slope",
    "workspace_doubling_slope", "decision"
  )
  expect_true(all(required %in% names(evidence)))
  expect_true(all(evidence$decision == "defer"))
  expect_true(all(evidence$minimum_simultaneous_dense_workspaces == 5L))
  expect_true(all(evidence$workspace_doubling_slope == 2))
  expect_equal(
    evidence$dense_workspace_lower_bound_bytes,
    5 * 8 * evidence$n_source * evidence$n_target
  )
  expect_equal(
    evidence$returned_plan_logical_bytes,
    8 * evidence$n_source * evidence$n_target
  )
  projected <- evidence[evidence$evidence_type == "source_audit_projection", ]
  expect_gte(min(projected$dense_workspace_lower_bound_bytes), 1e9)
})

test_that("current coordinate sampled path returns a dense uncertified plan", {
  X1 <- cbind(0:5, c(0, 1, 0, 1, 0, 1))
  X2 <- cbind(0:6, c(1, 0, 1, 0, 1, 0, 1))
  out <- sampled_gromov_wasserstein_coords(
    X1, X2,
    nb_samples_grad = c(2L, 1L), epsilon = 0.2,
    max_iter = 2L, sinkhorn_max_iter = 100L,
    random_state = 11L, log = TRUE
  )
  expect_true(is.matrix(out$plan))
  expect_identical(dim(out$plan), c(6L, 7L))
  expect_false(inherits(out$plan, "rfugw_transport_operator"))
  expect_false("certification" %in% names(out))
})
