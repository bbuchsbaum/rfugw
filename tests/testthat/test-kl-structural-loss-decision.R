test_that("KL structural prototype matches independent enumeration", {
  source(trust_test_resource(
    "numerical-trust", "kl-structural-loss-prototype.R"
  ), local = TRUE)
  C1 <- matrix(c(0.7, 0.3, 0.5, 0.3, 0.9, 0.4, 0.5, 0.4, 0.8), 3L)
  C2 <- matrix(c(0.6, 0.2, 0.45, 0.2, 0.75, 0.35, 0.45, 0.35, 0.95), 3L)
  T <- outer(c(0.2, 0.3, 0.5), c(0.4, 0.35, 0.25))

  for (floor in c(0, 1e-18, 1e-15, 1e-12)) {
    expect_equal(
      gw_kl_factorized(C1, C2, T, floor),
      gw_kl_enumerated(C1, C2, T, floor),
      tolerance = 1e-12
    )
  }

  C2[1, 2] <- C2[2, 1] <- 1e-16
  smoothed <- vapply(c(1e-18, 1e-15, 1e-12), function(floor) {
    expect_equal(
      gw_kl_factorized(C1, C2, T, floor),
      gw_kl_enumerated(C1, C2, T, floor),
      tolerance = 1e-11
    )
    gw_kl_factorized(C1, C2, T, floor)
  }, numeric(1))
  expect_gt(max(smoothed) - min(smoothed), 0.1)
})

test_that("zero-diagonal metric costs expose the KL domain failure", {
  source(trust_test_resource(
    "numerical-trust", "kl-structural-loss-prototype.R"
  ), local = TRUE)
  C1 <- abs(outer(c(0, 1, 3), c(0, 1, 3), "-"))
  C2 <- abs(outer(c(0, 2, 5), c(0, 2, 5), "-"))
  T <- outer(c(0.2, 0.3, 0.5), c(0.4, 0.35, 0.25))

  expect_identical(gw_kl_enumerated(C1, C2, T, 0), Inf)
  floored <- vapply(c(1e-18, 1e-15, 1e-12), function(floor) {
    expect_equal(
      gw_kl_factorized(C1, C2, T, floor),
      gw_kl_enumerated(C1, C2, T, floor),
      tolerance = 1e-10
    )
    gw_kl_factorized(C1, C2, T, floor)
  }, numeric(1))
  expect_true(all(is.finite(floored)))
  expect_gt(max(floored) - min(floored), 1)
})

test_that("public GW APIs reject KL structural loss with the decision boundary", {
  C <- matrix(c(0, 1, 1, 0), 2L)
  X <- matrix(c(0, 0, 1, 1, 2, 0), 3L, byrow = TRUE)
  calls <- list(
    entropic = function() entropic_gromov_wasserstein(C, C, loss_fun = "kl_loss"),
    exact = function() gromov_wasserstein(C, C, loss_fun = "kl_loss"),
    partial = function() partial_gromov_wasserstein(C, C, loss_fun = "kl_loss"),
    semirelaxed = function() semirelaxed_gromov_wasserstein(C, C, loss_fun = "kl_loss"),
    sampled = function() sampled_gromov_wasserstein(C, C, loss_fun = "kl_loss"),
    graphs = function() sampled_gw_from_graphs(C, C, loss_fun = "kl_loss"),
    barycenter = function() entropic_gromov_barycenters(
      2L, list(C, C), loss_fun = "kl_loss"
    ),
    coordinates = function() sampled_gromov_wasserstein_coords(
      X, X, loss_fun = "kl_loss"
    )
  )
  for (label in names(calls)) {
    expect_error(
      calls[[label]](),
      "deliberately unsupported.*zero-diagonal distance",
      info = label
    )
  }
})

test_that("KL structural decision and retained benchmark are complete", {
  decision <- trust_test_resource("kl-structural-loss-decision.md")
  text <- paste(readLines(decision, warn = FALSE), collapse = "\n")
  expect_match(text, "Decision: defer", fixed = TRUE)
  expect_match(text, "a log(a / b) - a + b", fixed = TRUE)
  expect_match(text, "zero-diagonal", fixed = TRUE)
  expect_match(text, "1e-18", fixed = TRUE)

  baseline <- utils::read.csv(
    bench_test_resource("kl-structural-loss-baseline.csv"),
    stringsAsFactors = FALSE
  )
  expect_setequal(
    unique(baseline$regime),
    c("positive", "near_zero", "zero_diagonal_distance")
  )
  finite <- is.finite(baseline$enumerated_objective) &
    is.finite(baseline$factorized_objective)
  expect_lt(max(baseline$absolute_difference[finite]), 1e-10)
})
