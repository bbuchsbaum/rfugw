one_dimensional_ot_cost <- function(x, y, a, b, power) {
  ix <- order(x)
  iy <- order(y)
  x <- x[ix]
  y <- y[iy]
  a <- a[ix] / sum(a)
  b <- b[iy] / sum(b)
  i <- j <- 1L
  value <- 0
  while (i <= length(x) && j <= length(y)) {
    moved <- min(a[[i]], b[[j]])
    value <- value + moved * abs(x[[i]] - y[[j]])^power
    a[[i]] <- a[[i]] - moved
    b[[j]] <- b[[j]] - moved
    if (a[[i]] <= 1e-15) i <- i + 1L
    if (b[[j]] <= 1e-15) j <- j + 1L
  }
  value
}

test_that("Wasserstein p-cost and p-distance match independent 1-D oracles", {
  x <- c(0, 2, 5)
  y <- c(1, 3)
  a <- c(0.2, 0.5, 0.3)
  b <- c(0.6, 0.4)

  for (power in c(1, 2)) {
    oracle <- one_dimensional_ot_cost(x, y, a, b, power)
    cost <- ot_wasserstein_cost(
      x, y, p = power, source_weights = a, target_weights = b
    )
    distance <- ot_wasserstein_distance(
      x, y, p = power, source_weights = a, target_weights = b
    )
    expect_equal(rfugw_value(cost), oracle, tolerance = 1e-12)
    expect_equal(rfugw_value(distance), oracle^(1 / power), tolerance = 1e-12)
    expect_equal(cost$wasserstein_p_cost, distance$wasserstein_p_cost)
    expect_identical(cost$value_certification, "certified_exact_transport")
    expect_true(cost$metric_certified)
    diagnostics <- rfugw_residuals(cost)
    expect_equal(diagnostics$wasserstein_power, power)
    expect_true(diagnostics$value_certified)
  }

  w1 <- ot_wasserstein_distance(x, y, p = 1, source_weights = a, target_weights = b)
  expect_equal(w1$wasserstein_p_distance, w1$wasserstein_p_cost)
  expect_equal(w1$value_root, 1)
})

test_that("Wasserstein helpers obey permutation, symmetry, and scale laws", {
  x <- matrix(c(0, 0, 1, 2, 3, 1), 3L, byrow = TRUE)
  y <- matrix(c(1, 0, 2, 2), 2L, byrow = TRUE)
  a <- c(0.2, 0.3, 0.5)
  b <- c(0.6, 0.4)
  base <- ot_wasserstein_distance(
    x, y, p = 2, source_weights = a, target_weights = b
  )
  reverse <- ot_wasserstein_distance(
    y, x, p = 2, source_weights = b, target_weights = a
  )
  permuted <- ot_wasserstein_distance(
    x[c(3, 1, 2), ], y[c(2, 1), ], p = 2,
    source_weights = a[c(3, 1, 2)], target_weights = b[c(2, 1)]
  )
  scaled <- ot_wasserstein_distance(
    3 * x, 3 * y, p = 2, source_weights = a, target_weights = b
  )
  expect_equal(rfugw_value(reverse), rfugw_value(base), tolerance = 1e-12)
  expect_equal(rfugw_value(permuted), rfugw_value(base), tolerance = 1e-12)
  expect_equal(rfugw_value(scaled), 3 * rfugw_value(base), tolerance = 1e-11)
  expect_equal(
    scaled$wasserstein_p_cost,
    9 * base$wasserstein_p_cost,
    tolerance = 1e-11
  )
})

test_that("declared supplied-cost power controls the root convention", {
  distances <- matrix(c(0, 2, 1, 3, 4, 0), 3L, 2L, byrow = TRUE)
  squared <- distances^2
  a <- c(0.2, 0.3, 0.5)
  b <- c(0.6, 0.4)
  from_distance <- ot_wasserstein_distance(
    cost = distances, cost_power = 1, p = 2,
    source_weights = a, target_weights = b
  )
  from_squared <- ot_wasserstein_distance(
    cost = squared, cost_power = 2, p = 2,
    source_weights = a, target_weights = b
  )
  expect_equal(rfugw_value(from_distance), rfugw_value(from_squared))
  expect_false(from_distance$metric_certified)
  expect_identical(from_distance$metric, "user_supplied")

  scaled <- ot_wasserstein_distance(
    cost = 9 * squared, cost_power = 2, p = 2,
    source_weights = a, target_weights = b
  )
  expect_equal(rfugw_value(scaled), 3 * rfugw_value(from_squared))
  expect_equal(scaled$wasserstein_p_cost, 9 * from_squared$wasserstein_p_cost)
})

test_that("identity, duplicates, zeros, and unequal supports are handled", {
  x <- c(0, 0, 2)
  weights <- c(0.25, 0.25, 0.5)
  identity <- ot_wasserstein_distance(
    x, x, p = 2, source_weights = weights, target_weights = weights
  )
  expect_equal(rfugw_value(identity), 0, tolerance = 1e-12)

  unequal <- ot_wasserstein_cost(
    c(0, 1, 2), c(0, 3), p = 1,
    source_weights = c(0, 0.4, 0.6), target_weights = c(0.5, 0.5)
  )
  expect_true(unequal$converged)
  expect_equal(dim(unequal$plan), c(3L, 2L))
  expect_equal(rowSums(unequal$plan), c(0, 0.4, 0.6), tolerance = 1e-12)

  rescaled_weights <- ot_wasserstein_cost(
    c(0, 1, 2), c(0, 3), p = 1,
    source_weights = 10 * c(0, 0.4, 0.6),
    target_weights = 4 * c(0.5, 0.5)
  )
  expect_equal(
    rescaled_weights$wasserstein_p_cost,
    unequal$wasserstein_p_cost,
    tolerance = 1e-12
  )
  expect_equal(rescaled_weights$original_source_mass, 10)
  expect_equal(rescaled_weights$original_target_mass, 4)
  expect_identical(
    rescaled_weights$measure_normalization,
    "separate_probability"
  )
})

test_that("entropic plan estimates are distinct and fail closed", {
  x <- c(0, 1, 4)
  y <- c(0.5, 2, 5)
  exact <- ot_wasserstein_cost(x, y, p = 2)
  entropic <- ot_wasserstein_cost(
    x, y, p = 2, solver = "sinkhorn", epsilon = 2,
    sinkhorn_method = "log", max_iter = 5000L, tol = 1e-9
  )
  expect_identical(exact$estimate_kind, "exact")
  expect_identical(entropic$estimate_kind, "entropic_plan")
  expect_match(entropic$value_certification, "not_exact_wasserstein")
  expect_true(entropic$value_certified)
  expect_false(identical(rfugw_value(exact), rfugw_value(entropic)))

  expect_error(
    ot_wasserstein_cost(
      x, rev(y), p = 2, solver = "sinkhorn", epsilon = 0.01,
      sinkhorn_method = "log", max_iter = 1L, tol = 1e-15
    ),
    "uncertified"
  )
  inspected <- ot_wasserstein_cost(
    x, rev(y), p = 2, solver = "sinkhorn", epsilon = 0.01,
    sinkhorn_method = "log", max_iter = 1L, tol = 1e-15,
    allow_uncertified = TRUE
  )
  expect_false(inspected$value_certified)
  expect_identical(
    inspected$value_certification,
    "uncertified_plan_value_requested_explicitly"
  )
})

test_that("Wasserstein input semantics fail before solving", {
  expect_error(ot_wasserstein_cost(c(0, 1), c(1, 2), p = 0), "p")
  expect_error(ot_wasserstein_cost(cost = matrix(1, 2, 2)), "cost_power")
  expect_error(
    ot_wasserstein_cost(c(0, 1), c(1, 2), cost = matrix(1, 2, 2), cost_power = 1),
    "either raw supports or `cost`"
  )
  expect_error(
    ot_wasserstein_cost(cost = matrix(c(0, -1, 1, 0), 2), cost_power = 1),
    "nonnegative"
  )
  expect_error(
    ot_wasserstein_cost(c(0, 1), c(1, 2), source_weights = c(0, 0)),
    "positive source"
  )
  expect_error(
    ot_wasserstein_cost(c(0, 1), c(1, 2), epsilon = 0.1),
    "only operational"
  )
})
