test_that("edge lists are canonical, sorted, and explicit about duplicates", {
  edges <- data.frame(
    source = c(3, 1, 2, 1, 3),
    target = c(1, 2, 3, 2, 1),
    weight = c(0.2, 0.1, 0.4, 0.3, 0)
  )
  plan <- as_transport_plan(edges, n_source = 3, n_target = 3)
  expected <- matrix(0, 3, 3)
  expected[1, 2] <- 0.4
  expected[2, 3] <- 0.4
  expected[3, 1] <- 0.2

  expect_s3_class(plan, "rfugw_transport_plan")
  expect_identical(plan$representation, "edge_list")
  expect_equal(as.matrix(plan), expected)
  expect_equal(plan$data$source, c(1L, 2L, 3L))
  expect_equal(plan$data$target, c(2L, 3L, 1L))
  expect_identical(plan$metadata$duplicate_policy, "sum")
  expect_error(
    as_transport_plan(edges, 3, 3, duplicates = "error"),
    "Duplicate"
  )
})

test_that("invalid edge support and weights fail before construction", {
  valid <- data.frame(source = 1, target = 1, weight = 0.5)
  expect_error(as_transport_plan(valid), "require")
  expect_error(
    as_transport_plan(transform(valid, source = 0), 2, 2),
    "source indices"
  )
  expect_error(
    as_transport_plan(transform(valid, target = 3), 2, 2),
    "target indices"
  )
  expect_error(
    as_transport_plan(transform(valid, source = 1.5), 2, 2),
    "source indices"
  )
  expect_error(
    as_transport_plan(transform(valid, weight = -1), 2, 2),
    "nonnegative"
  )
  expect_error(
    as_transport_plan(transform(valid, weight = Inf), 2, 2),
    "finite"
  )
  expect_error(
    as_transport_plan(list(source = 1, target = 1), 2, 2),
    "fields"
  )
})

test_that("dense, Matrix sparse, and edge plans agree without semantic drift", {
  skip_if_not_installed("Matrix")
  dense <- matrix(c(0.2, 0, 0.3, 0.5, 0, 0.4), 2L, 3L)
  sparse <- Matrix::Matrix(dense, sparse = TRUE)
  index <- which(dense > 0, arr.ind = TRUE)
  edges <- data.frame(
    source = index[, 1L], target = index[, 2L], weight = dense[index]
  )
  plans <- list(
    dense = as_transport_plan(dense),
    sparse = as_transport_plan(sparse),
    edges = as_transport_plan(edges, 2, 3)
  )
  target_values <- matrix(1:9, 3L, 3L)
  source_values <- matrix(1:6, 2L, 3L)
  target_points <- matrix(c(0, 0, 1, 0, 0, 2), 3L, 2L, byrow = TRUE)
  source_points <- matrix(c(0, 1, 2, 3), 2L, 2L, byrow = TRUE)
  cost <- matrix(c(0.1, 2, 1, 0.4, 3, 0.2), 2L, 3L)
  p <- rowSums(dense)
  q <- colSums(dense)

  for (plan in plans) {
    expect_equal(transport_plan_shape(plan), c(source = 2L, target = 3L))
    expect_equal(transport_plan_mass(plan), sum(dense))
    expect_equal(transport_plan_mass(plan, "source"), p)
    expect_equal(transport_plan_mass(plan, "target"), q)
    expect_equal(transport_plan_apply(plan, target_values), dense %*% target_values)
    expect_equal(transport_plan_adjoint(plan, source_values), t(dense) %*% source_values)
    expect_equal(
      transport_plan_barycentric(plan, target_points),
      ot_barycentric_project(dense, target_points)
    )
    expect_equal(
      transport_plan_barycentric(
        plan, source_points, orientation = "target_to_source",
        zero_mass = "zero"
      ),
      ot_barycentric_project(
        dense, source_points, orientation = "target_to_source",
        zero_mass = "zero"
      )
    )
    expect_equal(ot_linear_cost(cost, plan), sum(cost * dense))
    expect_equal(ot_entropy(plan), ot_entropy(dense))
    expect_equal(ot_kl(plan, p, q), ot_kl(dense, p, q))
    expect_no_error(
      ot_validate_plan(plan, p, q, mass = sum(dense), marginals = "balanced")
    )
    expect_equal(transport_plan_materialize(plan), dense)
  }
  expect_identical(
    transport_plan_representation(plans$dense)$representation,
    "dense_materialized"
  )
  expect_true(transport_plan_representation(plans$sparse)$sparse)
  expect_false(transport_plan_representation(plans$edges)$materialized)
})

test_that("zero transported rows and columns use explicit abstention", {
  plan <- as_transport_plan(
    data.frame(source = 1, target = 2, weight = 1), 3, 3
  )
  target_points <- matrix(c(0, 0, 2, 4, 8, 16), 3L, 2L, byrow = TRUE)
  source_points <- matrix(c(1, 2, 3, 4, 5, 6), 3L, 2L, byrow = TRUE)
  forward_nan <- transport_plan_barycentric(plan, target_points, zero_mass = "nan")
  forward_zero <- transport_plan_barycentric(plan, target_points, zero_mass = "zero")
  reverse_nan <- transport_plan_barycentric(
    plan, source_points, orientation = "target_to_source", zero_mass = "nan"
  )

  expect_equal(forward_nan[1, ], target_points[2, ])
  expect_true(all(is.nan(forward_nan[2:3, ])))
  expect_equal(forward_zero[2:3, ], matrix(0, 2L, 2L))
  expect_equal(reverse_nan[2, ], source_points[1, ])
  expect_true(all(is.nan(reverse_nan[c(1, 3), ])))
  expect_error(
    transport_plan_barycentric(plan, target_points, zero_mass = "error"),
    "zero transported mass"
  )
  expect_false(any(forward_zero[2:3, ] == colMeans(target_points)))
})

test_that("implicit operator distinguishes adjoint, barycentric reverse, and inverse", {
  dense <- matrix(c(0.3, 0, 0.2, 0.5, 0, 0.4), 2L, 3L)
  operator <- transport_operator(
    2, 3,
    apply = function(x) dense %*% x,
    adjoint = function(x) t(dense) %*% x,
    source_mass = rowSums(dense),
    target_mass = colSums(dense),
    materialize = function() dense,
    metadata = list(origin = "test_fixture")
  )
  x <- matrix(1:6, 3L, 2L)
  y <- matrix(1:4, 2L, 2L)

  expect_identical(
    transport_plan_representation(operator)$representation,
    "implicit_operator"
  )
  expect_equal(transport_plan_apply(operator, x), dense %*% x)
  expect_equal(transport_plan_adjoint(operator, y), t(dense) %*% y)
  expect_equal(
    transport_plan_barycentric(
      operator, y, orientation = "target_to_source", zero_mass = "zero"
    ),
    ot_barycentric_project(
      dense, y, orientation = "target_to_source", zero_mass = "zero"
    )
  )
  expect_equal(transport_plan_materialize(operator), dense)
  expect_no_error(
    ot_validate_plan(
      operator, rowSums(dense), colSums(dense),
      mass = sum(dense), marginals = "balanced"
    )
  )
  expect_error(ot_linear_cost(dense, operator), "enumerable support")

  no_materializer <- transport_operator(
    2, 3, function(x) dense %*% x, function(x) t(dense) %*% x,
    rowSums(dense), colSums(dense)
  )
  expect_error(transport_plan_materialize(no_materializer), "no materializer")
})

test_that("operator callbacks and mass contracts fail closed", {
  expect_error(
    transport_operator(2, 2, identity, NULL, c(0.5, 0.5), c(0.5, 0.5)),
    "functions"
  )
  expect_error(
    transport_operator(2, 2, identity, identity, c(1, 0), c(0.5, 0.4)),
    "equal totals"
  )
  bad <- transport_operator(
    2, 2, function(x) numeric(), function(x) numeric(),
    c(0.5, 0.5), c(0.5, 0.5)
  )
  expect_error(transport_plan_apply(bad, c(1, 2)), "incompatible")
  expect_error(transport_plan_adjoint(bad, c(1, 2)), "incompatible")
})

test_that("pruning reports lost mass and invalidates result certificates", {
  dense <- matrix(c(0.5, 0.01, 0.2, 0.29), 2L)
  result <- structure(
    list(
      plan = as_transport_plan(dense),
      status = "converged", converged = TRUE,
      feasible = TRUE, objective_consistent = TRUE,
      objective_components_consistent = TRUE,
      residual = 0, error = 0, formulation = "fixture", backend = "fixture"
    ),
    class = "rfugw_result"
  )
  pruned <- transport_plan_prune(result, min_weight = 0.05)
  info <- transport_plan_representation(rfugw_plan(pruned))

  expect_equal(pruned$pruning_lost_mass, 0.01)
  expect_true(pruned$certificate_invalidated_by_pruning)
  expect_identical(pruned$status, "requires_revalidation")
  expect_false(pruned$converged)
  expect_false(pruned$feasible)
  expect_false(pruned$objective_consistent)
  expect_false(pruned$objective_components_consistent)
  expect_true(is.infinite(pruned$residual))
  expect_true(info$pruned)
  expect_equal(info$lost_mass, 0.01)
  expect_equal(transport_plan_mass(rfugw_plan(pruned)), 0.99)
  expect_no_error(ot_validate_plan(rfugw_plan(pruned), marginals = "relaxed"))

  unchanged <- transport_plan_prune(result, min_weight = 0)
  expect_equal(unchanged$pruning_lost_mass, 0)
  expect_true(unchanged$converged)
})

test_that("result accessors and serialization preserve sparse behavior", {
  skip_if_not_installed("Matrix")
  dense <- matrix(c(0.4, 0, 0.1, 0.5), 2L)
  wrapped <- as_transport_plan(Matrix::Matrix(dense, sparse = TRUE))
  result <- structure(
    list(
      plan = wrapped, status = "converged", converged = TRUE,
      formulation = "fixture", backend = "fixture", residual = 0,
      max_iter = 1L, iterations = 1L, ot_dist = 0
    ),
    class = "rfugw_result"
  )
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(result, path)
  restored <- readRDS(path)

  expect_s3_class(rfugw_plan(restored), "rfugw_transport_plan")
  expect_true(transport_plan_representation(rfugw_plan(restored))$sparse)
  expect_equal(rfugw_plan(restored, materialize = TRUE), dense)
  expect_equal(
    transport_plan_apply(rfugw_plan(restored), c(1, 2)),
    as.numeric(dense %*% c(1, 2))
  )
  expect_equal(dim(rfugw_plan(restored)), c(2L, 2L))
  expect_match(paste(capture.output(print(rfugw_plan(restored))), collapse = "\n"), "sparse")

  legacy <- result
  legacy$plan <- dense
  expect_true(is.matrix(rfugw_plan(legacy)))
  expect_identical(rfugw_plan(legacy), dense)
})

test_that("representative edge operations do not allocate a dense n-by-m plan", {
  n <- 20000L
  edges <- data.frame(
    source = rep(seq_len(n), each = 2L),
    target = c(rbind(seq_len(n), c(seq.int(2L, n), 1L))),
    weight = rep(0.5 / n, 2L * n)
  )
  plan <- as_transport_plan(edges, n, n)
  dense_bytes <- as.numeric(n) * as.numeric(n) * 8
  profile <- tempfile(fileext = ".mem")
  on.exit(unlink(profile), add = TRUE)
  Rprofmem(profile)
  applied <- transport_plan_apply(plan, rep(1, n))
  source_mass <- transport_plan_mass(plan, "source")
  Rprofmem(NULL)
  allocations <- readLines(profile, warn = FALSE)
  bytes <- suppressWarnings(as.numeric(sub(" .*", "", allocations)))

  expect_equal(applied, rep(1 / n, n), tolerance = 1e-15)
  expect_equal(source_mass, rep(1 / n, n), tolerance = 1e-15)
  expect_lt(as.numeric(object.size(plan)), dense_bytes / 100)
  expect_lt(max(bytes[is.finite(bytes)], 0), dense_bytes / 100)
  expect_false(transport_plan_representation(plan)$materialized)
})
