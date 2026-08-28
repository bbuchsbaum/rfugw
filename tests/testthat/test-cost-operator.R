test_that("affine-bilinear cost constructors reproduce dense definitions", {
  set.seed(731)
  source <- matrix(rnorm(24), 6, 4)
  target <- matrix(rnorm(20), 5, 4)
  squared <- sqeuclidean_cost(source, target)
  squared_dense <- outer(rowSums(source^2), rowSums(target^2), "+") -
    2 * tcrossprod(source, target)

  expect_s3_class(squared, "rfugw_cost_operator")
  expect_equal(cost_shape(squared), c(source = 6L, target = 5L))
  expect_equal(cost_block(squared), squared_dense, tolerance = 1e-13)
  expect_equal(
    cost_block(squared, c(6, 2, 2), c(5, 1)),
    squared_dense[c(6, 2, 2), c(5, 1), drop = FALSE],
    tolerance = 1e-13
  )

  cosine <- cosine_cost(source, target)
  source_norm <- source / sqrt(rowSums(source^2))
  target_norm <- target / sqrt(rowSums(target^2))
  expect_equal(
    cost_block(cosine), 1 - tcrossprod(source_norm, target_norm),
    tolerance = 1e-13
  )

  correlation <- correlation_cost(source, target)
  source_centered <- source - rowMeans(source)
  target_centered <- target - rowMeans(target)
  source_centered <- source_centered / sqrt(rowSums(source_centered^2))
  target_centered <- target_centered / sqrt(rowSums(target_centered^2))
  expect_equal(
    cost_block(correlation),
    1 - tcrossprod(source_centered, target_centered),
    tolerance = 1e-13
  )
  expect_true(isTRUE(cost_provenance(squared)$exact))
})

test_that("weighted reductions and algebraic composition stay matrix-free", {
  set.seed(732)
  left <- matrix(rnorm(21), 7, 3)
  right <- matrix(rnorm(15), 5, 3)
  row <- rnorm(7)
  column <- rnorm(5)
  base <- factorized_cost(left, right, row, column, constant = 0.4)
  dense <- outer(row + 0.4, column, "+") + tcrossprod(left, right)
  target_weight <- runif(5)
  source_weight <- runif(7)

  expect_equal(
    cost_row_reduce(base, target_weight),
    as.numeric(dense %*% target_weight),
    tolerance = 1e-13
  )
  expect_equal(
    cost_col_reduce(base, source_weight),
    as.numeric(crossprod(dense, source_weight)),
    tolerance = 1e-13
  )
  expect_equal(cost_row_reduce(base), rowSums(dense), tolerance = 1e-13)
  expect_equal(cost_col_reduce(base), colSums(dense), tolerance = 1e-13)

  combined <- sum_costs(scaled_cost(base, 0.25), scaled_cost(base, -0.5))
  expect_equal(cost_block(combined), -0.25 * dense, tolerance = 1e-13)
  expect_identical(cost_provenance(combined)$representation,
                   "affine_bilinear_sum")
})

test_that("generic callbacks evaluate only requested tiles", {
  calls <- list()
  evaluator <- function(rows, columns) {
    calls[[length(calls) + 1L]] <<- c(length(rows), length(columns))
    outer(rows^2, columns, "+")
  }
  operator <- block_cost(
    9, 8, evaluator, exact = TRUE,
    provenance = list(origin = "test_callback")
  )
  dense <- outer(seq_len(9)^2, seq_len(8), "+")

  expect_equal(cost_block(operator, 2:4, c(1, 7)), dense[2:4, c(1, 7)])
  expect_equal(tail(calls, 1)[[1L]], c(3, 2))
  expect_equal(
    cost_row_reduce(operator, seq_len(8), block_size = 3),
    as.numeric(dense %*% seq_len(8))
  )
  expect_equal(
    cost_col_reduce(operator, seq_len(9), block_size = 4),
    as.numeric(crossprod(dense, seq_len(9)))
  )
  expect_true(all(vapply(calls, function(x) all(x <= 4), logical(1))))
})

test_that("cost operators fail closed on malformed or degenerate inputs", {
  expect_error(
    factorized_cost(matrix(1, 2, 1), matrix(1, 3, 2)),
    "same number of columns"
  )
  expect_error(
    factorized_cost(matrix(c(1, Inf), 2, 1), matrix(1, 3, 1)),
    "finite numeric matrix"
  )
  expect_error(cost_block(sqeuclidean_cost(matrix(1:4, 2)), 0, 1),
               "invalid indices")
  expect_error(cosine_cost(matrix(c(0, 0, 1, 2), 2, 2, byrow = TRUE)),
               "zero-norm")
  expect_error(correlation_cost(matrix(1:3, 3, 1)), "at least two")
  expect_no_error(block_cost(2, 2, function(rows, columns) numeric()))
  malformed <- block_cost(2, 2, function(rows, columns) numeric())
  expect_error(cost_block(malformed), "requested shape")
})

test_that("large factorized costs have linear storage", {
  n <- 20000L
  coordinates <- cbind(seq_len(n) / n, (seq_len(n) %% 17L) / 17)
  operator <- sqeuclidean_cost(coordinates)
  dense_bytes <- as.double(n) * as.double(n) * 8

  expect_lt(as.numeric(object.size(operator)), dense_bytes / 1000)
  expect_equal(cost_shape(operator), c(source = n, target = n))
  expect_equal(dim(cost_block(operator, 1:8, 1:9)), c(8L, 9L))
  expect_equal(cost_provenance(operator)$rank, 2L)
})
