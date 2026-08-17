test_that("ot_sinkhorn scaling and log agree on a well-conditioned problem", {
  M <- matrix(c(0, 1, 1, 0), 2, 2)
  p <- c(0.4, 0.6)
  q <- c(0.55, 0.45)
  sc <- ot_sinkhorn(M, p, q, epsilon = 0.1, method = "scaling", max_iter = 400L)
  lg <- ot_sinkhorn(M, p, q, epsilon = 0.1, method = "log", max_iter = 400L)
  expect_s3_class(sc, "rfugw_result")
  expect_equal(sc$ot_dist, lg$ot_dist, tolerance = 1e-6)
  expect_equal(rowSums(sc$plan), p, tolerance = 1e-6)
  expect_equal(colSums(sc$plan), q, tolerance = 1e-6)
  expect_equal(rfugw_value(sc), sc$ot_dist)
})

test_that("ot_emd recovers a permutation on a square assignment problem", {
  M <- matrix(c(
    0, 2, 3,
    2, 0, 4,
    3, 4, 0
  ), 3, 3, byrow = TRUE)
  p <- rep(1 / 3, 3)
  q <- p
  out <- ot_emd(M, p, q)
  expect_s3_class(out, "rfugw_result")
  expect_true(out$converged)
  expect_equal(rowSums(out$plan), p, tolerance = 1e-8)
  expect_equal(colSums(out$plan), q, tolerance = 1e-8)
  expect_equal(out$ot_dist, 0, tolerance = 1e-10)
  expect_equal(diag(out$plan), p, tolerance = 1e-8)
})

test_that("ot_sinkhorn_unbalanced relaxes mass and stays finite", {
  M <- matrix(c(0, 1, 1, 0), 2, 2)
  p <- c(0.9, 0.1)
  q <- c(0.1, 0.9)
  bal <- ot_sinkhorn(M, p, q, epsilon = 0.05, max_iter = 300L)
  unb <- ot_sinkhorn_unbalanced(M, p, q, epsilon = 0.05, rho = c(0.2, 0.2), max_iter = 300L)
  expect_s3_class(unb, "rfugw_result")
  expect_true(is.finite(unb$ot_dist))
  expect_true(unb$mass < 1)
  expect_lt(abs(sum(unb$plan) - 1), 1)
  expect_true(is.finite(bal$ot_dist))
})

test_that("linear OT rejects invalid inputs before compute", {
  M <- matrix(c(0, 1, 1, 0), 2, 2)
  expect_error(ot_sinkhorn(M, epsilon = 0), "epsilon")
  expect_error(ot_emd(M + NA), "finite")
  expect_error(ot_sinkhorn_unbalanced(M, rho = -1), "rho")
})
