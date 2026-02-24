test_that("fgw_entropic matches POT fixture (square loss)", {
  fx <- read_fixture("fgw_entropic_square_fixture.json")

  out <- rfugw::fgw_entropic(
    M = fx$inputs$M,
    C1 = fx$inputs$C1,
    C2 = fx$inputs$C2,
    p = fx$inputs$p,
    q = fx$inputs$q,
    alpha = fx$params$alpha,
    epsilon = fx$params$epsilon,
    max_iter = fx$params$max_iter,
    tol = fx$params$tol,
    sinkhorn_max_iter = fx$params$sinkhorn_numItermax,
    sinkhorn_tol = fx$params$sinkhorn_stopThr,
    symmetric = TRUE,
    solver = fx$params$solver
  )

  expect_equal(dim(out$plan), dim(fx$outputs$plan))
  expect_equal(out$plan, fx$outputs$plan, tolerance = 1e-5)
  expect_equal(out$fgw_dist, fx$outputs$fgw_dist, tolerance = 1e-6)
  expect_equal(rowSums(out$plan), fx$inputs$p, tolerance = 1e-7)
  expect_equal(colSums(out$plan), fx$inputs$q, tolerance = 1e-7)
})
