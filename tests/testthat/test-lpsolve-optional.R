test_that("lpSolve-backed exact CG fails with an install route when missing", {
  skip_if(requireNamespace("lpSolve", quietly = TRUE))
  d <- list(
    C1 = matrix(c(0, 1, 1, 0), 2, 2),
    C2 = matrix(c(0, 1, 1, 0), 2, 2),
    M = matrix(c(0, 1, 1, 0), 2, 2)
  )
  expect_error(
    fgw_exact_cg(d$M, d$C1, d$C2, lp_solver = "lp_transport", max_iter = 1L),
    "install.packages\\(\"lpSolve\"\\)"
  )
  expect_no_error(
    fgw_exact_cg(d$M, d$C1, d$C2, lp_solver = "cpp_transport", max_iter = 2L)
  )
})

test_that("default partial solvers do not need lpSolve; lp_matrix still does", {
  C1 <- matrix(c(0, 1, 1, 0), 2, 2)
  C2 <- matrix(c(0, 1, 1, 0), 2, 2)
  out <- partial_gromov_wasserstein(C1, C2, m = 0.5, log = TRUE)
  expect_true(is.finite(out$partial_gw_dist))
  expect_equal(sum(out$plan), 0.5, tolerance = 5e-6)
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    expect_error(
      partial_gromov_wasserstein(C1, C2, m = 0.5, log = TRUE, lp_solver = "lp_matrix"),
      "lpSolve"
    )
  } else {
    lp <- partial_gromov_wasserstein(
      C1, C2, m = 0.5, log = TRUE, lp_solver = "lp_matrix", numItermax = 40L
    )
    expect_true(is.finite(lp$partial_gw_dist))
  }
  ent <- entropic_partial_gromov_wasserstein(C1, C2, m = 0.5, reg = 0.1, log = TRUE)
  expect_true(is.finite(ent$partial_gw_dist))
})
