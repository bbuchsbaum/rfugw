test_that("fgw_exact_cg matches POT exact fixture", {
  fx <- read_fixture("fgw_exact_square_fixture.json")

  out <- rfugw::fgw_exact_cg(
    M = fx$inputs$M,
    C1 = fx$inputs$C1,
    C2 = fx$inputs$C2,
    p = fx$inputs$p,
    q = fx$inputs$q,
    alpha = fx$params$alpha,
    max_iter = fx$params$max_iter,
    tol_rel = fx$params$tol_rel,
    tol_abs = fx$params$tol_abs,
    lp_solver = "cpp_transport"
  )

  expect_equal(out$fgw_dist, fx$outputs$fgw_dist, tolerance = 1e-5)
  expect_equal(out$plan, fx$outputs$plan, tolerance = 1e-4)
  expect_equal(rowSums(out$plan), fx$inputs$p, tolerance = 2e-6)
  expect_equal(colSums(out$plan), fx$inputs$q, tolerance = 2e-6)
})
