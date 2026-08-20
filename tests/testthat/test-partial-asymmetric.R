make_partial_asymmetric_fixture <- function() {
  C1 <- matrix(c(0, 0.2, 0.8, 0), 2, 2, byrow = TRUE)
  C2 <- matrix(c(0, 0.7, 0.1, 0), 2, 2, byrow = TRUE)
  list(
    C1 = C1,
    C2 = C2,
    M = matrix(c(0.1, 0.8, 0.6, 0.2), 2, 2),
    p = c(0.4, 0.6),
    q = c(0.55, 0.45),
    m = 0.7
  )
}

test_that("public exact and entropic partial paths accept asymmetric costs", {
  d <- make_partial_asymmetric_fixture()

  outputs <- list(
    partial_gromov_wasserstein(
      d$C1, d$C2, d$p, d$q, d$m, symmetric = FALSE, log = TRUE,
      numItermax = 20L
    ),
    partial_fused_gromov_wasserstein(
      d$M, d$C1, d$C2, d$p, d$q, d$m, symmetric = FALSE, log = TRUE,
      numItermax = 20L
    ),
    entropic_partial_gromov_wasserstein(
      d$C1, d$C2, d$p, d$q, m = d$m, symmetric = FALSE, log = TRUE,
      numItermax = 20L
    ),
    entropic_partial_fused_gromov_wasserstein(
      d$M, d$C1, d$C2, d$p, d$q, m = d$m, symmetric = FALSE, log = TRUE,
      numItermax = 20L
    )
  )
  for (out in outputs) {
    expect_true(all(is.finite(out$plan)))
    expect_true(all(out$plan >= 0))
    expect_equal(sum(out$plan), d$m, tolerance = 2e-7)
    expect_true(all(rowSums(out$plan) <= d$p + 2e-7))
    expect_true(all(colSums(out$plan) <= d$q + 2e-7))
  }
})

test_that("native exact and entropic partial entry points use asymmetric algebra", {
  d <- make_partial_asymmetric_fixture()
  G0 <- (d$p %o% d$q) * d$m

  exact <- rfugw:::cpp_partial_fgw_exact_square(
    d$M, d$C1, d$C2, d$p, d$q, d$m, 0.5, FALSE,
    G0, 2L, 1e-8, 1L, 100L, 1e-10
  )
  entropic <- rfugw:::cpp_partial_fgw_entropic_square(
    d$M, d$C1, d$C2, d$p, d$q, d$m, 0.1, 0.5, FALSE,
    G0, 2L, 1e-8, 20L, 1e-10, 1L
  )

  expect_true(all(is.finite(exact$plan)))
  expect_true(all(is.finite(entropic$plan)))
  expect_equal(sum(exact$plan), d$m, tolerance = 2e-8)
  expect_equal(sum(entropic$plan), d$m, tolerance = 2e-7)
})

test_that("symmetric partial paths remain available after general-path repair", {
  C <- matrix(c(0, 1, 1, 0), 2, 2)
  p <- q <- c(0.5, 0.5)

  exact <- partial_gromov_wasserstein(
    C, C, p, q, m = 0.5, numItermax = 10L, symmetric = TRUE, log = TRUE
  )
  entropic <- entropic_partial_gromov_wasserstein(
    C, C, p, q, m = 0.5, reg = 0.1, numItermax = 10L,
    symmetric = TRUE, log = TRUE
  )

  expect_equal(sum(exact$plan), 0.5, tolerance = 1e-8)
  expect_equal(sum(entropic$plan), 0.5, tolerance = 1e-7)
  expect_true(all(is.finite(exact$plan)))
  expect_true(all(is.finite(entropic$plan)))
})
