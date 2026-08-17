# Linear-OT certification.
#
# Tolerances are |a - e| <= atol + rtol * |e|.
#   analytic / EMD assignment : atol = 1e-10, rtol = 0
#     Exact simplex / permutation identities should match in double.
#   Sinkhorn vs POT           : atol = 1e-6,  rtol = 1e-5
#     Iterative scaling; POT 0.9.6 and rfugw stop on different residuals.
#   scaling vs log            : atol = 1e-7,  rtol = 1e-6
#     Same kernel, different numerics at moderate epsilon.
#   metamorphic (perm/scale)  : atol = 1e-8,  rtol = 1e-7
#
# Conservative vs fast compile flags are different artifacts
# (`inst/build-profiles.md`). This suite certifies the default portable
# build. `RFUGW_FAST_FLAGS` is not compared in-process.

.expect_close <- function(actual, expected, atol = 1e-8, rtol = 1e-6, label = NULL) {
  actual <- as.numeric(actual)
  expected <- as.numeric(expected)
  ok <- abs(actual - expected) <= atol + rtol * abs(expected)
  expect_true(
    all(ok),
    info = paste0(
      if (is.null(label)) "values" else label,
      " max|a-e|=", format(max(abs(actual - expected)), digits = 4),
      " atol=", atol, " rtol=", rtol
    )
  )
}

test_that("fixture provenance records POT version and seed", {
  fx <- read_fixture("linear_ot_fixture.json")
  expect_identical(fx$meta$pot_version, "0.9.6.post1")
  expect_equal(fx$meta$seed, 20260816)
  expect_true(nzchar(fx$meta$script))
})

test_that("analytic 2x2 EMD matches the closed-form assignment", {
  M <- matrix(c(0, 2, 2, 0), 2, 2)
  p <- c(0.4, 0.6)
  q <- c(0.4, 0.6)
  out <- ot_emd(M, p, q)
  expect_equal(out$ot_dist, 0, tolerance = 1e-12)
  expect_equal(diag(out$plan), p, tolerance = 1e-10)
  ot_validate_plan(out, p, q)
  expect_equal(ot_linear_cost(M, out), out$ot_dist, tolerance = 1e-12)
})

test_that("ot_emd and ot_sinkhorn match the POT linear-OT fixture", {
  fx <- read_fixture("linear_ot_fixture.json")
  M <- fx$balanced$inputs$M
  p <- fx$balanced$inputs$p
  q <- fx$balanced$inputs$q
  pr <- fx$balanced$params

  sink <- ot_sinkhorn(M, p, q, epsilon = pr$epsilon, method = "scaling",
                      max_iter = as.integer(pr$max_iter), tol = pr$tol)
  slog <- ot_sinkhorn(M, p, q, epsilon = pr$epsilon, method = "log",
                      max_iter = as.integer(pr$max_iter), tol = pr$tol)
  emd <- ot_emd(M, p, q)

  .expect_close(sink$ot_dist, fx$balanced$outputs$sinkhorn_cost, 1e-6, 1e-5, "sinkhorn cost")
  .expect_close(sink$plan, fx$balanced$outputs$sinkhorn_plan, 1e-6, 1e-5, "sinkhorn plan")
  .expect_close(slog$ot_dist, fx$balanced$outputs$sinkhorn_log_cost, 1e-6, 1e-5, "log cost")
  .expect_close(emd$ot_dist, fx$balanced$outputs$emd_cost, 1e-10, 0, "emd cost")
  .expect_close(emd$plan, fx$balanced$outputs$emd_plan, 1e-8, 1e-7, "emd plan")
  .expect_close(sink$ot_dist, slog$ot_dist, 1e-7, 1e-6, "scaling vs log")
})

test_that("unbalanced Sinkhorn matches the POT UOT fixture", {
  fx <- read_fixture("linear_ot_fixture.json")
  M <- fx$unbalanced$inputs$M
  p <- fx$unbalanced$inputs$p
  q <- fx$unbalanced$inputs$q
  pr <- fx$unbalanced$params
  unb <- ot_sinkhorn_unbalanced(
    M, p, q,
    epsilon = pr$epsilon,
    rho = pr$rho,
    max_iter = as.integer(pr$max_iter),
    tol = pr$tol
  )
  .expect_close(unb$ot_dist, fx$unbalanced$outputs$ot_dist, 1e-6, 1e-5, "uot cost")
  .expect_close(unb$mass, fx$unbalanced$outputs$mass, 1e-6, 1e-5, "uot mass")
  .expect_close(unb$plan, fx$unbalanced$outputs$plan, 1e-6, 1e-5, "uot plan")
  expect_lt(unb$mass, 1)
})

test_that("linear OT is permutation- and scale-equivariant", {
  set.seed(16)
  M <- matrix(runif(12), 3, 4)
  p <- c(0.2, 0.5, 0.3)
  q <- c(0.1, 0.2, 0.4, 0.3)
  perm_i <- c(2L, 1L, 3L)
  perm_j <- c(4L, 1L, 3L, 2L)
  base <- ot_sinkhorn(M, p, q, epsilon = 0.1, max_iter = 400L)
  perm <- ot_sinkhorn(M[perm_i, perm_j], p[perm_i], q[perm_j], epsilon = 0.1, max_iter = 400L)
  .expect_close(base$ot_dist, perm$ot_dist, 1e-8, 1e-7, "perm cost")
  .expect_close(base$plan[perm_i, perm_j], perm$plan, 1e-8, 1e-7, "perm plan")

  scaled <- ot_sinkhorn(3 * M, p, q, epsilon = 0.3, max_iter = 400L)
  .expect_close(scaled$ot_dist, 3 * base$ot_dist, 1e-7, 1e-6, "cost scale")
})

test_that("adversarial linear-OT inputs stay finite or fail cleanly", {
  M <- matrix(c(0, 1e6, 1e6, 0), 2, 2)
  p <- c(1e-8, 1 - 1e-8)
  q <- c(1 - 1e-8, 1e-8)
  out <- ot_sinkhorn(M, p, q, epsilon = 0.05, method = "log", max_iter = 400L)
  expect_true(is.finite(out$ot_dist))
  expect_true(all(is.finite(out$plan)))
  expect_error(ot_sinkhorn(matrix(Inf, 2, 2)), "finite")
  expect_error(ot_emd(matrix(numeric(0), 0, 2)), "nonempty|empty|positive")
  expect_error(ot_sinkhorn_unbalanced(matrix(c(0, 1, 1, 0), 2, 2), rho = 0), "rho")
  tiny <- ot_sinkhorn(matrix(c(0, 1, 1, 0), 2, 2), epsilon = 1e-4, method = "log", max_iter = 800L)
  expect_true(is.finite(tiny$ot_dist))
})

test_that("linear OT primitives are serial (no OpenMP dependence)", {
  skip_on_cran()
  M <- matrix(c(0, 1, 2, 1, 0, 1, 2, 1, 0), 3, 3)
  p <- rep(1 / 3, 3)
  q <- p
  old <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
  on.exit({
    if (is.na(old)) Sys.unsetenv("OMP_NUM_THREADS") else Sys.setenv(OMP_NUM_THREADS = old)
  })
  Sys.setenv(OMP_NUM_THREADS = "1")
  a <- ot_sinkhorn(M, p, q, epsilon = 0.1, max_iter = 300L)
  Sys.setenv(OMP_NUM_THREADS = "4")
  b <- ot_sinkhorn(M, p, q, epsilon = 0.1, max_iter = 300L)
  .expect_close(a$ot_dist, b$ot_dist, 1e-12, 0, "thread cost")
  .expect_close(a$plan, b$plan, 1e-12, 0, "thread plan")
})
