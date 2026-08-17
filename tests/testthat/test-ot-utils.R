independent_gw_n4 <- function(C1, C2, G) {
  ns <- nrow(C1)
  nt <- nrow(C2)
  s <- 0
  for (i in seq_len(ns)) {
    for (j in seq_len(ns)) {
      for (k in seq_len(nt)) {
        for (l in seq_len(nt)) {
          s <- s + (C1[i, j] - C2[k, l])^2 * G[i, k] * G[j, l]
        }
      }
    }
  }
  s
}

test_that("ot_validate_plan accepts results and rejects bad plans", {
  M <- matrix(c(0, 1, 1, 0), 2, 2)
  p <- c(0.4, 0.6)
  q <- c(0.55, 0.45)
  out <- ot_emd(M, p, q)
  expect_identical(ot_validate_plan(out, p, q, marginals = "balanced"), out$plan)
  expect_error(ot_validate_plan(out$plan + 0.2, p, q), "balanced")
  expect_error(ot_validate_plan(matrix(-1, 2, 2), p, q), "nonnegative")
  expect_error(ot_validate_plan(matrix(NA_real_, 2, 2), p, q), "finite")
  partial <- out$plan * 0.5
  expect_identical(
    ot_validate_plan(partial, p, q, mass = sum(partial), marginals = "partial"),
    partial
  )
  expect_error(
    ot_validate_plan(out$plan, p, q, mass = 0.1, marginals = "partial"),
    "mass"
  )
})

test_that("objective helpers match analytic and permutation identities", {
  C1 <- matrix(c(0, 1, 1, 0), 2, 2)
  C2 <- C1
  G <- diag(c(0.5, 0.5))
  M <- matrix(c(0, 2, 2, 0), 2, 2)
  expect_equal(ot_linear_cost(M, G), 0)
  expect_equal(ot_gw_square(C1, C2, G), independent_gw_n4(C1, C2, G), tolerance = 1e-12)
  expect_equal(ot_fgw_square(M, C1, C2, G, alpha = 0.25), 0.25 * ot_gw_square(C1, C2, G))

  perm <- c(2L, 1L)
  expect_equal(
    ot_gw_square(C1[perm, perm], C2, G[perm, ]),
    ot_gw_square(C1, C2, G),
    tolerance = 1e-12
  )
  expect_equal(ot_linear_cost(2 * M, G), 2 * ot_linear_cost(M, G))
  expect_equal(ot_entropy(G), 2 * 0.5 * log(0.5), tolerance = 1e-12)
  expect_equal(ot_kl(G, c(0.5, 0.5), c(0.5, 0.5)), ot_entropy(G) - log(0.25), tolerance = 1e-12)
})

test_that("barycentric projection handles orientation and zero mass", {
  plan <- matrix(c(1, 0, 0, 0), 2, 2)
  tgt <- matrix(c(10, 0, 0, 20), 2, 2, byrow = TRUE)
  src <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
  fwd <- ot_barycentric_project(plan, tgt, "source_to_target")
  expect_equal(fwd[1, ], c(10, 0))
  expect_true(all(is.nan(fwd[2, ])))
  fwd0 <- ot_barycentric_project(plan, tgt, "source_to_target", zero_mass = "zero")
  expect_equal(fwd0[2, ], c(0, 0))
  back <- ot_barycentric_project(plan, src, "target_to_source")
  expect_equal(back[1, ], c(1, 2))
  expect_true(all(is.nan(back[2, ])))
})

test_that("helpers accept solver results and match reported objectives", {
  set.seed(2)
  X <- matrix(rnorm(8), 4, 2)
  Y <- matrix(rnorm(8), 4, 2)
  C1 <- as.matrix(dist(X)); C1 <- C1 / max(C1)
  C2 <- as.matrix(dist(Y)); C2 <- C2 / max(C2)
  M <- as.matrix(dist(rbind(X, Y)))[1:4, 5:8]
  M <- M / max(M)
  exact <- fgw_exact_cg(M, C1, C2, alpha = 0.4, max_iter = 30L)
  expect_equal(ot_fgw_square(M, C1, C2, exact, alpha = 0.4), exact$fgw_dist, tolerance = 1e-7)
  sink <- ot_sinkhorn(M, epsilon = 0.1, max_iter = 200L)
  expect_equal(ot_linear_cost(M, sink), sink$ot_dist, tolerance = 1e-8)
  ot_validate_plan(sink, rowSums(sink$plan), colSums(sink$plan))
})
