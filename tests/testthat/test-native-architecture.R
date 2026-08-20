test_that("canonical native contracts have R-independent module boundaries", {
  root <- rfugw_test_source_root()
  files <- file.path(root, "src", c(
    "gw_square.h", "gw_square.cpp", "transport_simplex.h",
    "approximation_cache.h", "batch_worker.h", "thread_probe.cpp"
  ))
  expect_true(all(file.exists(files)))

  gw_header <- readLines(files[[1]], warn = FALSE)
  gw_adapter <- readLines(files[[2]], warn = FALSE)
  transport <- readLines(files[[3]], warn = FALSE)
  approximation <- readLines(files[[4]], warn = FALSE)
  worker <- readLines(files[[5]], warn = FALSE)

  expect_false(any(grepl("Rcpp", gw_header, fixed = TRUE)))
  expect_false(any(grepl("Rcpp", transport, fixed = TRUE)))
  expect_false(any(grepl("Rcpp", approximation, fixed = TRUE)))
  expect_false(any(grepl("Rcpp", worker, fixed = TRUE)))
  expect_true(any(grepl("gw_square_reverse_tensor", gw_adapter, fixed = TRUE)))
  expect_true(any(grepl("duality_gap", transport, fixed = TRUE)))
  expect_true(any(grepl("template <typename T>", approximation, fixed = TRUE)))
  expect_true(any(grepl("run_worker_guarded", worker, fixed = TRUE)))
})

test_that("exact and partial paths use canonical general square-loss algebra", {
  root <- rfugw_test_source_root()
  core <- readLines(file.path(root, "src", "fgw_core.cpp"), warn = FALSE)
  canonical_calls <- grep("rfugw::gw_square_terms", core, fixed = TRUE)
  expect_gte(length(canonical_calls), 4L)
  expect_false(any(grepl("C1 * G * C2", core, fixed = TRUE)))

  optimized_general <- grep("tensor_product_asym_blas_scaled", core, fixed = TRUE)
  expect_gte(length(optimized_general), 4L)
  architecture <- readLines(file.path(root, "inst", "native-architecture.md"), warn = FALSE)
  expect_true(any(grepl("deliberate second implementation", architecture, fixed = TRUE)))
  expect_true(any(grepl("test-gw-gradient-laws.R", architecture, fixed = TRUE)))
})

test_that("canonical square-loss adapter still matches brute force", {
  set.seed(20260820L)
  n <- 3L
  m <- 4L
  C1 <- matrix(runif(n * n), n, n)
  C2 <- matrix(runif(m * m), m, m)
  G <- matrix(rexp(n * m), n, m)
  G <- G / sum(G)
  out <- rfugw:::cpp_gw_square_terms_square(C1, C2, G, symmetric = FALSE)

  brute_forward <- matrix(0, n, m)
  brute_reverse <- matrix(0, n, m)
  for (i in seq_len(n)) for (j in seq_len(m)) {
    for (k in seq_len(n)) for (l in seq_len(m)) {
      brute_forward[i, j] <- brute_forward[i, j] +
        (C1[i, k] - C2[j, l])^2 * G[k, l]
      brute_reverse[i, j] <- brute_reverse[i, j] +
        (C1[k, i] - C2[l, j])^2 * G[k, l]
    }
  }
  expect_equal(out$forward_tensor, brute_forward, tolerance = 1e-12)
  expect_equal(out$reverse_tensor, brute_reverse, tolerance = 1e-12)
  expect_equal(out$grad, brute_forward + brute_reverse, tolerance = 1e-12)
})
