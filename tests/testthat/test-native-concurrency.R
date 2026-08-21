test_that("pure native probe is deterministic across requested thread counts", {
  runs <- lapply(c(1L, 2L, 4L), function(threads) {
    rfugw:::cpp_thread_kernel_probe(
      n_jobs = 24L,
      n_threads = threads,
      nested_threads = 4L,
      inject_failure = FALSE
    )
  })

  expect_equal(runs[[1]]$checksums, runs[[2]]$checksums, tolerance = 1e-12)
  expect_equal(runs[[1]]$checksums, runs[[3]]$checksums, tolerance = 1e-12)
  expect_equal(runs[[1]]$requested_threads, 1L)
  expect_equal(runs[[2]]$requested_threads, 2L)
  expect_equal(runs[[3]]$requested_threads, 4L)
  expect_lte(runs[[3]]$used_threads, runs[[3]]$max_threads)
  if (runs[[3]]$used_threads > 1L) {
    expect_true(all(as.integer(runs[[3]]$nested_suppressed) == 1L))
  }
})

test_that("worker exceptions are isolated and reported after the region", {
  out <- rfugw:::cpp_thread_kernel_probe(
    n_jobs = 12L,
    n_threads = 4L,
    nested_threads = 2L,
    inject_failure = TRUE
  )
  failed <- as.integer(out$failed)
  expect_equal(sum(failed), 1L)
  expect_equal(which(failed == 1L), 7L)
  expect_identical(out$errors[[7]], "injected worker failure")
  expect_true(all(is.finite(out$checksums[-7])))
})

test_that("feature-cost workers match at 1, 2, and 4 requested threads", {
  set.seed(20260820L)
  F1 <- list(
    matrix(rnorm(36L), 12L, 3L),
    matrix(rnorm(30L), 10L, 3L),
    matrix(rnorm(24L), 8L, 3L),
    matrix(rnorm(18L), 6L, 3L)
  )
  F2 <- matrix(rnorm(27L), 9L, 3L)
  run <- function(threads) {
    rfugw:::cpp_feature_cost_batch(F1, F2, TRUE, TRUE, threads)
  }
  one <- run(1L)
  two <- run(2L)
  four <- run(4L)

  expect_equal(one, two, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(one, four, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(attr(one, "requested_threads"), 1L)
  expect_equal(attr(two, "requested_threads"), 2L)
  expect_equal(attr(four, "requested_threads"), 4L)
  expect_false(attr(four, "nested_parallel"))
})

test_that("OpenMP batch regions contain no R or Rcpp access", {
  source <- readLines(
    file.path(rfugw_test_source_root(), "src", "fgw_core.cpp"),
    warn = FALSE
  )
  pragmas <- grep("#pragma omp parallel for", source, fixed = TRUE)
  batch_pragmas <- pragmas[pragmas > 4000L]
  expect_length(batch_pragmas, 4L)

  worker_lines <- unlist(lapply(batch_pragmas, function(start) {
    finish <- start + which(grepl("^#else$", source[(start + 1L):length(source)]))[1L]
    source[start:finish]
  }), use.names = FALSE)
  expect_false(any(grepl("Rcpp::|Rf_", worker_lines)))
  expect_true(any(grepl("run_worker_guarded", worker_lines, fixed = TRUE)))
  expect_true(any(grepl("_owned", worker_lines, fixed = TRUE)))
})
