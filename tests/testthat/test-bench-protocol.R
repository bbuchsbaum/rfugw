source(testthat::test_path("..", "..", "inst", "bench", "protocol.R"))

test_that("quality gate rejects a failed status and accepts a valid plan", {
  M <- matrix(c(0, 1, 1, 0), 2, 2)
  p <- c(0.5, 0.5)
  q <- p
  good <- ot_sinkhorn(M, p, q, epsilon = 0.1, max_iter = 400L)
  prob <- list(M = M, p = p, q = q)
  ok <- bench_check_quality("ot_sinkhorn", good, prob)
  expect_true(ok$valid)

  bad <- good
  bad$status <- "numerical_failure"
  bad$residual <- 1
  no <- bench_check_quality("ot_sinkhorn", bad, prob)
  expect_false(no$valid)
  expect_match(no$reject_reason, "status")
})

test_that("invalid-quality timings cannot be treated as a baseline row", {
  M <- matrix(c(0, 1, 1, 0), 2, 2)
  timing <- bench_run_split(
    prepare_fn = function() list(M = M, p = c(0.5, 0.5), q = c(0.5, 0.5)),
    solve_fn = function(d) {
      out <- ot_sinkhorn(d$M, d$p, d$q, epsilon = 0.1, max_iter = 1L, tol = 1e-20)
      out$status <- "numerical_failure"
      out
    },
    warmup = 1L,
    reps = 1L,
    method = "ot_sinkhorn"
  )
  expect_false(timing$valid)
  expect_true(is.na(timing$solve_ms))
})

test_that("split timing records prepare, solve, and e2e separately", {
  M <- matrix(c(0, 1, 1, 0), 2, 2)
  timing <- bench_run_split(
    prepare_fn = function() {
      Sys.sleep(0.01)
      list(M = M, p = c(0.5, 0.5), q = c(0.5, 0.5))
    },
    solve_fn = function(d) ot_sinkhorn(d$M, d$p, d$q, epsilon = 0.1, max_iter = 200L),
    warmup = 1L,
    reps = 1L,
    method = "ot_sinkhorn"
  )
  expect_true(timing$valid)
  expect_gt(timing$prepare_ms, 5)
  expect_true(is.finite(timing$solve_ms))
  expect_gte(timing$solve_ms, 0)
  expect_equal(timing$e2e_ms, timing$prepare_ms + timing$solve_ms, tolerance = 1e-8)
})

test_that("environment capture records commit, seed, threads, and hardware", {
  meta <- bench_capture_env(seed = 7L, threads = 2L, warmup = 1L, reps = 3L)
  expect_identical(meta$package, "rfugw")
  expect_true(nzchar(meta$version))
  expect_identical(meta$seed, 7L)
  expect_identical(meta$threads, 2L)
  expect_true(nzchar(meta$sysname))
  expect_true(nzchar(meta$machine))
  expect_true(nzchar(meta$timestamp))
})

test_that("partial quality uses independent reconstruction and mass bounds", {
  skip_if_not_installed("lpSolve")
  d <- bench_make_problem("fgw", 6L, 20260816L)
  out <- partial_gromov_wasserstein(
    d$C1, d$C2, p = d$p, q = d$q, m = 0.7, numItermax = 40L, tol = 1e-6, log = TRUE
  )
  ok <- bench_check_quality("partial_gromov_wasserstein", out, d)
  expect_true(ok$valid)
})

test_that("UCOOT quality accepts a finite independent-mode solve", {
  d <- bench_make_problem("ucoot", 6L, 20260816L)
  out <- unbalanced_co_optimal_transport(
    d$X, d$Y, max_iter = 20L, tol = 1e-6, max_iter_ot = 80L, log = TRUE
  )
  ok <- bench_check_quality("unbalanced_co_optimal_transport", out, d)
  expect_true(ok$valid)
})

test_that("semirelaxed quality uses independent GW reconstruction", {
  d <- bench_make_problem("fgw", 8L, 20260816L)
  out <- semirelaxed_gromov_wasserstein(
    d$C1, d$C2, p = d$p, max_iter = 80L, tol_rel = 1e-6, tol_abs = 1e-6
  )
  ok <- bench_check_quality("semirelaxed_gromov_wasserstein", out, d)
  expect_true(ok$valid)
})

test_that("shared problems are determined by kind, n, and seed", {
  a <- bench_make_problem("linear", 5L, 11L)
  b <- bench_make_problem("linear", 5L, 11L)
  expect_equal(a$M, b$M)
  c <- bench_make_problem("linear", 5L, 12L)
  expect_false(isTRUE(all.equal(a$M, c$M)))
})

test_that("current baselines are archived into scratch, not mixed in", {
  root <- tempfile("bench-results-")
  dir.create(file.path(root, "current"), recursive = TRUE)
  writeLines("old", file.path(root, "current", "runs.csv"))
  bench_archive_current(root)
  expect_false(dir.exists(file.path(root, "current")))
  archived <- list.files(file.path(root, "scratch"), recursive = TRUE)
  expect_true(any(grepl("runs.csv$", archived)))
})
