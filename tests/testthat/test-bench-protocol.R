source(bench_test_resource("protocol.R"))

test_that("benchmark resources resolve from an installed-package-safe path", {
  protocol <- bench_test_resource("protocol.R")
  gate <- bench_test_resource("gate_protocol.R")
  caps <- bench_test_resource("ci_time_caps.json")

  expect_true(all(file.exists(protocol, gate, caps)))
  expect_match(normalizePath(protocol), "bench[/\\\\]protocol[.]R$")
})

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
  expect_true(timing$certified)
  expect_true(timing$comparison_eligible)
  expect_false(timing$performance_regression_eligible)
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

test_that("a finite but uncertified UCOOT solve is not comparison eligible", {
  d <- bench_make_problem("ucoot", 6L, 20260816L)
  out <- unbalanced_co_optimal_transport(
    d$X, d$Y, max_iter = 20L, tol = 1e-6, max_iter_ot = 80L, log = TRUE
  )
  ok <- bench_check_quality("unbalanced_co_optimal_transport", out, d)
  expect_false(ok$valid)
  expect_false(ok$certified)
  expect_false(ok$comparison_eligible)
  expect_match(ok$reject_reason, "status|not_converged|certificate")
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

test_that("protocol suite filter accepts known families and rejects unknown", {
  old <- Sys.getenv("RFUGW_PROTOCOL_SUITES", unset = NA)
  on.exit({
    if (is.na(old)) Sys.unsetenv("RFUGW_PROTOCOL_SUITES") else Sys.setenv(RFUGW_PROTOCOL_SUITES = old)
  })
  Sys.unsetenv("RFUGW_PROTOCOL_SUITES")
  expect_equal(
    bench_parse_suites("fgw,fugw"),
    c("fgw", "fugw")
  )
  expect_error(bench_parse_suites("nope"), "Unknown protocol suite")
  expect_true("sampled" %in% bench_parse_suites(""))
})

test_that("FUGW and sampled fixed budgets never masquerade as certification", {
  d <- bench_make_problem("fgw", 6L, 20260816L)
  fugw <- fugw_kl(
    Cx = d$C1, Cy = d$C2, wx = d$p, wy = d$q, M = d$M,
    epsilon = 0.05, max_iter = 10L, max_iter_ot = 120L
  )
  fugw_cert <- bench_check_quality("fugw_kl", fugw, d)
  fugw_perf <- bench_check_quality(
    "fugw_kl", fugw, d, evidence_class = "fixed_budget_performance"
  )
  expect_false(fugw_cert$valid)
  expect_true(fugw_perf$valid)
  expect_false(fugw_perf$certified)
  expect_false(fugw_perf$comparison_eligible)
  expect_true(fugw_perf$performance_regression_eligible)

  samp <- sampled_gromov_wasserstein(
    d$C1, d$C2, p = d$p, q = d$q,
    nb_samples_grad = c(3L, 2L),
    epsilon = 0.15,
    max_iter = 15L,
    random_state = 20260816L,
    log = TRUE
  )
  sampled_cert <- bench_check_quality("sampled_gromov_wasserstein", samp, d)
  sampled_perf <- bench_check_quality(
    "sampled_gromov_wasserstein", samp, d,
    evidence_class = "fixed_budget_performance"
  )
  expect_false(sampled_cert$valid)
  expect_true(sampled_perf$valid)
  expect_false(sampled_perf$certified)
  expect_false(sampled_perf$comparison_eligible)
})

test_that("nightly boundary configurations satisfy their declared evidence class", {
  seed <- 20260816L

  fugw_stop <- bench_read_thresholds("fugw_kl")$performance_stop
  for (n in c(8L, 12L)) {
    d <- bench_make_problem("fgw", n, seed)
    out <- fugw_kl(
      Cx = d$C1, Cy = d$C2, wx = d$p, wy = d$q, M = d$M,
      alpha = fugw_stop$alpha, epsilon = fugw_stop$epsilon,
      reg_marginals = unlist(fugw_stop$reg_marginals),
      max_iter = as.integer(fugw_stop$max_iter), tol = fugw_stop$tol,
      max_iter_ot = as.integer(fugw_stop$max_iter_ot),
      tol_ot = fugw_stop$tol_ot
    )
    quality <- bench_check_quality(
      "fugw_kl", out, d, evidence_class = "fixed_budget_performance"
    )
    expect_true(quality$performance_regression_eligible, info = paste("FUGW n", n))
    expect_false(quality$certified, info = paste("FUGW n", n))
  }

  sr_stop <- bench_read_thresholds("entropic_semirelaxed_gromov_wasserstein")$stop
  d <- bench_make_problem("fgw", 12L, seed)
  sr <- entropic_semirelaxed_gromov_wasserstein(
    d$C1, d$C2, p = d$p, epsilon = sr_stop$epsilon,
    max_iter = as.integer(sr_stop$max_iter), tol = sr_stop$tol,
    check_every = as.integer(sr_stop$check_every)
  )
  expect_true(
    bench_check_quality("entropic_semirelaxed_gromov_wasserstein", sr, d)$comparison_eligible
  )

  for (n in c(8L, 12L)) {
    d <- bench_make_problem("ucoot", n, seed)
    for (method in c(
      "unbalanced_co_optimal_transport",
      "fused_unbalanced_across_spaces_divergence"
    )) {
      stop <- bench_read_thresholds(method)$stop
      args <- list(
        X = d$X, Y = d$Y,
        wx_samp = d$wx_samp, wy_samp = d$wy_samp,
        wx_feat = d$wx_feat, wy_feat = d$wy_feat,
        reg_marginals = unlist(stop$reg_marginals),
        epsilon = unlist(stop$epsilon),
        max_iter = as.integer(stop$max_iter), tol = stop$tol,
        max_iter_ot = as.integer(stop$max_iter_ot), tol_ot = stop$tol_ot,
        log = TRUE
      )
      if (identical(method, "fused_unbalanced_across_spaces_divergence")) {
        args$reg_type <- "joint"
      }
      out <- do.call(get(method, mode = "function"), args)
      quality <- bench_check_quality(method, out, d)
      expect_true(quality$comparison_eligible, info = paste(method, "n", n))
      expect_true(out$inner_converged, info = paste(method, "n", n))
      expect_true(
        out$max_inner_residual <= stop$tol_ot,
        info = paste(method, "n", n)
      )
    }
  }
})

test_that("threshold changes require a retained matching evidence entry", {
  expect_no_error(bench_validate_threshold_history())
  history <- jsonlite::fromJSON(bench_threshold_history_path(), simplifyVector = FALSE)
  expect_identical(
    history$entries[[length(history$entries)]]$thresholds_md5,
    unname(tools::md5sum(bench_thresholds_path())[[1]])
  )
})

test_that("sampled protocol records answer quality against budget", {
  script <- readLines(bench_test_resource("run_protocol.R"), warn = FALSE)
  expect_true(any(grepl("c(2L, 1L)", script, fixed = TRUE)))
  expect_true(any(grepl("c(4L, 2L)", script, fixed = TRUE)))
  expect_true(any(grepl("c(8L, 4L)", script, fixed = TRUE)))
  expect_true(any(grepl("approximation_error_abs", script, fixed = TRUE)))
  expect_true(any(grepl("reference_value", script, fixed = TRUE)))
})

test_that("protocol gate distinguishes infrastructure from solver failures", {
  gate <- bench_test_resource("gate_protocol.R")
  caps <- bench_test_resource("ci_time_caps.json")
  missing <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(gate, tempfile("no-such-gate-"), caps, "pr"),
    stdout = TRUE,
    stderr = TRUE
  ))
  expect_identical(attr(missing, "status"), 2L)

  good <- tempfile("gate-ok-")
  dir.create(good)
  jsonlite::write_json(
    list(
      commit = "test", seed = 1L, threads = 1L,
      profile = "conservative", suites = "fgw",
      thresholds_md5 = unname(tools::md5sum(bench_test_resource("thresholds.json"))[[1]]),
      threshold_policy_version = 2L
    ),
    file.path(good, "meta.json"),
    auto_unbox = TRUE
  )
  utils::write.csv(
    data.frame(
      suite = "fgw", method = "fgw_entropic", n = 12L, valid = TRUE,
      certified = TRUE, comparison_eligible = TRUE,
      performance_regression_eligible = FALSE,
      evidence_class = "certified_comparison", reject_reason = "",
      prepare_ms = 1, setup_ms = 1, solve_ms = 10, e2e_ms = 11,
      mem_solve_bytes = NA_real_, mem_e2e_bytes = NA_real_,
      requested_precision = "double", effective_precision = "double",
      requested_threads = 1L, used_threads = 1L,
      stringsAsFactors = FALSE
    ),
    file.path(good, "runs.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(
      suite = "fgw", method = "fgw_entropic", n = 12L, valid = TRUE,
      certified = TRUE, comparison_eligible = TRUE,
      performance_regression_eligible = FALSE,
      evidence_class = "certified_comparison", status = "converged",
      stringsAsFactors = FALSE
    ),
    file.path(good, "quality.csv"),
    row.names = FALSE
  )
  ok <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(gate, good, caps, "pr"),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_true(is.null(attr(ok, "status")) || identical(attr(ok, "status"), 0L))

  bad <- tempfile("gate-bad-")
  dir.create(bad)
  file.copy(file.path(good, "meta.json"), file.path(bad, "meta.json"))
  file.copy(file.path(good, "quality.csv"), file.path(bad, "quality.csv"))
  utils::write.csv(
    transform(
      utils::read.csv(file.path(good, "runs.csv"), stringsAsFactors = FALSE),
      valid = FALSE,
      certified = FALSE,
      comparison_eligible = FALSE,
      reject_reason = "status=numerical_failure"
    ),
    file.path(bad, "runs.csv"),
    row.names = FALSE
  )
  failed <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(gate, bad, caps, "pr"),
    stdout = TRUE,
    stderr = TRUE
  ))
  expect_identical(attr(failed, "status"), 1L)
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
