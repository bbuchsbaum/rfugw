source(bench_test_resource("protocol.R"))
source(bench_test_resource("benchmark_moment_fugw.R"))

.moment_benchmark_fit <- function(max_iter = 30L, max_iter_ot = 800L) {
  problem <- bench_make_moment_fugw_problem(
    8L, structure_rank = 4L, feature_dimension = 3L, seed = 90210L
  )
  fit <- bench_solve_moment_fugw(
    problem,
    max_iter = max_iter,
    max_iter_ot = max_iter_ot,
    tol = 1e-6,
    tol_ot = 1e-8,
    epsilon = 0.1,
    reg_marginals = c(1, 1),
    block_size = 4L
  )
  list(problem = problem, fit = fit)
}

test_that("Moment-FUGW certification uses structured stationary evidence", {
  fixture <- .moment_benchmark_fit()
  expect_identical(fixture$fit$status, "converged_stationary")
  quality <- bench_check_quality(
    "fugw_factorized", fixture$fit, fixture$problem,
    evidence_class = "certified_comparison"
  )
  expect_true(quality$valid)
  expect_true(quality$certified)
  expect_true(quality$comparison_eligible)
  expect_true(quality$timing_eligible)
  expect_false(quality$threshold_update_eligible)
  expect_identical(quality$contract_status, "candidate")
  expect_identical(quality$stationarity_classification, "stationary")
})

test_that("Moment-FUGW quality gate fails closed for independent failures", {
  fixture <- .moment_benchmark_fit()
  mutations <- list(
    nonfinite = function(x) {
      x$fugw_cost <- Inf
      x
    },
    nonstationary = function(x) {
      x$certificate$outer_stationarity$certified <- FALSE
      x
    },
    inner = function(x) {
      x$certificate$inner_uot$certified <- FALSE
      x
    },
    sample_block = function(x) {
      x$certificate$outer_stationarity$final_blocks$sample$certified <- FALSE
      x
    },
    geometry = function(x) {
      x$geometry_certified <- FALSE
      x$certificate$geometry$certified <- FALSE
      x
    },
    support = function(x) {
      x$certificate$support$complete <- FALSE
      x
    },
    dense = function(x) {
      x$runtime_provenance$memory_contract$dense_allocations_requested <- TRUE
      x
    }
  )
  for (name in names(mutations)) {
    bad <- mutations[[name]](fixture$fit)
    quality <- bench_check_quality(
      "fugw_factorized", bad, fixture$problem,
      evidence_class = "certified_comparison"
    )
    expect_false(quality$valid, info = name)
    expect_false(quality$timing_eligible, info = name)
    expect_false(quality$threshold_update_eligible, info = name)
    expect_true(nzchar(quality$reject_reason), info = name)
  }
})

test_that("fixed-work scaling is visible but never a timing comparison", {
  fixture <- .moment_benchmark_fit(max_iter = 1L, max_iter_ot = 1L)
  quality <- bench_check_quality(
    "fugw_factorized", fixture$fit,
    list(expected_work = list(outer_iterations = 1L)),
    evidence_class = "fixed_work_scaling"
  )
  expect_true(quality$valid)
  expect_true(quality$scaling_eligible)
  expect_false(quality$certified)
  expect_false(quality$comparison_eligible)
  expect_false(quality$timing_eligible)
  expect_false(quality$threshold_update_eligible)
  expect_match(quality$limitation_reason, "nonstationary_fixed_work")

  fixture$fit$status <- "numerical_failure"
  rejected <- bench_check_quality(
    "fugw_factorized", fixture$fit,
    list(expected_work = list(outer_iterations = 1L)),
    evidence_class = "fixed_work_scaling"
  )
  expect_false(rejected$valid)
  expect_false(rejected$scaling_eligible)
})

test_that("multiscale timing requires every hierarchy certificate", {
  fixture <- .moment_benchmark_fit()
  multi <- fixture$fit
  multi$multiscale_certified <- TRUE
  multi$hierarchy_transfer_certified <- TRUE
  multi$certificate$multiscale <- list(
    certified = TRUE,
    all_levels_stationary = TRUE,
    all_inner_uot_certified = TRUE,
    hierarchy_certified = TRUE,
    all_geometry_certified = TRUE
  )
  quality <- bench_check_quality(
    "fugw_multiscale", multi, fixture$problem,
    evidence_class = "certified_comparison"
  )
  expect_true(quality$timing_eligible)

  multi$certificate$multiscale$hierarchy_certified <- FALSE
  rejected <- bench_check_quality(
    "fugw_multiscale", multi, fixture$problem,
    evidence_class = "certified_comparison"
  )
  expect_false(rejected$valid)
  expect_match(rejected$reject_reason, "hierarchy_transfer")

  multi$certificate$multiscale$hierarchy_certified <- TRUE
  multi$status <- "converged_final_level_with_budgeted_warm_starts"
  multi$multiscale_certified <- FALSE
  multi$certificate$multiscale$certified <- FALSE
  multi$certificate$multiscale$all_levels_stationary <- FALSE
  budgeted <- bench_check_quality(
    "fugw_multiscale", multi, fixture$problem,
    evidence_class = "certified_comparison"
  )
  expect_false(budgeted$valid)
  expect_false(budgeted$timing_eligible)
  expect_match(budgeted$reject_reason, "budgeted_warm_starts|multiscale")
})

test_that("benchmark design and threshold history are fail closed", {
  design <- bench_validate_moment_fugw_design()
  expect_true(design$passed)
  expect_true(all(design$checks))
  contract <- bench_read_moment_fugw_thresholds()
  expect_equal(as.integer(contract$scaling$domain_sizes),
               c(625L, 1250L, 2500L, 5000L))
  expect_equal(as.integer(contract$scaling$structure_ranks),
               c(4L, 16L, 64L, 256L))
  expect_equal(as.integer(contract$scaling$heldout_batches),
               c(1L, 10L, 100L, 1000L))
  expect_true(contract$scaling$require_r_visible_allocation)
  expect_setequal(
    contract$scaling$r_visible_allocation_actions,
    c("solve", "objective", "apply", "adjoint")
  )
  expect_no_error(bench_validate_moment_fugw_threshold_history())
  expect_error(
    bench_validate_moment_fugw_threshold_history(TRUE),
    "candidate-only"
  )
})

test_that("R-visible allocation profiling is positive and separate", {
  profile <- bench_profile_r_visible_allocation(function() {
    matrix(runif(4096), 64L, 64L)
  })
  expect_true(is.matrix(profile$value))
  expect_gt(profile$bytes, 4096 * 8)
  expect_gt(profile$events, 0L)
})

test_that("release provenance fails closed when Git identity is unavailable", {
  outside <- tempfile("moment-fugw-no-git-")
  dir.create(outside)
  provenance <- bench_git_provenance(outside)

  expect_true(is.na(provenance$commit))
  expect_true(is.na(provenance$git_dirty))
  expect_true(is.na(provenance$git_status_entry_count))
  expect_false(provenance$git_provenance_complete)

  missing <- list(
    rfugw_fast_flags = "",
    git_provenance_complete = FALSE,
    git_dirty = NA,
    commit = NA_character_
  )
  expect_false(.bench_moment_fugw_release_hardware(
    "full", "installed", missing
  ))

  clean <- list(
    rfugw_fast_flags = "",
    git_provenance_complete = TRUE,
    git_dirty = FALSE,
    commit = paste(rep("a", 40L), collapse = "")
  )
  expect_true(.bench_moment_fugw_release_hardware(
    "full", "installed", clean
  ))
  clean$git_dirty <- TRUE
  expect_false(.bench_moment_fugw_release_hardware(
    "full", "installed", clean
  ))
})

test_that("retained kernel evidence satisfies its frozen admission gates", {
  evidence <- utils::read.csv(
    bench_test_resource("moment-fugw-kernel-evidence.csv"),
    stringsAsFactors = FALSE
  )
  endpoints <- evidence[
    evidence$evidence_group == "component_profile" &
      evidence$component == "whole_fixed_work_solver",
    , drop = FALSE
  ]
  smaller <- evidence[
    evidence$evidence_group == "fresh_process_scaling" &
      evidence$n_source %in% c(625L, 1250L),
    , drop = FALSE
  ]
  rss <- evidence[evidence$evidence_group == "fresh_process_rss", , drop = FALSE]

  expect_equal(sort(endpoints$n_source), c(1000L, 2000L))
  expect_true(all(endpoints$improvement_percent >= 15))
  expect_true(all(smaller$after_value <= 1.05 * smaller$before_value))
  expect_equal(sort(rss$n_source), c(625L, 1250L, 2500L, 5000L))
  expect_true(all(rss$after_value <= rss$before_value))
  expect_true(all(evidence$passed))
  expect_true(all(evidence$admission == "candidate_source_tree"))
  expect_false(any(grepl("release", evidence$admission, fixed = TRUE)))
})

test_that("Moment-FUGW admission ledger is fail-closed and reproducible", {
  ledger <- jsonlite::fromJSON(
    bench_test_resource("moment-fugw-admission-evidence.json"),
    simplifyVector = TRUE
  )
  expect_equal(ledger$schema_version, 1L)
  expect_identical(ledger$status, "experimental_candidate")
  expect_false(ledger$supported)
  expect_identical(ledger$global_optimality, "not_claimed")
  expect_true(ledger$candidate_binding$git_dirty)
  expect_false(ledger$candidate_binding$exact_candidate_commit_bound)
  expect_false(ledger$release_decision$eligible)
  expect_identical(ledger$release_decision$public_status, "experimental")
  expect_setequal(
    ledger$release_decision$unmet_gates,
    c(
      "exact_candidate_commit_binding",
      "clean_same_commit_installed_full_performance_artifact",
      "hosted_cross_platform_release_replay",
      "independent_numerical_performance_and_scientific_review"
    )
  )

  artifact_gates <- c(
    "independent_exact_factor_oracle", "multiscale_efficacy",
    "candidate_performance", "installed_full_candidate_performance",
    "installed_full_rprofmem_candidate_performance",
    "planted_image_validation",
    "planted_cortical_validation", "validation_protocol",
    "core_certification_protocol", "robustness_protocol"
  )
  for (gate_name in artifact_gates) {
    gate <- ledger$gates[[gate_name]]
    expect_true(gate$passed, info = gate_name)
    expect_identical(
      unname(tools::md5sum(bench_test_resource(gate$artifact))[[1L]]),
      gate$artifact_md5,
      info = gate_name
    )
  }

  accuracy_gate <- ledger$gates$independent_exact_factor_oracle
  expect_identical(
    unname(tools::md5sum(bench_test_resource(
      accuracy_gate$design_artifact
    ))[[1L]]),
    accuracy_gate$design_artifact_md5
  )
  expect_equal(accuracy_gate$zero_weight_rows, 6L)
  expect_equal(accuracy_gate$near_zero_weight_rows, 3L)

  core_protocol <- ledger$gates$core_certification_protocol
  expect_equal(core_protocol$local_accuracy_rows, 12L)
  expect_equal(core_protocol$local_multiscale_rows, 10L)
  expect_equal(core_protocol$local_required_cases, 27L)
  expect_false(core_protocol$release_receipts_generated)
  robustness_protocol <- ledger$gates$robustness_protocol
  expect_equal(robustness_protocol$named_properties, 15L)
  expect_false(robustness_protocol$release_receipts_generated)

  installed <- ledger$gates$installed_full_candidate_performance
  archive <- bench_test_resource(installed$raw_worker_archive)
  expect_identical(
    unname(tools::md5sum(archive)[[1L]]),
    installed$raw_worker_archive_md5
  )
  attestation <- jsonlite::fromJSON(
    bench_test_resource(installed$artifact), simplifyVector = TRUE
  )
  expect_identical(
    attestation$status,
    "candidate_installed_artifact_provenance_incomplete"
  )
  expect_false(attestation$release_baseline)
  expect_equal(attestation$retained_run_artifact$raw_artifact_count, 228L)
  expect_equal(attestation$retained_run_artifact$invalid_rows, 0L)
  expect_false(attestation$gate_audit$exact_commit_provenance)
  expect_false(attestation$gate_audit$admission_evaluable)
  expect_true(attestation$provenance_correction$regression_test_added)

  current <- ledger$gates$installed_full_rprofmem_candidate_performance
  current_archive <- bench_test_resource(current$raw_worker_archive)
  expect_identical(
    unname(tools::md5sum(current_archive)[[1L]]),
    current$raw_worker_archive_md5
  )
  current_attestation <- jsonlite::fromJSON(
    bench_test_resource(current$artifact), simplifyVector = TRUE
  )
  expect_identical(
    current_attestation$status,
    "candidate_installed_artifact_current_contract_provenance_incomplete"
  )
  expect_equal(current_attestation$retained_run_artifact$raw_artifact_count, 244L)
  expect_equal(current_attestation$retained_run_artifact$measured_and_warmup_rows, 60L)
  expect_true(current_attestation$gate_audit$r_visible_allocation_complete)
  expect_false(current_attestation$gate_audit$clean_exact_source_tree)
  expect_false(current_attestation$gate_audit$admission_evaluable)
  expect_false(current_attestation$release_baseline)

  receipt_protocol <- ledger$gates$admission_receipt_protocol
  expect_true(receipt_protocol$passed)
  expect_identical(
    unname(tools::md5sum(bench_test_resource(
      receipt_protocol$artifact
    ))[[1L]]),
    receipt_protocol$artifact_md5
  )
  expect_identical(
    unname(tools::md5sum(bench_test_resource(
      receipt_protocol$review_template
    ))[[1L]]),
    receipt_protocol$review_template_md5
  )
  expect_false(receipt_protocol$technical_receipt_generated)
  expect_false(receipt_protocol$promotion_receipt_generated)

  expect_identical(
    unname(tools::md5sum(bench_test_resource(
      ledger$threshold_contract$artifact
    ))[[1L]]),
    ledger$threshold_contract$artifact_md5
  )
  expect_identical(
    unname(tools::md5sum(bench_test_resource(
      ledger$threshold_contract$history_artifact
    ))[[1L]]),
    ledger$threshold_contract$history_artifact_md5
  )

  performance <- ledger$gates$candidate_performance
  exponent <- function(size, value) {
    unname(stats::coef(stats::lm(log(value) ~ log(size)))[[2L]])
  }
  incremental_rss <- performance$absolute_peak_rss_bytes -
    performance$package_load_baseline_rss_bytes
  expect_equal(
    exponent(performance$domain_sizes, incremental_rss),
    performance$incremental_peak_rss_exponent,
    tolerance = 1e-14
  )
  expect_equal(
    exponent(
      tail(performance$domain_sizes, 3L),
      tail(performance$fixed_work_seconds, 3L)
    ),
    performance$fixed_work_runtime_exponent_largest_three,
    tolerance = 1e-14
  )
  expect_lt(
    performance$incremental_peak_rss_exponent,
    performance$incremental_peak_rss_exponent_max
  )
  expect_true(
    performance$fixed_work_runtime_exponent_largest_three >=
      performance$fixed_work_runtime_exponent_interval[[1L]] &&
      performance$fixed_work_runtime_exponent_largest_three <=
      performance$fixed_work_runtime_exponent_interval[[2L]]
  )
  expect_false(performance$release_baseline)
  expect_false(ledger$gates$hosted_cross_platform_release_replay$passed)
  expect_false(ledger$gates$independent_review$passed)
})

test_that("native peak RSS parsers preserve platform units", {
  expect_equal(
    bench_parse_peak_rss(
      "  123456 maximum resident set size", sysname = "Darwin"
    ),
    123456
  )
  expect_equal(
    bench_parse_peak_rss(
      "Maximum resident set size (kbytes): 321", sysname = "Linux"
    ),
    321 * 1024
  )
  expect_true(is.na(bench_parse_peak_rss("unrelated", "Darwin")))
  expect_equal(
    bench_loglog_exponent(c(100, 200, 400, 800),
                          c(1, 4, 16, 64)),
    2,
    tolerance = 1e-12
  )
})

test_that("invalid raw timings remain visible but are not reportable", {
  fake <- list(
    process_status = 0L,
    rss_backend = "macos_time_l_bytes",
    peak_rss_bytes = 200,
    result = list(
      valid = FALSE,
      action = "solve",
      operation_seconds = 1.25,
      timing_eligible = FALSE,
      scaling_eligible = FALSE,
      reject_reason = "certificate_outer_stationarity"
    ),
    artifacts = list()
  )
  row <- .bench_worker_row(
    fake, "end_to_end_convergence", "measured", 1L, 100, "smoke"
  )
  expect_equal(row$operation_seconds_raw, 1.25)
  expect_true(is.na(row$reportable_seconds))
  expect_match(row$reject_reason, "outer_stationarity")
})

test_that("fresh workers measure RSS and retain matrix-free action evidence", {
  skip_if(bench_detect_rss_backend() %in% c("unavailable", "unsupported"))
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  source_mode <- file.exists(file.path(root, "DESCRIPTION")) &&
    dir.exists(file.path(root, "R"))
  if (source_mode) skip_if_not_installed("pkgload")
  artifacts <- tempfile("moment-worker-")
  dir.create(artifacts)
  base <- list(
    package_mode = if (source_mode) "source" else "installed",
    package_root = root,
    protocol_path = bench_test_resource("protocol.R"),
    benchmark_path = bench_test_resource("benchmark_moment_fugw.R"),
    threads = 1L,
    seed = 20260827L,
    feature_dimension = 3L,
    epsilon = 0.1,
    reg_marginals = c(1, 1),
    tol = 1e-6,
    tol_ot = 1e-8,
    block_size = 4L
  )
  baseline <- bench_run_fresh_process(
    utils::modifyList(base, list(action = "baseline")),
    artifacts, "baseline"
  )
  if (!identical(baseline$process_status, 0L) ||
      !is.finite(baseline$peak_rss_bytes)) {
    skip(paste0(
      "Fresh-process peak RSS is not operational in this sandbox; ",
      "dedicated benchmark artifacts retain the native-memory gate."
    ))
  }
  expect_identical(baseline$process_status, 0L)
  expect_true(is.finite(baseline$peak_rss_bytes))

  fit_path <- file.path(artifacts, "fit.rds")
  solve <- bench_run_fresh_process(
    utils::modifyList(base, list(
      action = "solve", n = 8L, structure_rank = 4L,
      max_iter = 1L, max_iter_ot = 1L,
      evidence_class = "fixed_work_scaling", fit_path = fit_path
    )),
    artifacts, "solve"
  )
  expect_identical(solve$process_status, 0L)
  expect_true(solve$result$valid)
  expect_true(solve$result$scaling_eligible)
  expect_true(solve$result$dense_sentinel$passed)
  expect_equal(solve$result$dense_sentinel$full_block_calls, 0L)

  objective <- bench_run_fresh_process(
    utils::modifyList(base, list(
      action = "objective", fit_path = fit_path, batch = 1L,
      source_evidence_class = "fixed_work_scaling"
    )),
    artifacts, "objective"
  )
  expect_identical(objective$process_status, 0L)
  expect_true(objective$result$valid)
  expect_true(objective$result$memory_contract_eligible)
  expect_true(objective$result$dense_sentinel$passed)
  expect_true(all(file.exists(unlist(objective$artifacts))))
})
