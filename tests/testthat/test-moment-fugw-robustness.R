test_that("Moment-FUGW robustness contract names every AC8 property", {
  skip_if_not_installed("jsonlite")
  contract_path <- bench_test_resource(
    "moment-fugw-robustness-contract.json"
  )
  contract <- jsonlite::fromJSON(contract_path, simplifyVector = FALSE)
  expected_properties <- c(
    "independent_node_permutations",
    "coordinate_translation",
    "orthogonal_transforms",
    "factor_gauge_changes",
    "zero_column_factor_padding",
    "cost_composition",
    "block_size_invariance",
    "split_weight_duplicates",
    "serialization",
    "deterministic_thread_settings",
    "small_epsilon",
    "extreme_rho_epsilon",
    "constant_correlation_profiles",
    "extreme_coordinate_scales",
    "underflow_diagnostics"
  )
  required_cases <- unique(unlist(lapply(
    contract$test_files,
    function(entry) unlist(entry$required_cases, use.names = FALSE)
  ), use.names = FALSE))

  expect_identical(as.integer(contract$schema_version), 1L)
  expect_identical(as.character(contract$method_family), "moment_fugw")
  expect_identical(as.character(contract$profile), "release")
  expect_identical(as.character(contract$required_package_mode), "installed")
  expect_setequal(
    unlist(contract$required_platforms, use.names = FALSE),
    c("Darwin", "Linux", "Windows")
  )
  expect_setequal(names(contract$property_cases), expected_properties)
  expect_length(required_cases, 7L)
  expect_true(all(unlist(contract$property_cases, use.names = FALSE) %in%
    required_cases))
  expect_true(all(unlist(contract$thread_environment, use.names = FALSE) ==
    "1"))
})

test_that("local robustness replay emits recomputable fail-closed evidence", {
  skip_if_not_installed("devtools")
  skip_if_not_installed("digest")
  skip_if_not_installed("jsonlite")
  source_runner <- testthat::test_path(
    "..", "..", "inst", "bench", "run_moment_fugw_robustness.R"
  )
  runner <- if (file.exists(source_runner)) {
    source_runner
  } else {
    bench_test_resource("run_moment_fugw_robustness.R")
  }
  if (!grepl("[/\\\\]inst[/\\\\]bench[/\\\\]", runner)) {
    skip("source-tree robustness runner is unavailable")
  }
  output <- tempfile("moment-fugw-robustness-")
  on.exit(unlink(output, recursive = TRUE, force = TRUE), add = TRUE)
  command_output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(runner, "local", output, "source"),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(command_output, "status", exact = TRUE)
  if (is.null(status)) status <- 0L
  expect_identical(as.integer(status), 0L, info = paste(command_output))

  summary_path <- file.path(output, "summary.json")
  cases_path <- file.path(output, "cases.csv")
  expect_true(file.exists(summary_path))
  expect_true(file.exists(cases_path))
  summary <- jsonlite::fromJSON(summary_path, simplifyVector = TRUE)
  cases <- utils::read.csv(cases_path, stringsAsFactors = FALSE)
  expect_identical(summary$status, "passed")
  expect_false(summary$admission_evaluable)
  expect_identical(summary$profile, "local")
  expect_identical(summary$package_mode, "source")
  expect_identical(summary$global_optimality, "not_claimed")
  expect_identical(summary$tests$observed_cases, nrow(cases))
  expect_identical(summary$tests$required_cases, 7L)
  expect_true(summary$tests$all_required_cases_present)
  expect_identical(summary$tests$failures, 0L)
  expect_identical(summary$tests$errors, 0L)
  expect_identical(summary$tests$warnings, 0L)
  expect_identical(summary$tests$skips, 0L)
  expect_true(all(unlist(summary$properties, use.names = FALSE)))
  expect_length(summary$properties, 15L)
  expect_true(all(cases$passed_expectations > 0L))
  expect_true(all(
    cases$failures + cases$errors + cases$warnings + cases$skips == 0L
  ))
})
