test_that("Moment-FUGW core contract binds AC1 through AC5 replays", {
  skip_if_not_installed("jsonlite")
  contract_path <- bench_test_resource("moment-fugw-core-contract.json")
  contract <- jsonlite::fromJSON(contract_path, simplifyVector = FALSE)
  required_cases <- unique(unlist(lapply(
    contract$test_files,
    function(entry) unlist(entry$required_cases, use.names = FALSE)
  ), use.names = FALSE))
  expected_test_files <- c(
    "tests/testthat/test-fugw-factorized-oracle.R",
    "tests/testthat/test-fugw-factorized-certificate.R",
    "tests/testthat/test-fugw-geometry-audit.R",
    "tests/testthat/test-fugw-multiscale.R"
  )

  expect_identical(as.integer(contract$schema_version), 1L)
  expect_identical(as.character(contract$method_family), "moment_fugw")
  expect_identical(as.character(contract$profile), "release")
  expect_identical(as.character(contract$required_package_mode), "installed")
  expect_setequal(
    unlist(contract$required_platforms, use.names = FALSE),
    c("Darwin", "Linux", "Windows")
  )
  expect_setequal(
    vapply(
      contract$test_files,
      function(entry) as.character(entry$path[[1L]]),
      character(1)
    ),
    expected_test_files
  )
  expect_length(required_cases, 27L)
  expect_identical(as.integer(contract$accuracy$rows), 12L)
  expect_identical(as.integer(contract$accuracy$seed), 20260827L)
  design <- utils::read.csv(
    bench_test_resource(basename(as.character(contract$accuracy$design))),
    stringsAsFactors = FALSE
  )
  expect_identical(nrow(design), 12L)
  expect_identical(design$case_id, sprintf("mfugw-%02d", 1:12))
  expect_setequal(
    paste0(design$ns, "x", design$nt),
    unlist(contract$accuracy$expected_pairs, use.names = FALSE)
  )
  expect_true(any(design$zero_weight))
  expect_true(any(design$near_zero_weight))
  expect_true(all(!design$zero_weight | !design$near_zero_weight))
  expect_true(isTRUE(contract$accuracy$require_zero_weight))
  expect_true(isTRUE(contract$accuracy$require_near_zero_weight))
  expect_identical(as.integer(contract$multiscale$rows), 10L)
  expect_equal(
    as.numeric(contract$multiscale$minimum_improved_fraction), 0.8
  )
  expect_equal(
    as.numeric(contract$multiscale$minimum_median_inner_iteration_reduction),
    0.2
  )
  expect_true(all(unlist(contract$thread_environment, use.names = FALSE) ==
    "1"))
})

test_that("retained core evidence includes positive near-zero weights", {
  accuracy <- utils::read.csv(
    bench_test_resource("moment-fugw-accuracy-evidence.csv"),
    stringsAsFactors = FALSE
  )
  multiscale <- utils::read.csv(
    bench_test_resource("moment-fugw-multiscale-evidence.csv"),
    stringsAsFactors = FALSE
  )

  expect_equal(nrow(accuracy), 12L)
  expect_true(all(accuracy$pass))
  expect_true(any(accuracy$zero_weight))
  expect_true(any(accuracy$near_zero_weight))
  expect_true(all(!accuracy$zero_weight | !accuracy$near_zero_weight))
  expect_lte(max(accuracy$objective_relative_error), 1e-8)
  expect_lte(max(accuracy$component_relative_error), 1e-8)
  expect_lte(max(accuracy$action_relative_error), 1e-6)
  expect_lte(max(accuracy$adjoint_relative_error), 1e-6)
  expect_lte(max(accuracy$mass_relative_error), 1e-7)
  expect_lte(max(accuracy$moment_relative_error), 1e-7)
  expect_equal(nrow(multiscale), 10L)
  expect_true(all(multiscale$efficacy_gate))
  expect_gte(mean(multiscale$inner_iteration_reduction > 0), 0.8)
  expect_gte(stats::median(multiscale$inner_iteration_reduction), 0.2)
  expect_lte(max(abs(multiscale$objective_relative_difference)), 1e-4)
  expect_gte(min(multiscale$heldout_score_ratio), 0.99)
  expect_lte(max(multiscale$plan_action_relative_difference), 1e-7)
})
