test_that("installed live-POT contract names its authority and exceptions", {
  contract <- system.file("pot-oracle-contract.md", package = "rfugw")
  expect_true(nzchar(contract))
  text <- paste(readLines(contract, warn = FALSE), collapse = "\n")

  expect_match(text, "live differential\\s+oracle, not as its sole source")
  expect_match(text, "0\\.9\\.7\\.post1")
  expect_match(text, "Entropic partial FGW gradient")
  expect_match(text, "Translation-invariant KL-UOT")
  expect_match(text, "Sampled GW")
  expect_match(text, "Entropic barycenter return semantics")
  expect_match(text, "canonical source artifact")
})

test_that("source checkout keeps baseline and latest POT policies separate", {
  source_root <- testthat::test_path("..", "..")
  baseline <- file.path(source_root, "tools/pot-oracle/requirements-baseline.txt")
  latest <- file.path(source_root, "tools/pot-oracle/requirements-latest.txt")
  skip_if_not(file.exists(baseline))
  skip_if_not(file.exists(latest))

  baseline_lines <- readLines(baseline, warn = FALSE)
  latest_lines <- readLines(latest, warn = FALSE)
  expect_true("POT==0.9.7.post1" %in% baseline_lines)
  expect_true(any(grepl("^numpy==", baseline_lines)))
  expect_true(any(grepl("^scipy==", baseline_lines)))
  expect_true(any(grepl("^POT>=0.9.7.post1,<1.0$", latest_lines)))
  expect_false(identical(baseline_lines, latest_lines))
})

test_that("an explicitly supplied live POT receipt passes the installed comparator", {
  receipt <- Sys.getenv("RFUGW_POT_ORACLE_JSON", unset = "")
  skip_if(!nzchar(receipt), "set RFUGW_POT_ORACLE_JSON to exercise the live comparator")
  skip_if_not(file.exists(receipt))
  comparator <- file.path(
    testthat::test_path("..", ".."), "tools/pot-oracle/compare.R"
  )
  skip_if_not(file.exists(comparator))

  output <- tempfile(fileext = ".json")
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      comparator,
      paste0("--input=", normalizePath(receipt, mustWork = TRUE)),
      paste0("--output=", output),
      "--scope=testthat"
    )
  )
  expect_identical(status, 0L)
  expect_true(file.exists(output))
})
