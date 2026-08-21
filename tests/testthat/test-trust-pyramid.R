.repo_text <- function(...) {
  path <- testthat::test_path("..", "..", ...)
  if (!file.exists(path)) {
    path <- trust_test_resource("numerical-trust", "workflow-contract.txt")
  }
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

test_that("PR, nightly, and release trust scopes are distinct and replayable", {
  pr <- .repo_text(".github", "workflows", "numerical-trust.yml")
  nightly <- .repo_text(".github", "workflows", "numerical-trust-nightly.yml")
  release <- .repo_text(".github", "workflows", "numerical-trust-release.yml")
  expect_match(pr, "--family=pr --scope=pr")
  expect_match(pr, "Clean source tarball and installed-package check")
  expect_match(pr, "run-mutation-proof")
  expect_match(nightly, "--family=all --scope=nightly")
  expect_match(nightly, "Extended fuzz, differential, scale, approximation")
  expect_match(release, "ubuntu-latest, macos-latest, windows-latest", fixed = TRUE)
  expect_match(release, "--family=all --scope=release --installed")
  expect_match(release, "fsanitize=address,undefined", fixed = TRUE)
  for (workflow in list(pr, nightly, release)) {
    expect_match(workflow, "collect-evidence.R")
    expect_match(workflow, "upload-artifact@v4", fixed = TRUE)
  }
})

test_that("evidence collection separates hosted and publication status", {
  collector <- .repo_text("tools", "numerical-trust", "collect-evidence.R")
  gate <- .repo_text("tools", "release-gate.R")
  expect_match(collector, "evidence_channel")
  expect_match(collector, "publication_status")
  expect_match(collector, "not_evaluated_by_numerical_trust_job")
  expect_match(collector, "representative_certificates")
  expect_match(collector, "sha256")
  expect_match(collector, "exact_commit_evidence")
  expect_match(collector, "experimental_boundaries")
  expect_match(collector, "release-dossier")
  expect_match(gate, "run-mutation-proof.R", fixed = TRUE)
  expect_match(gate, "RFUGW_TRUST_SCOPE = \"release\"")
})

test_that("performance workflows run correctness gates first", {
  fast <- .repo_text(
    ".github", "workflows", "sparse-sampled-perf-gate.yml"
  )
  nightly <- .repo_text(
    ".github", "workflows", "sparse-sampled-perf-gate-nightly.yml"
  )
  for (workflow in list(fast, nightly)) {
    accuracy <- regexpr("Accuracy gate", workflow, fixed = TRUE)[[1L]]
    benchmark <- regexpr("Run sparse sampled benchmark", workflow, fixed = TRUE)[[1L]]
    expect_gt(accuracy, 0L)
    expect_gt(benchmark, accuracy)
  }
})
