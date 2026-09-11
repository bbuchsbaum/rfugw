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
  expect_match(release, "LD_PRELOAD", fixed = TRUE)
  expect_match(release, "gcc -print-file-name=libasan.so", fixed = TRUE)
  expect_match(release, "build-artifact:", fixed = TRUE)
  expect_match(release, "release-artifact.R --mode=manifest", fixed = TRUE)
  expect_match(release, "release-artifact.R --mode=verify", fixed = TRUE)
  expect_match(release, "needs: build-artifact", fixed = TRUE)
  expect_match(release, "R CMD INSTALL rfugw_*.tar.gz", fixed = TRUE)
  expect_match(release, "run_moment_fugw_robustness.R", fixed = TRUE)
  expect_match(release, ".gate/moment-fugw-robustness", fixed = TRUE)
  expect_match(release, ".gate/moment-fugw-core", fixed = TRUE)
  expect_match(
    release, "run_moment_fugw_core_certification.R", fixed = TRUE
  )
  expect_match(
    nightly, "run_moment_fugw_core_certification.R", fixed = TRUE
  )
  expect_match(release, "moment-fugw-admission:", fixed = TRUE)
  expect_match(release, "needs: [build-artifact, release-matrix]", fixed = TRUE)
  expect_match(release, "verify-moment-fugw-admission.R", fixed = TRUE)
  expect_match(release, "--mode=candidate", fixed = TRUE)
  expect_match(release, "pattern: numerical-trust-release-*", fixed = TRUE)
  expect_match(release, "merge-multiple: false", fixed = TRUE)
  expect_match(release, "include-hidden-files: true", fixed = TRUE)
  expect_match(release, "candidate-receipt.json", fixed = TRUE)
  expect_false(grepl("R CMD INSTALL .\\n", release))
  for (workflow in list(pr, nightly, release)) {
    expect_match(workflow, "collect-evidence.R")
    expect_match(workflow, "upload-artifact@v4", fixed = TRUE)
    expect_match(workflow, "include-hidden-files: true", fixed = TRUE)
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
  expect_match(collector, "git_provenance_complete")
  expect_match(collector, "commit_matches_host")
  expect_match(collector, "sys_info")
  expect_match(collector, "experimental_boundaries")
  expect_match(collector, "release-dossier")
  expect_match(gate, "run-mutation-proof.R", fixed = TRUE)
  expect_match(gate, "RFUGW_TRUST_SCOPE = \"release\"")
})

test_that("sanitizer CI preloads ASan and keeps package warnings fatal", {
  sanitizer <- .repo_text(".github", "workflows", "sanitizer.yml")
  expect_match(sanitizer, "fsanitize=address,undefined", fixed = TRUE)
  expect_match(sanitizer, "LD_PRELOAD", fixed = TRUE)
  expect_match(sanitizer, "gcc -print-file-name=libasan.so", fixed = TRUE)
  expect_match(sanitizer, "-Wall -Wextra -Werror", fixed = TRUE)
  expect_match(sanitizer, "-Wno-cast-function-type", fixed = TRUE)
})

test_that("canonical release artifact verification kills digest and commit drift", {
  script <- testthat::test_path(
    "..", "..", "tools", "numerical-trust", "release-artifact.R"
  )
  skip_if_not(file.exists(script), "release tooling is excluded from tarballs")
  scratch <- tempfile("rfugw-artifact-proof-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE, force = TRUE), add = TRUE)
  artifact <- file.path(scratch, "rfugw_0.1.0.tar.gz")
  manifest <- file.path(scratch, "release-artifact.json")
  writeBin(charToRaw("canonical artifact bytes"), artifact)
  rscript <- file.path(R.home("bin"), "Rscript")
  run <- function(extra) {
    suppressWarnings(system2(
      rscript,
      c(script, extra),
      stdout = FALSE, stderr = FALSE
    ))
  }
  commit <- system2(
    "git", c("-C", testthat::test_path("..", ".."), "rev-parse", "HEAD"),
    stdout = TRUE
  )[[1L]]
  wrong_commit <- paste(rep("b", 40L), collapse = "")

  expect_identical(run(c(
    "--mode=manifest", paste0("--artifact=", artifact),
    paste0("--manifest=", manifest), paste0("--commit=", commit),
    "--allow-dirty"
  )), 0L)
  expect_identical(run(c(
    "--mode=verify", paste0("--artifact=", artifact),
    paste0("--manifest=", manifest), paste0("--commit=", commit),
    "--allow-dirty"
  )), 0L)
  expect_false(identical(run(c(
    "--mode=verify", paste0("--artifact=", artifact),
    paste0("--manifest=", manifest), paste0("--commit=", wrong_commit),
    "--allow-dirty"
  )), 0L))
  expect_false(identical(run(c(
    "--mode=verify", paste0("--artifact=", artifact),
    paste0("--manifest=", manifest), "--commit=abc123", "--allow-dirty"
  )), 0L))

  writeBin(charToRaw("tampered artifact bytes"), artifact)
  expect_false(identical(run(c(
    "--mode=verify", paste0("--artifact=", artifact),
    paste0("--manifest=", manifest), paste0("--commit=", commit),
    "--allow-dirty"
  )), 0L))
})

test_that("Git provenance cannot infer clean state from failed commands", {
  source(
    trust_test_resource("numerical-trust", "provenance-lib.R"),
    local = TRUE
  )
  outside <- tempfile("rfugw-no-git-")
  dir.create(outside)
  provenance <- rfugw_git_provenance(outside)

  expect_true(is.na(provenance$commit))
  expect_true(is.na(provenance$working_tree_dirty))
  expect_false(provenance$git_provenance_complete)
  expect_false(provenance$exact_commit_evidence)
  expect_true(rfugw_valid_commit(paste(rep("a", 40L), collapse = "")))
  expect_false(rfugw_valid_commit("abc123"))
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
