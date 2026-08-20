with_restored_env <- function(values, code) {
  keys <- names(values)
  old <- Sys.getenv(keys, unset = NA_character_)
  on.exit({
    present <- !is.na(old)
    if (any(present)) do.call(Sys.setenv, as.list(old[present]))
    if (any(!present)) Sys.unsetenv(keys[!present])
  }, add = TRUE)
  do.call(Sys.setenv, as.list(values))
  force(code)
}

test_that("runtime provenance enumerates every behavior-changing environment control", {
  provenance <- rfugw:::.solver_runtime_provenance()
  expected <- c(
    "RFUGW_ENABLE_WARM_START", "RFUGW_ADAPTIVE_INNER_TOL",
    "RFUGW_ADAPTIVE_INNER_TOL_STAGE1", "RFUGW_ADAPTIVE_INNER_TOL_STAGE2",
    "RFUGW_SINKHORN_TARGET_TOL_ABS_MULT", "RFUGW_SINKHORN_TARGET_TOL_MULT",
    "RFUGW_SINKHORN_TARGET_TOL_CAP_D", "RFUGW_SINKHORN_TARGET_TOL_CAP_F",
    "RFUGW_SINKHORN_REL_TOL_MULT", "RFUGW_SINKHORN_REL_TOL_FLOOR_D",
    "RFUGW_SINKHORN_REL_TOL_FLOOR_F", "RFUGW_SINKHORN_COL_REL_TOL_MULT",
    "RFUGW_SINKHORN_COL_REL_TOL_FLOOR_D",
    "RFUGW_SINKHORN_COL_REL_TOL_FLOOR_F",
    "RFUGW_MATVEC_BLOCKED_MIN_WORK_D", "RFUGW_MATVEC_GEMV_MIN_WORK_D",
    "RFUGW_MATVEC_BLOCKED_MIN_WORK_F", "RFUGW_MATVEC_GEMV_MIN_WORK_F",
    "RFUGW_AUTOTUNE_MATVEC", "RFUGW_SAMPLED_MIXED",
    "RFUGW_SEMIRELAXED_MIXED", "RFUGW_PIN_BLAS_THREADS",
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"
  )
  expect_setequal(names(provenance$environment), expected)
  expect_true(all(c(
    "openmp_available", "openmp_max_threads", "matvec_thresholds_initialized"
  ) %in% names(provenance$native_effective)))
})

test_that("invalid environment values fall back with a recorded warning", {
  key <- "RFUGW_ADAPTIVE_INNER_TOL_STAGE1"
  before <- Sys.getenv(key, unset = NA_character_)
  inside <- with_restored_env(
    stats::setNames("not-a-number", key),
    rfugw:::.solver_runtime_provenance()
  )
  after <- Sys.getenv(key, unset = NA_character_)

  expect_false(inside$environment[[key]]$valid)
  expect_identical(
    inside$environment[[key]]$source,
    "invalid_environment_fallback"
  )
  expect_match(inside$warnings, key)
  expect_identical(after, before)

  semirelaxed_default <- with_restored_env(
    c(RFUGW_SEMIRELAXED_MIXED = "invalid"),
    rfugw:::.runtime_env_bool("RFUGW_SEMIRELAXED_MIXED", TRUE)
  )
  expect_true(semirelaxed_default)
})

test_that("results explain automatic precision and Sinkhorn transitions", {
  M <- matrix(c(0, 1, 1, 0), 2, 2)
  out <- fgw_entropic(
    M, M, M, epsilon = 0.1, max_iter = 5L,
    tol = 1e-9, sinkhorn_tol = 1e-9,
    precision = "mixed", sinkhorn_method = "auto"
  )
  automatic <- out$runtime_provenance$automatic
  expect_identical(automatic$precision$requested, "mixed")
  expect_identical(automatic$precision$effective, "strict_double")
  expect_identical(automatic$precision$source, "automatic_policy")
  expect_identical(automatic$sinkhorn_method$requested, "auto")
  expect_identical(automatic$sinkhorn_method$effective, "log")
  expect_identical(automatic$sinkhorn_method$transition, "auto_to_log")
  expect_identical(
    out$sinkhorn_dispatch_reason,
    "dynamic_range_exceeds_scaling_threshold"
  )
  expect_true(length(out$runtime_provenance$environment) >= 27L)
})

test_that("runtime environment helper restores absent variables", {
  key <- "RFUGW_MATVEC_BLOCKED_MIN_WORK_D"
  old <- Sys.getenv(key, unset = NA_character_)
  if (!is.na(old)) Sys.unsetenv(key)
  on.exit({
    if (!is.na(old)) Sys.setenv(RFUGW_MATVEC_BLOCKED_MIN_WORK_D = old)
  }, add = TRUE)

  observed <- with_restored_env(
    stats::setNames("12345", key),
    rfugw:::.solver_runtime_provenance()$environment[[key]]
  )
  expect_identical(observed$source, "environment")
  expect_identical(observed$parsed, 12345)
  expect_identical(Sys.getenv(key, unset = NA_character_), NA_character_)
})

test_that("different effective native behavior is explained by diagnostics", {
  off <- with_restored_env(
    c(RFUGW_ENABLE_WARM_START = "0"),
    rfugw:::cpp_runtime_controls()
  )
  on <- with_restored_env(
    c(RFUGW_ENABLE_WARM_START = "1"),
    rfugw:::cpp_runtime_controls()
  )
  expect_false(off$enable_warm_start)
  expect_true(on$enable_warm_start)

  defaulted <- fgw_entropic(
    matrix(c(0, 1, 1, 0), 2, 2),
    matrix(c(0, 1, 1, 0), 2, 2),
    matrix(c(0, 1, 1, 0), 2, 2),
    max_iter = 2L
  )
  expect_identical(
    defaulted$runtime_provenance$automatic$precision$requested_source,
    "default"
  )
  expect_identical(
    defaulted$runtime_provenance$automatic$structure_rank$requested_source,
    "default"
  )
})
