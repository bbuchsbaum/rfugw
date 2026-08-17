test_that("flagship solvers return rfugw_result with stable accessors", {
  d <- list(
    C1 = matrix(c(0, 1, 1, 0), 2, 2),
    C2 = matrix(c(0, 1, 1, 0), 2, 2),
    M = matrix(c(0, 0.5, 0.5, 0), 2, 2)
  )
  out <- fgw_entropic(d$M, d$C1, d$C2, epsilon = 0.1, max_iter = 20L)
  expect_s3_class(out, "rfugw_result")
  expect_true(is.matrix(rfugw_plan(out)))
  expect_true(is.finite(rfugw_value(out)))
  expect_true(rfugw_status(out) %in% c("converged", "max_iter"))
  res <- rfugw_residuals(out)
  expect_true(is.finite(res$residual))
  printed <- paste(capture.output(print(out)), collapse = "\n")
  expect_false(grepl("0.0000000 0.0000000 0.0000000", printed))
  expect_true(grepl("formulation", printed))
  sm <- summary(out)
  expect_s3_class(sm, "summary.rfugw_result")
  expect_equal(sm$value, rfugw_value(out))
})

test_that("legacy fields remain available on the result object", {
  d <- list(
    C1 = matrix(c(0, 1, 1, 0), 2, 2),
    C2 = matrix(c(0, 1, 1, 0), 2, 2),
    M = matrix(0, 2, 2)
  )
  out <- fgw_exact_cg(d$M, d$C1, d$C2, max_iter = 10L)
  expect_true(is.matrix(out$plan))
  expect_true(is.finite(out$fgw_dist))
  expect_equal(rfugw_value(out), out$fgw_dist)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(out, tmp)
  restored <- readRDS(tmp)
  expect_s3_class(restored, "rfugw_result")
  expect_equal(rfugw_value(restored), rfugw_value(out))
})
