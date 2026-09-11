moment_fugw_oracle_path <- system.file(
  "bench", "moment_fugw_oracle.R", package = "rfugw"
)
if (!nzchar(moment_fugw_oracle_path) ||
    !file.exists(moment_fugw_oracle_path)) {
  moment_fugw_oracle_path <- testthat::test_path(
    "..", "..", "inst", "bench", "moment_fugw_oracle.R"
  )
}
if (!file.exists(moment_fugw_oracle_path)) {
  stop("Cannot locate the independent Moment-FUGW oracle.", call. = FALSE)
}
source(moment_fugw_oracle_path, local = TRUE)
