bench_test_resource <- function(name) {
  installed <- system.file("bench", name, package = "rfugw")
  if (nzchar(installed) && file.exists(installed)) {
    return(installed)
  }

  source_path <- testthat::test_path("..", "..", "inst", "bench", name)
  if (file.exists(source_path)) {
    return(source_path)
  }

  stop("Benchmark resource not found: ", name, call. = FALSE)
}
