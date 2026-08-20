trust_test_resource <- function(...) {
  parts <- c(...)
  installed <- do.call(system.file, c(as.list(parts), list(package = "rfugw")))
  if (nzchar(installed) && file.exists(installed)) {
    return(installed)
  }
  source_path <- do.call(
    testthat::test_path,
    c(list("..", "..", "inst"), as.list(parts))
  )
  if (file.exists(source_path)) return(source_path)
  stop("Numerical trust resource not found: ", file.path(parts), call. = FALSE)
}
