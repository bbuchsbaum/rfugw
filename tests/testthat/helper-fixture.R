fixture_path <- function(name) {
  installed <- system.file("extdata", "fixtures", name, package = "rfugw")
  if (nzchar(installed) && file.exists(installed)) {
    return(installed)
  }
  src <- file.path(testthat::test_path("..", "..", "inst", "extdata", "fixtures"), name)
  if (file.exists(src)) {
    return(src)
  }
  stop("Fixture not found: ", name, call. = FALSE)
}

read_fixture <- function(name) {
  jsonlite::fromJSON(fixture_path(name), simplifyVector = TRUE)
}
