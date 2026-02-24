fixture_path <- function(name) {
  file.path(testthat::test_path("..", "..", "inst", "extdata", "fixtures"), name)
}

read_fixture <- function(name) {
  jsonlite::fromJSON(fixture_path(name), simplifyVector = TRUE)
}
