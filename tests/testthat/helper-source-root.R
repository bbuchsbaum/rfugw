rfugw_test_source_root <- function() {
  candidates <- c(
    test_path("../.."),
    test_path("../../00_pkg_src/rfugw")
  )
  is_source_root <- vapply(candidates, function(path) {
    file.exists(file.path(path, "DESCRIPTION")) &&
      dir.exists(file.path(path, "src"))
  }, logical(1))
  if (!any(is_source_root)) {
    stop(
      "Cannot locate the rfugw source tree required by native contract tests.",
      call. = FALSE
    )
  }
  normalizePath(candidates[which(is_source_root)[[1L]]], mustWork = TRUE)
}
