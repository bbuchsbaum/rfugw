#!/usr/bin/env Rscript
# Fail if Rcpp attributes drifted from the tracked generated files.

if (!file.exists("DESCRIPTION")) {
  stop("Run tools/check-generated.R from the package root.", call. = FALSE)
}

tmp <- tempfile("rfugw-gen-")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
file.copy("R/RcppExports.R", file.path(tmp, "RcppExports.R"))
file.copy("src/RcppExports.cpp", file.path(tmp, "RcppExports.cpp"))

Rcpp::compileAttributes()

ok <- TRUE
for (rel in c("R/RcppExports.R", "src/RcppExports.cpp")) {
  old <- readLines(file.path(tmp, basename(rel)))
  new <- readLines(rel)
  if (!identical(old, new)) {
    message("Generated file drifted: ", rel)
    ok <- FALSE
  }
}

ns <- readLines("NAMESPACE")
export_names <- sub("^export\\(([^)]+)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
man_names <- sub("\\.Rd$", "", list.files("man", pattern = "\\.Rd$"))
missing_rd <- setdiff(export_names, man_names)
if (length(missing_rd)) {
  message("Exports without Rd topics: ", paste(missing_rd, collapse = ", "))
  ok <- FALSE
}

if (!ok) {
  quit(status = 1L)
}
cat("Generated Rcpp files and Rd coverage are consistent.\n")
