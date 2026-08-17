#!/usr/bin/env Rscript
# Check that pkgdown reference coverage matches public exports.

if (!file.exists("DESCRIPTION") || !file.exists("_pkgdown.yml")) {
  stop("Run tools/check-docs.R from the package root with _pkgdown.yml present.", call. = FALSE)
}

ns <- readLines("NAMESPACE")
exports <- sub("^export\\(([^)]+)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))

yaml <- readLines("_pkgdown.yml")
listed <- unique(trimws(sub("^\\s*-\\s*", "", grep("^\\s*-\\s+[A-Za-z._0-9]+\\s*$", yaml, value = TRUE))))
listed <- listed[!grepl("^(title|desc|contents):", listed)]

# Allow section selectors such as starts_with("ot_")
selectors <- grep("starts_with|ends_with|matches|has_keyword", yaml, value = TRUE)
covered <- listed
if (any(grepl("starts_with\\(\"ot_\"\\)", selectors))) {
  covered <- union(covered, grep("^ot_", exports, value = TRUE))
}
if (any(grepl("starts_with\\(\"rfugw_\"\\)", selectors))) {
  covered <- union(covered, grep("^rfugw_", exports, value = TRUE))
}
if (any(grepl("starts_with\\(\"multialign_\"\\)", selectors))) {
  covered <- union(covered, grep("^multialign_", exports, value = TRUE))
}

missing <- setdiff(exports, covered)
if (length(missing)) {
  stop(
    "Public exports missing from _pkgdown.yml reference:\n  ",
    paste(missing, collapse = "\n  "),
    call. = FALSE
  )
}
cat("pkgdown reference covers all exported functions.\n")
