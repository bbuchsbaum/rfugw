#!/usr/bin/env Rscript
# Check that pkgdown reference coverage matches public exports.

if (!file.exists("DESCRIPTION") || !file.exists("_pkgdown.yml")) {
  stop("Run tools/check-docs.R from the package root with _pkgdown.yml present.", call. = FALSE)
}

ns <- readLines("NAMESPACE")
exports <- sub("^export\\(([^)]+)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
s3_lines <- grep("^S3method\\(", ns, value = TRUE)
s3_methods <- gsub(
  "\\s+", "",
  sub("^S3method\\(([^,]+),([^)]+)\\)$", "\\1.\\2", s3_lines)
)
public_topics <- union(exports, s3_methods)

yaml <- readLines("_pkgdown.yml")
listed <- unique(trimws(sub("^\\s*-\\s*", "", grep("^\\s*-\\s+[A-Za-z._0-9]+\\s*$", yaml, value = TRUE))))
listed <- listed[!grepl("^(title|desc|contents):", listed)]

# Allow section selectors such as starts_with("ot_")
selectors <- grep("starts_with|ends_with|matches|has_keyword", yaml, value = TRUE)
covered <- listed
if (any(grepl("starts_with\\(\"ot_\"\\)", selectors))) {
  covered <- union(covered, grep("^ot_", public_topics, value = TRUE))
}
if (any(grepl("starts_with\\(\"rfugw_\"\\)", selectors))) {
  covered <- union(covered, grep("^rfugw_", public_topics, value = TRUE))
}
if (any(grepl("starts_with\\(\"multialign_\"\\)", selectors))) {
  covered <- union(covered, grep("^multialign_", public_topics, value = TRUE))
}

missing <- setdiff(public_topics, covered)
if (length(missing)) {
  stop(
    "Public exports missing from _pkgdown.yml reference:\n  ",
    paste(missing, collapse = "\n  "),
    call. = FALSE
  )
}
cat("pkgdown reference covers all exported functions and registered S3 methods.\n")
