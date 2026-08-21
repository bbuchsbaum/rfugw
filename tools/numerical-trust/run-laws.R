#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_arg <- function(prefix, default) {
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}
family <- value_arg("--family=", "pr")
scope <- value_arg("--scope=", Sys.getenv("RFUGW_TRUST_SCOPE", "pr"))
installed <- "--installed" %in% args
Sys.setenv(RFUGW_TRUST_SCOPE = scope)

filters <- c(
  gw = "gw-gradient-laws",
  transport = "transport-simplex-certificates",
  convergence = "convergence-certificates|nested-solver-status",
  precision = "precision-tolerance-contract",
  validation = "validation-boundaries|ot-utils",
  paths = "advertised-path-laws|solver-invariants|threading|sampled-budget",
  mutation = "mutation-proof",
  pr = paste(
    c(
      "audit-regressions", "gw-gradient-laws", "transport-simplex-certificates",
      "convergence-certificates", "nested-solver-status",
      "precision-tolerance-contract", "validation-boundaries",
      "advertised-path-laws", "mutation-proof", "trust-pyramid"
    ),
    collapse = "|"
  ),
  all = ""
)
if (!family %in% names(filters)) {
  stop("Unknown numerical-trust family: ", family, call. = FALSE)
}

cat(sprintf("rfugw numerical laws: family=%s scope=%s installed=%s\n",
            family, scope, installed))
cat(sprintf("replay: RFUGW_TRUST_SCOPE=%s Rscript tools/numerical-trust/run-laws.R --family=%s\n",
            scope, family))

if (installed) {
  library(rfugw)
  testthat::test_dir(
    "tests/testthat",
    filter = filters[[family]],
    reporter = "summary",
    stop_on_failure = TRUE,
    package = "rfugw",
    load_package = "installed"
  )
} else {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Source-law execution requires package `devtools`.", call. = FALSE)
  }
  devtools::test(
    filter = filters[[family]],
    reporter = "summary",
    stop_on_failure = TRUE
  )
}
