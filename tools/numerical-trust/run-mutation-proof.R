#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1L]] else file.path(".gate", "numerical-trust")
source(file.path("inst", "numerical-trust", "mutation-proof-lib.R"))
report <- run_rfugw_mutation_proof(output_dir)
print(report[c("schema_version", "seed", "all_killed", "outputs")])
if (!isTRUE(report$all_killed)) quit(status = 1L)
