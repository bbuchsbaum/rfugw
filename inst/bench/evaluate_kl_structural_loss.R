#!/usr/bin/env Rscript
# Reproduce the evidence behind inst/kl-structural-loss-decision.md.

args <- commandArgs(trailingOnly = TRUE)
out_csv <- if (length(args)) args[[1L]] else
  "inst/bench/kl-structural-loss-baseline.csv"
if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the rfugw package root.", call. = FALSE)
}
source("inst/numerical-trust/kl-structural-loss-prototype.R")

p <- c(0.2, 0.3, 0.5)
q <- c(0.4, 0.35, 0.25)
T <- outer(p, q)
positive <- list(
  C1 = matrix(c(0.7, 0.3, 0.5, 0.3, 0.9, 0.4, 0.5, 0.4, 0.8), 3L),
  C2 = matrix(c(0.6, 0.2, 0.45, 0.2, 0.75, 0.35, 0.45, 0.35, 0.95), 3L)
)
near_zero <- list(
  C1 = positive$C1,
  C2 = positive$C2
)
near_zero$C2[1, 2] <- near_zero$C2[2, 1] <- 1e-16
points1 <- c(0, 1, 3)
points2 <- c(0, 2, 5)
zero_diagonal <- list(
  C1 = abs(outer(points1, points1, "-")),
  C2 = abs(outer(points2, points2, "-"))
)

cases <- list(
  positive = positive,
  near_zero = near_zero,
  zero_diagonal_distance = zero_diagonal
)
floors <- c(0, 1e-18, 1e-15, 1e-12)
elapsed_ms <- function(expr, repetitions = 25L) {
  expr <- substitute(expr)
  1000 * system.time({
    for (i in seq_len(repetitions)) eval(expr, envir = parent.frame())
  })[["elapsed"]] / repetitions
}

rows <- list()
index <- 1L
for (case_name in names(cases)) {
  case <- cases[[case_name]]
  for (floor in floors) {
    enumerated <- gw_kl_enumerated(case$C1, case$C2, T, floor)
    factorized <- gw_kl_factorized(case$C1, case$C2, T, floor)
    rows[[index]] <- data.frame(
      regime = case_name,
      log_floor = floor,
      enumerated_objective = enumerated,
      factorized_objective = factorized,
      absolute_difference = if (is.finite(enumerated) && is.finite(factorized))
        abs(enumerated - factorized) else NA_real_,
      enumerated_ms = elapsed_ms(
        gw_kl_enumerated(case$C1, case$C2, T, floor)
      ),
      factorized_ms = elapsed_ms(
        gw_kl_factorized(case$C1, case$C2, T, floor)
      ),
      r_version = paste(R.version$major, R.version$minor, sep = "."),
      platform = R.version$platform,
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }
}
result <- do.call(rbind, rows)
utils::write.csv(result, out_csv, row.names = FALSE, na = "")
print(result)
cat("wrote ", out_csv, "\n", sep = "")
