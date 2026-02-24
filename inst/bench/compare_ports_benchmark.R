suppressPackageStartupMessages({
  library(utils)
})

args <- commandArgs(trailingOnly = TRUE)
rfugw_csv <- if (length(args) >= 1) args[[1]] else "inst/bench/results/benchmark_ports_latest.csv"
pot_csv <- if (length(args) >= 2) args[[2]] else "inst/bench/results/benchmark_ports_pot_latest.csv"
out_md <- if (length(args) >= 3) args[[3]] else "inst/bench/results/benchmark_ports_report_latest.md"

rf <- read.csv(rfugw_csv, stringsAsFactors = FALSE)
pot <- read.csv(pot_csv, stringsAsFactors = FALSE)

map <- data.frame(
  rfugw = c(
    "rfugw_entropic_partial_gw",
    "rfugw_entropic_partial_fgw",
    "rfugw_semirelaxed_gw",
    "rfugw_semirelaxed_fgw",
    "rfugw_sampled_gw",
    "rfugw_ucoot_kl"
  ),
  pot = c(
    "pot_entropic_partial_gw",
    "pot_entropic_partial_fgw",
    "pot_semirelaxed_gw",
    "pot_semirelaxed_fgw",
    "pot_sampled_gw",
    "pot_ucoot_kl"
  ),
  stringsAsFactors = FALSE
)

rows <- list()
k <- 1L
for (i in seq_len(nrow(map))) {
  sub_rf <- rf[rf$method == map$rfugw[[i]], c("suite", "n", "median_ms")]
  sub_p <- pot[pot$method == map$pot[[i]], c("suite", "n", "median_ms")]
  if (nrow(sub_rf) == 0L || nrow(sub_p) == 0L) next

  merged <- merge(sub_rf, sub_p, by = c("suite", "n"), suffixes = c("_rfugw", "_pot"))
  if (nrow(merged) == 0L) next

  merged$method_rfugw <- map$rfugw[[i]]
  merged$method_pot <- map$pot[[i]]
  merged$speedup_rfugw_vs_pot <- merged$median_ms_pot / merged$median_ms_rfugw

  rows[[k]] <- merged
  k <- k + 1L
}

out <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame()

fmt_num <- function(x) ifelse(is.finite(x), sprintf("%.3f", x), "NA")

lines <- c(
  "# Benchmark Ports Report",
  "",
  sprintf("- rfugw CSV: `%s`", rfugw_csv),
  sprintf("- POT CSV: `%s`", pot_csv),
  "",
  "## Matched Methods",
  "",
  "| suite | n | rfugw_method | pot_method | rfugw_median_ms | pot_median_ms | speedup_pot_over_rfugw |",
  "|---|---:|---|---|---:|---:|---:|"
)

if (nrow(out) > 0L) {
  for (i in seq_len(nrow(out))) {
    lines <- c(lines, sprintf(
      "| %s | %d | %s | %s | %s | %s | %s |",
      out$suite[[i]],
      out$n[[i]],
      out$method_rfugw[[i]],
      out$method_pot[[i]],
      fmt_num(out$median_ms_rfugw[[i]]),
      fmt_num(out$median_ms_pot[[i]]),
      fmt_num(out$speedup_rfugw_vs_pot[[i]])
    ))
  }

  lines <- c(lines, "", "## Geometric Mean Speed Ratio by Suite", "")
  suites <- unique(out$suite)
  lines <- c(lines, "| suite | geom_mean(pot_ms / rfugw_ms) |", "|---|---:|")
  for (s in suites) {
    vals <- out$speedup_rfugw_vs_pot[out$suite == s]
    vals <- vals[is.finite(vals) & vals > 0]
    gm <- if (length(vals) == 0L) NA_real_ else exp(mean(log(vals)))
    lines <- c(lines, sprintf("| %s | %s |", s, fmt_num(gm)))
  }
} else {
  lines <- c(lines, "| (none) | | | | | | |")
}

writeLines(lines, out_md)
cat(sprintf("Wrote benchmark ports report: %s\n", out_md))
