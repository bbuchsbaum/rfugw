suppressPackageStartupMessages({
  library(utils)
})

args <- commandArgs(trailingOnly = TRUE)
in_csv <- if (length(args) >= 1) args[[1]] else "inst/bench/results/benchmark_sparse_sampled_latest.csv"
out_md <- if (length(args) >= 2) args[[2]] else "inst/bench/results/benchmark_sparse_sampled_report_latest.md"

tbl <- read.csv(in_csv, stringsAsFactors = FALSE)

fmt <- function(x) ifelse(is.finite(x), sprintf("%.3f", x), "NA")

lookup_ms <- function(sub, method) {
  hit <- sub$median_ms[sub$method == method]
  if (length(hit) == 0L) return(NA_real_)
  as.numeric(hit[[1]])
}

lines <- c(
  "# Sparse Sampled GW Benchmark Report",
  "",
  sprintf("- input CSV: `%s`", in_csv),
  "",
  "## Raw Results",
  "",
  "| suite | n | method | median_ms | mem_bytes | objective_gap_dense |",
  "|---|---:|---|---:|---:|---:|"
)

for (i in seq_len(nrow(tbl))) {
  lines <- c(lines, sprintf(
    "| %s | %d | %s | %s | %.0f | %s |",
    tbl$suite[[i]],
    tbl$n[[i]],
    tbl$method[[i]],
    fmt(tbl$median_ms[[i]]),
    tbl$mem_bytes[[i]],
    fmt(tbl$objective_gap_dense[[i]])
  ))
}

lines <- c(lines, "", "## Speed Ratios", "")
lines <- c(lines, "| suite | n | ratio | value |", "|---|---:|---|---:|")

for (s in unique(tbl$suite)) {
  sub_s <- tbl[tbl$suite == s, , drop = FALSE]
  for (nn in unique(sub_s$n)) {
    sub <- sub_s[sub_s$n == nn, , drop = FALSE]
    med <- setNames(sub$median_ms, sub$method)

    if ("rfugw_sampled_dense_solver" %in% names(med) &&
        "rfugw_sampled_coords_solver" %in% names(med)) {
      lines <- c(lines, sprintf(
        "| %s | %d | dense_solver_over_coords_solver | %s |",
        s, nn, fmt(med[["rfugw_sampled_dense_solver"]] / med[["rfugw_sampled_coords_solver"]])
      ))
    }

    if ("pot_sampled_dense_solver_from_cost" %in% names(med) &&
        "rfugw_sampled_coords_solver" %in% names(med)) {
      lines <- c(lines, sprintf(
        "| %s | %d | pot_solver_over_rfugw_coords_solver | %s |",
        s, nn, fmt(med[["pot_sampled_dense_solver_from_cost"]] / med[["rfugw_sampled_coords_solver"]])
      ))
    }

    if ("rfugw_sampled_dense_end_to_end" %in% names(med) &&
        "rfugw_sampled_coords_end_to_end" %in% names(med)) {
      lines <- c(lines, sprintf(
        "| %s | %d | dense_e2e_over_coords_e2e | %s |",
        s, nn, fmt(med[["rfugw_sampled_dense_end_to_end"]] / med[["rfugw_sampled_coords_end_to_end"]])
      ))
    }

    if ("pot_sampled_dense_end_to_end_from_coords" %in% names(med) &&
        "rfugw_sampled_coords_end_to_end" %in% names(med)) {
      lines <- c(lines, sprintf(
        "| %s | %d | pot_e2e_over_rfugw_coords_e2e | %s |",
        s, nn, fmt(med[["pot_sampled_dense_end_to_end_from_coords"]] / med[["rfugw_sampled_coords_end_to_end"]])
      ))
    }

    if ("rfugw_sparse_graph_pipeline" %in% names(med) &&
        "rfugw_sampled_coords_end_to_end" %in% names(med)) {
      lines <- c(lines, sprintf(
        "| %s | %d | graph_pipeline_over_coords_e2e | %s |",
        s, nn, fmt(med[["rfugw_sparse_graph_pipeline"]] / med[["rfugw_sampled_coords_end_to_end"]])
      ))
    }
  }
}

lines <- c(lines, "", "## Crossover Map (POT vs rfugw coords)", "")
lines <- c(lines, "| suite | n | rfugw_coords_ms | pot_ms | pot_over_rfugw | winner |", "|---|---:|---:|---:|---:|---|")

for (nn in sort(unique(tbl$n))) {
  sub_solver <- tbl[tbl$suite == "sampled_sparse_solver_only" & tbl$n == nn, , drop = FALSE]
  r_solver <- lookup_ms(sub_solver, "rfugw_sampled_coords_solver")
  p_solver <- lookup_ms(sub_solver, "pot_sampled_dense_solver_from_cost")
  ratio_solver <- if (is.finite(r_solver) && is.finite(p_solver) && r_solver > 0) p_solver / r_solver else NA_real_
  winner_solver <- if (is.finite(ratio_solver)) {
    if (ratio_solver >= 1) "rfugw" else "pot"
  } else {
    "NA"
  }
  lines <- c(lines, sprintf(
    "| sampled_sparse_solver_only | %d | %s | %s | %s | %s |",
    nn, fmt(r_solver), fmt(p_solver), fmt(ratio_solver), winner_solver
  ))

  sub_e2e <- tbl[tbl$suite == "sampled_sparse_end_to_end" & tbl$n == nn, , drop = FALSE]
  r_e2e <- lookup_ms(sub_e2e, "rfugw_sampled_coords_end_to_end")
  p_e2e <- lookup_ms(sub_e2e, "pot_sampled_dense_end_to_end_from_coords")
  ratio_e2e <- if (is.finite(r_e2e) && is.finite(p_e2e) && r_e2e > 0) p_e2e / r_e2e else NA_real_
  winner_e2e <- if (is.finite(ratio_e2e)) {
    if (ratio_e2e >= 1) "rfugw" else "pot"
  } else {
    "NA"
  }
  lines <- c(lines, sprintf(
    "| sampled_sparse_end_to_end | %d | %s | %s | %s | %s |",
    nn, fmt(r_e2e), fmt(p_e2e), fmt(ratio_e2e), winner_e2e
  ))
}

writeLines(lines, out_md)
cat(sprintf("Wrote sparse sampled benchmark report: %s\n", out_md))
