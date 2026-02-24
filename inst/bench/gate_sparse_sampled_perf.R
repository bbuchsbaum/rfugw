suppressPackageStartupMessages({
  library(utils)
})

args <- commandArgs(trailingOnly = TRUE)
in_csv <- if (length(args) >= 1) args[[1]] else "inst/bench/results/benchmark_sparse_sampled_latest.csv"
min_n <- if (length(args) >= 2) as.integer(args[[2]]) else 400L
min_ratio_solver <- if (length(args) >= 3) as.numeric(args[[3]]) else 1.00
min_ratio_e2e <- if (length(args) >= 4) as.numeric(args[[4]]) else 1.00
max_obj_gap <- if (length(args) >= 5) as.numeric(args[[5]]) else 5e-3

if (!file.exists(in_csv)) {
  stop("Benchmark CSV not found: ", in_csv, call. = FALSE)
}

tbl <- read.csv(in_csv, stringsAsFactors = FALSE)
if (!all(c("suite", "n", "method", "median_ms", "objective_gap_dense") %in% names(tbl))) {
  stop("Benchmark CSV missing required columns.", call. = FALSE)
}

fmt <- function(x) ifelse(is.finite(x), sprintf("%.4f", x), "NA")
lookup_ms <- function(sub, method) {
  hit <- sub$median_ms[sub$method == method]
  if (length(hit) == 0L) return(NA_real_)
  as.numeric(hit[[1]])
}

targets <- sort(unique(tbl$n[tbl$n >= min_n]))
if (length(targets) == 0L) {
  stop(sprintf("No benchmark rows found with n >= %d", min_n), call. = FALSE)
}

fail_msgs <- character(0)
out_rows <- list()
k <- 1L

for (nn in targets) {
  sub_solver <- tbl[tbl$suite == "sampled_sparse_solver_only" & tbl$n == nn, , drop = FALSE]
  r_solver <- lookup_ms(sub_solver, "rfugw_sampled_coords_solver")
  p_solver <- lookup_ms(sub_solver, "pot_sampled_dense_solver_from_cost")
  ratio_solver <- if (is.finite(r_solver) && is.finite(p_solver) && r_solver > 0) p_solver / r_solver else NA_real_
  if (!is.finite(ratio_solver)) {
    fail_msgs <- c(fail_msgs, sprintf("missing solver-only comparison at n=%d", nn))
  } else if (ratio_solver < min_ratio_solver) {
    fail_msgs <- c(
      fail_msgs,
      sprintf("solver-only ratio below floor at n=%d: %.4f < %.4f", nn, ratio_solver, min_ratio_solver)
    )
  }
  out_rows[[k]] <- data.frame(
    suite = "sampled_sparse_solver_only",
    n = nn,
    rfugw_ms = r_solver,
    pot_ms = p_solver,
    pot_over_rfugw = ratio_solver,
    floor = min_ratio_solver,
    stringsAsFactors = FALSE
  )
  k <- k + 1L

  sub_e2e <- tbl[tbl$suite == "sampled_sparse_end_to_end" & tbl$n == nn, , drop = FALSE]
  r_e2e <- lookup_ms(sub_e2e, "rfugw_sampled_coords_end_to_end")
  p_e2e <- lookup_ms(sub_e2e, "pot_sampled_dense_end_to_end_from_coords")
  ratio_e2e <- if (is.finite(r_e2e) && is.finite(p_e2e) && r_e2e > 0) p_e2e / r_e2e else NA_real_
  if (!is.finite(ratio_e2e)) {
    fail_msgs <- c(fail_msgs, sprintf("missing end-to-end comparison at n=%d", nn))
  } else if (ratio_e2e < min_ratio_e2e) {
    fail_msgs <- c(
      fail_msgs,
      sprintf("end-to-end ratio below floor at n=%d: %.4f < %.4f", nn, ratio_e2e, min_ratio_e2e)
    )
  }
  out_rows[[k]] <- data.frame(
    suite = "sampled_sparse_end_to_end",
    n = nn,
    rfugw_ms = r_e2e,
    pot_ms = p_e2e,
    pot_over_rfugw = ratio_e2e,
    floor = min_ratio_e2e,
    stringsAsFactors = FALSE
  )
  k <- k + 1L
}

max_gap <- suppressWarnings(max(tbl$objective_gap_dense, na.rm = TRUE))
if (is.finite(max_gap) && max_gap > max_obj_gap) {
  fail_msgs <- c(
    fail_msgs,
    sprintf("objective_gap_dense exceeded max: %.6f > %.6f", max_gap, max_obj_gap)
  )
}

out <- do.call(rbind, out_rows)
cat("== Sparse Sampled Performance Gate ==\n")
print(within(out, {
  rfugw_ms <- round(rfugw_ms, 3)
  pot_ms <- round(pot_ms, 3)
  pot_over_rfugw <- round(pot_over_rfugw, 4)
  floor <- round(floor, 4)
}))
cat(sprintf("max objective_gap_dense: %s (limit %s)\n", fmt(max_gap), fmt(max_obj_gap)))

if (length(fail_msgs) > 0L) {
  stop(paste(c("Performance gate failed:", fail_msgs), collapse = "\n- "), call. = FALSE)
}

cat("Performance gate passed.\n")
