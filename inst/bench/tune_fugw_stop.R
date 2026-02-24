suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) {
    .libPaths(c(rlib, .libPaths()))
  }
  library(rfugw)
  library(jsonlite)
  library(bench)
})

parse_numeric_list <- function(text, fallback) {
  if (!nzchar(text)) {
    return(fallback)
  }
  vals <- suppressWarnings(as.numeric(strsplit(text, ",", fixed = TRUE)[[1L]]))
  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) == 0L) {
    return(fallback)
  }
  unique(vals)
}

with_env <- function(env, expr) {
  keys <- names(env)
  old <- Sys.getenv(keys, unset = NA_character_)
  on.exit({
    restore <- old[!is.na(old)]
    cleared <- names(old)[is.na(old)]
    if (length(restore) > 0L) {
      do.call(Sys.setenv, as.list(restore))
    }
    if (length(cleared) > 0L) {
      Sys.unsetenv(cleared)
    }
  }, add = TRUE)
  do.call(Sys.setenv, as.list(env))
  force(expr)
}

resolve_fixture_dir <- function(fixture_dir) {
  candidates <- unique(c(fixture_dir, "inst/extdata/fixtures", "rfugw/inst/extdata/fixtures"))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop("Could not locate fixture directory")
  }
  hit[[1L]]
}

read_fixture <- function(name, fixture_dir) {
  jsonlite::fromJSON(file.path(fixture_dir, name), simplifyVector = TRUE)
}

make_fugw_problem <- function(n) {
  X1 <- matrix(rnorm(n * 3), n, 3)
  X2 <- matrix(rnorm(n * 3), n, 3)
  Cx <- as.matrix(dist(X1))
  Cy <- as.matrix(dist(X2))
  Cx <- Cx / max(Cx)
  Cy <- Cy / max(Cy)

  M <- matrix(1, nrow = n, ncol = n)
  idx <- seq_len(n)
  M[cbind(idx, rev(idx))] <- 0
  list(Cx = Cx, Cy = Cy, M = M)
}

bench_one <- function(fn, iters) {
  b <- bench::mark(
    run = fn(),
    iterations = iters,
    check = FALSE,
    memory = FALSE,
    filter_gc = TRUE
  )
  as.numeric(tibble::as_tibble(b)[1, "median"][[1]]) * 1000
}

dominance_frontier <- function(speed, err) {
  n <- length(speed)
  dom <- rep(FALSE, n)
  for (i in seq_len(n)) {
    if (dom[i]) next
    for (j in seq_len(n)) {
      if (i == j) next
      if (speed[j] <= speed[i] && err[j] <= err[i] &&
          (speed[j] < speed[i] || err[j] < err[i])) {
        dom[i] <- TRUE
        break
      }
    }
  }
  !dom
}

args <- commandArgs(trailingOnly = TRUE)
iters <- if (length(args) >= 1L) as.integer(args[[1L]]) else 2L
out_csv <- if (length(args) >= 2L) args[[2L]] else "inst/bench/results/fugw_stop_sweep.csv"
out_frontier_csv <- if (length(args) >= 3L) args[[3L]] else "inst/bench/results/fugw_stop_sweep_frontier.csv"
seed <- if (length(args) >= 4L) as.integer(args[[4L]]) else 42L
validate_top_n <- if (length(args) >= 5L) as.integer(args[[5L]]) else 8L

set.seed(seed)
fixture_dir <- resolve_fixture_dir(Sys.getenv("RFUGW_FIXTURE_DIR", unset = "inst/extdata/fixtures"))
fu <- read_fixture("fugw_kl_sinkhorn_fixture.json", fixture_dir)

sizes <- c(30L, 50L, 70L, 90L)
problems <- lapply(sizes, make_fugw_problem)

call_fugw <- function(problem) {
  rfugw::fugw_kl(
    Cx = problem$Cx,
    Cy = problem$Cy,
    reg_marginals = c(100, 50),
    epsilon = 1e-2,
    alpha = 0.5,
    M = problem$M,
    precision = "double",
    max_iter = 80,
    max_iter_ot = 300,
    tol = 1e-8,
    tol_ot = 1e-8
  )
}

call_fixture <- function() {
  rfugw::fugw_kl(
    Cx = fu$inputs$Cx,
    Cy = fu$inputs$Cy,
    wx = fu$inputs$wx,
    wy = fu$inputs$wy,
    reg_marginals = unlist(fu$params$reg_marginals),
    epsilon = fu$params$epsilon,
    alpha = fu$params$alpha,
    M = fu$inputs$M,
    max_iter = fu$params$max_iter,
    tol = fu$params$tol,
    max_iter_ot = fu$params$max_iter_ot,
    tol_ot = fu$params$tol_ot,
    precision = "double"
  )
}

target_mult <- parse_numeric_list(
  Sys.getenv("RFUGW_SWEEP_TARGET_MULT", unset = ""),
  c(20, 28, 36, 40, 44, 52)
)
target_cap_d <- parse_numeric_list(
  Sys.getenv("RFUGW_SWEEP_TARGET_CAP_D", unset = ""),
  c(7e-5, 1e-4, 1.5e-4)
)
rel_floor_d <- parse_numeric_list(
  Sys.getenv("RFUGW_SWEEP_REL_FLOOR_D", unset = ""),
  c(1e-5, 2e-5, 4e-5)
)
col_rel_floor_d <- parse_numeric_list(
  Sys.getenv("RFUGW_SWEEP_COL_REL_FLOOR_D", unset = ""),
  c(7e-5, 1e-4, 1.5e-4)
)

settings <- expand.grid(
  target_mult = target_mult,
  target_cap_d = target_cap_d,
  rel_floor_d = rel_floor_d,
  col_rel_floor_d = col_rel_floor_d,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Ensure the current default is explicitly present.
default_row <- data.frame(
  target_mult = 52,
  target_cap_d = 1.5e-4,
  rel_floor_d = 4e-5,
  col_rel_floor_d = 1.5e-4,
  stringsAsFactors = FALSE
)
settings <- unique(rbind(default_row, settings))
settings <- settings[order(settings$target_mult, settings$target_cap_d, settings$rel_floor_d, settings$col_rel_floor_d), ]
settings$setting_id <- seq_len(nrow(settings))

run_one_setting <- function(row) {
  env <- c(
    RFUGW_SINKHORN_TARGET_TOL_MULT = format(row$target_mult, scientific = TRUE, digits = 16),
    RFUGW_SINKHORN_TARGET_TOL_CAP_D = format(row$target_cap_d, scientific = TRUE, digits = 16),
    RFUGW_SINKHORN_REL_TOL_FLOOR_D = format(row$rel_floor_d, scientific = TRUE, digits = 16),
    RFUGW_SINKHORN_COL_REL_TOL_FLOOR_D = format(row$col_rel_floor_d, scientific = TRUE, digits = 16)
  )

  with_env(env, {
    fx <- call_fixture()
    obj_abs <- abs(fx$fugw_cost - fu$outputs$fugw_cost)
    pi_samp_abs <- max(abs(fx$pi_samp - fu$outputs$pi_samp))
    pi_feat_abs <- max(abs(fx$pi_feat - fu$outputs$pi_feat))

    tol_obj <- 5e-4
    tol_pi_samp <- 2e-4
    tol_pi_feat <- 3e-4
    fixture_ratio <- max(
      obj_abs / tol_obj,
      pi_samp_abs / tol_pi_samp,
      pi_feat_abs / tol_pi_feat
    )

    med <- numeric(length(sizes))
    inner <- numeric(length(sizes))
    outer <- numeric(length(sizes))
    for (i in seq_along(sizes)) {
      p <- problems[[i]]
      warm <- call_fugw(p)
      inner[[i]] <- warm$inner_iters_total
      outer[[i]] <- warm$iterations
      med[[i]] <- bench_one(function() call_fugw(p), iters = iters)
    }

    data.frame(
      setting_id = row$setting_id,
      target_mult = row$target_mult,
      target_cap_d = row$target_cap_d,
      rel_floor_d = row$rel_floor_d,
      col_rel_floor_d = row$col_rel_floor_d,
      fixture_obj_abs = obj_abs,
      fixture_pi_samp_abs = pi_samp_abs,
      fixture_pi_feat_abs = pi_feat_abs,
      fixture_ratio = fixture_ratio,
      fixture_ok = fixture_ratio <= 1,
      median_ms_n30 = med[[1L]],
      median_ms_n50 = med[[2L]],
      median_ms_n70 = med[[3L]],
      median_ms_n90 = med[[4L]],
      median_ms_mean = mean(med),
      median_ms_geom = exp(mean(log(pmax(med, 1e-12)))),
      inner_total_n30 = inner[[1L]],
      inner_total_n50 = inner[[2L]],
      inner_total_n70 = inner[[3L]],
      inner_total_n90 = inner[[4L]],
      outer_iters_n30 = outer[[1L]],
      outer_iters_n50 = outer[[2L]],
      outer_iters_n70 = outer[[3L]],
      outer_iters_n90 = outer[[4L]],
      stringsAsFactors = FALSE
    )
  })
}

cat(sprintf("Sweeping %d settings (iters=%d)...\n", nrow(settings), iters))
rows <- vector("list", nrow(settings))
for (i in seq_len(nrow(settings))) {
  rows[[i]] <- run_one_setting(settings[i, , drop = FALSE])
  if (i %% 10L == 0L || i == nrow(settings)) {
    cat(sprintf("  completed %d/%d\n", i, nrow(settings)))
  }
}
res <- do.call(rbind, rows)

res$is_pareto <- dominance_frontier(res$median_ms_geom, res$fixture_ratio)
res$is_pareto_fixture_ok <- FALSE
ok_idx <- which(res$fixture_ok)
if (length(ok_idx) > 0L) {
  po <- dominance_frontier(res$median_ms_geom[ok_idx], res$fixture_ratio[ok_idx])
  res$is_pareto_fixture_ok[ok_idx] <- po
}

# Optional full-gate validation for top frontier settings.
res$gate_ok <- NA
res$gate_fail_count <- NA_integer_
res$gate_max_ratio <- NA_real_

front_ok <- res[res$is_pareto_fixture_ok, , drop = FALSE]
front_ok <- front_ok[order(front_ok$median_ms_geom), , drop = FALSE]

if (nrow(front_ok) > 0L && validate_top_n > 0L) {
  gate_script <- if (file.exists("inst/bench/accuracy_gate.R")) {
    "inst/bench/accuracy_gate.R"
  } else {
    "rfugw/inst/bench/accuracy_gate.R"
  }
  source(gate_script)

  top_n <- min(validate_top_n, nrow(front_ok))
  cat(sprintf("Running full accuracy gate on top %d frontier settings...\n", top_n))
  for (i in seq_len(top_n)) {
    sid <- front_ok$setting_id[[i]]
    row <- res[res$setting_id == sid, , drop = FALSE]
    env <- c(
      RFUGW_SINKHORN_TARGET_TOL_MULT = format(row$target_mult, scientific = TRUE, digits = 16),
      RFUGW_SINKHORN_TARGET_TOL_CAP_D = format(row$target_cap_d, scientific = TRUE, digits = 16),
      RFUGW_SINKHORN_REL_TOL_FLOOR_D = format(row$rel_floor_d, scientific = TRUE, digits = 16),
      RFUGW_SINKHORN_COL_REL_TOL_FLOOR_D = format(row$col_rel_floor_d, scientific = TRUE, digits = 16)
    )
    tbl <- with_env(env, run_accuracy_gate(
      fixture_dir = fixture_dir,
      stop_on_fail = FALSE,
      verbose = FALSE
    ))
    fail_n <- sum(!tbl$ok)
    max_ratio <- max(tbl$max_abs / pmax(tbl$tol, .Machine$double.eps))
    res$gate_ok[res$setting_id == sid] <- fail_n == 0L
    res$gate_fail_count[res$setting_id == sid] <- fail_n
    res$gate_max_ratio[res$setting_id == sid] <- max_ratio
    cat(sprintf("  setting_id=%d gate_ok=%s max_ratio=%.3f\n", sid, ifelse(fail_n == 0L, "TRUE", "FALSE"), max_ratio))
  }
}

res <- res[order(res$median_ms_geom, res$fixture_ratio), ]
frontier <- res[res$is_pareto, ]
frontier <- frontier[order(frontier$median_ms_geom, frontier$fixture_ratio), ]

if (nrow(frontier) > 0L) {
  # Prefer safe + validated settings first in printed summary.
  ord_safe <- order(
    -(frontier$fixture_ok),
    -(ifelse(is.na(frontier$gate_ok), FALSE, frontier$gate_ok)),
    frontier$median_ms_geom,
    frontier$fixture_ratio
  )
  frontier <- frontier[ord_safe, ]
}

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(res, out_csv, row.names = FALSE)
write.csv(frontier, out_frontier_csv, row.names = FALSE)

cat("\nTop frontier settings:\n")
print(utils::head(frontier[, c(
  "setting_id",
  "target_mult",
  "target_cap_d",
  "rel_floor_d",
  "col_rel_floor_d",
  "median_ms_geom",
  "fixture_ratio",
  "fixture_ok",
  "gate_ok"
)], 12L), row.names = FALSE)

cat(sprintf("\nWrote sweep CSV: %s\n", out_csv))
cat(sprintf("Wrote frontier CSV: %s\n", out_frontier_csv))
