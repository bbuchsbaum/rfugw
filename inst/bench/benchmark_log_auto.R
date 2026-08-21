#!/usr/bin/env Rscript

if (file.exists("DESCRIPTION") && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(rfugw)
}

time_case <- function(label, fun, reps = 3L) {
  runs <- vector("list", reps)
  elapsed <- numeric(reps)
  for (i in seq_len(reps)) {
    timing <- system.time(runs[[i]] <- fun())
    elapsed[[i]] <- unname(timing[["elapsed"]])
  }
  out <- runs[[which.min(abs(elapsed - stats::median(elapsed)))]]
  data.frame(
    case = label,
    median_elapsed_s = stats::median(elapsed),
    result_bytes = as.numeric(object.size(out)),
    status = out$status,
    effective_method = out$effective_sinkhorn_method,
    dynamic_range = out$sinkhorn_dynamic_range,
    threshold = out$sinkhorn_scaling_threshold,
    stringsAsFactors = FALSE
  )
}

set.seed(20260820)
n_bal <- 64L
M_bal <- matrix(runif(n_bal * n_bal), n_bal, n_bal)
p_bal <- stats::rexp(n_bal); p_bal <- p_bal / sum(p_bal)
q_bal <- stats::rexp(n_bal); q_bal <- q_bal / sum(q_bal)

moderate_scaling <- ot_sinkhorn_unbalanced(
  M_bal, p_bal, q_bal, epsilon = 0.1, rho = 2,
  method = "scaling", max_iter = 2000L, tol = 1e-6
)
moderate_log <- ot_sinkhorn_unbalanced(
  M_bal, p_bal, q_bal, epsilon = 0.1, rho = 2,
  method = "log", max_iter = 2000L, tol = 1e-6
)
if (max(abs(moderate_scaling$plan - moderate_log$plan)) > 5e-7 ||
    abs(moderate_scaling$ot_dist - moderate_log$ot_dist) > 5e-6) {
  stop("Correctness gate failed: scaling/log UOT plans disagree.")
}

rows <- list(
  time_case("uot_n64_scaling", function() ot_sinkhorn_unbalanced(
    M_bal, p_bal, q_bal, epsilon = 0.1, rho = 2,
    method = "scaling", max_iter = 2000L, tol = 1e-6
  )),
  time_case("uot_n64_log", function() ot_sinkhorn_unbalanced(
    M_bal, p_bal, q_bal, epsilon = 0.1, rho = 2,
    method = "log", max_iter = 2000L, tol = 1e-6
  )),
  time_case("uot_n64_auto_safe", function() ot_sinkhorn_unbalanced(
    M_bal, p_bal, q_bal, epsilon = 0.1, rho = 2,
    method = "auto", max_iter = 2000L, tol = 1e-6
  )),
  time_case("uot_n64_auto_adversarial", function() ot_sinkhorn_unbalanced(
    1000 * M_bal, p_bal, q_bal, epsilon = 0.01, rho = 2,
    method = "auto", max_iter = 2000L, tol = 1e-6
  ))
)

result <- do.call(rbind, rows)
output <- if (length(commandArgs(trailingOnly = TRUE))) {
  commandArgs(trailingOnly = TRUE)[[1]]
} else {
  "inst/bench/log-auto-baseline.csv"
}
utils::write.csv(result, output, row.names = FALSE)
print(result, row.names = FALSE)
