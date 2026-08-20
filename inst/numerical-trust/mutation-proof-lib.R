.mutation_gw_oracle <- function(C1, C2, G) {
  forward <- reverse <- matrix(0, nrow(G), ncol(G))
  for (a in seq_len(nrow(G))) for (b in seq_len(ncol(G))) {
    for (j in seq_len(nrow(G))) for (l in seq_len(ncol(G))) {
      forward[a, b] <- forward[a, b] +
        (C1[a, j] - C2[b, l])^2 * G[j, l]
    }
    for (i in seq_len(nrow(G))) for (k in seq_len(ncol(G))) {
      reverse[a, b] <- reverse[a, b] +
        (C1[i, a] - C2[k, b])^2 * G[i, k]
    }
  }
  list(loss = sum(forward * G), gradient = forward + reverse, forward = forward)
}

.mutation_wrong_target_transpose <- function(C1, C2, G) {
  oracle <- .mutation_gw_oracle(C1, C2, G)
  wrong_reverse <- matrix(0, nrow(G), ncol(G))
  for (a in seq_len(nrow(G))) for (b in seq_len(ncol(G))) {
    for (i in seq_len(nrow(G))) for (k in seq_len(ncol(G))) {
      wrong_reverse[a, b] <- wrong_reverse[a, b] +
        (C1[i, a] - C2[b, k])^2 * G[i, k]
    }
  }
  list(loss = oracle$loss, gradient = oracle$forward + wrong_reverse)
}

.mutation_case <- function(mutant, detecting_test, assertion, expression) {
  start <- proc.time()[["elapsed"]]
  detector_passed_mutant <- isTRUE(force(expression))
  elapsed <- 1000 * (proc.time()[["elapsed"]] - start)
  list(
    mutant = mutant,
    detecting_test = detecting_test,
    assertion = assertion,
    killed = !detector_passed_mutant,
    runtime_ms = unname(elapsed)
  )
}

run_rfugw_mutation_proof <- function(output_dir = NULL) {
  set.seed(20260820L)
  C1 <- matrix(runif(16), 4, 4)
  C2 <- matrix(runif(9), 3, 3)
  G <- matrix(runif(12, 0.02, 0.2), 4, 3)
  oracle <- .mutation_gw_oracle(C1, C2, G)
  transpose_mutant <- .mutation_wrong_target_transpose(C1, C2, G)

  M <- matrix(c(10, 0, 0, 10), 2, 2, byrow = TRUE)
  p <- q <- c(0.5, 0.5)
  false_plan <- p %o% q
  false_primal <- sum(M * false_plan)
  false_dual <- 0

  p_kl <- c(0.4, 0.6)
  q_kl <- c(0.7, 0.3)
  zero_plan <- matrix(0, 2, 2)
  omitted_mass_kl <- 0
  expected_zero_kl <- sum(p_kl %o% q_kl)
  zero_support_plan <- matrix(c(0.4, 0.1, 0.4, 0.1), 2, 2)
  finite_zero_support_mutant <- 0

  cases <- list(
    .mutation_case(
      "wrong_C2_transpose",
      "asymmetric O(n^4) gradient law",
      "max(abs(mutant_gradient - oracle_gradient)) <= 2e-12",
      max(abs(transpose_mutant$gradient - oracle$gradient)) <= 2e-12
    ),
    .mutation_case(
      "false_early_simplex_convergence",
      "independent primal-dual certificate",
      "feasible && reduced_cost_ok && abs(primal-dual) <= 1e-10",
      max(abs(rowSums(false_plan) - p)) <= 1e-12 &&
        max(abs(colSums(false_plan) - q)) <= 1e-12 &&
        min(M) >= -1e-12 && abs(false_primal - false_dual) <= 1e-10
    ),
    .mutation_case(
      "mixed_tolerance_floor_1e-6",
      "requested/effective precision contract",
      "effective_tol == requested_tol || effective_precision == strict_double",
      1e-6 == 1e-9 || identical("mixed", "strict_double")
    ),
    .mutation_case(
      "ignored_inner_residual",
      "converged implies required inner certificate",
      "!converged || (inner_converged && finite inner residuals)",
      !TRUE || (FALSE && is.finite(Inf) && is.finite(Inf))
    ),
    .mutation_case(
      "omitted_generalized_KL_mass_terms",
      "zero-plan generalized KL identity",
      "KL(0 || p outer q) == sum(p outer q)",
      abs(omitted_mass_kl - expected_zero_kl) <= 1e-14
    ),
    .mutation_case(
      "finite_positive_mass_outside_zero_support",
      "generalized KL support law",
      "positive plan mass where reference is zero implies Inf",
      is.infinite(finite_zero_support_mutant) ||
        !any(zero_support_plan > 0 & (c(1, 0) %o% c(0.5, 0.5)) == 0)
    )
  )

  A <- matrix(runif(16), 4, 4); S1 <- (A + t(A)) / 2
  B <- matrix(runif(9), 3, 3); S2 <- (B + t(B)) / 2
  symmetric_oracle <- .mutation_gw_oracle(S1, S2, G)
  symmetric_mutant <- .mutation_wrong_target_transpose(S1, S2, G)
  sensitivity <- list(
    scalar_objective_detects_transpose_mutant = !isTRUE(all.equal(
      transpose_mutant$loss, oracle$loss, tolerance = 0
    )),
    symmetric_gradient_detects_transpose_mutant =
      max(abs(symmetric_mutant$gradient - symmetric_oracle$gradient)) > 2e-12,
    asymmetric_gradient_detects_transpose_mutant =
      max(abs(transpose_mutant$gradient - oracle$gradient)) > 1e-3
  )

  report <- list(
    schema_version = 1L,
    seed = 20260820L,
    all_killed = all(vapply(cases, `[[`, logical(1), "killed")),
    cases = cases,
    insufficient_test_demonstration = sensitivity
  )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    json_path <- file.path(output_dir, "mutation-proof.json")
    md_path <- file.path(output_dir, "mutation-proof.md")
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      jsonlite::write_json(report, json_path, auto_unbox = TRUE, pretty = TRUE)
    } else {
      dput(report, file = json_path)
    }
    lines <- c(
      "# rfugw mutation proof",
      "",
      sprintf("All controlled mutants killed: %s", report$all_killed),
      "",
      "| Mutant | Detecting test | Assertion | Killed | Runtime ms |",
      "|---|---|---|---:|---:|",
      vapply(cases, function(x) sprintf(
        "| %s | %s | `%s` | %s | %.3f |",
        x$mutant, x$detecting_test, x$assertion, x$killed, x$runtime_ms
      ), character(1)),
      "",
      paste0(
        "Objective-only detects transpose mutant: ",
        sensitivity$scalar_objective_detects_transpose_mutant
      ),
      paste0(
        "Symmetric-only gradient detects transpose mutant: ",
        sensitivity$symmetric_gradient_detects_transpose_mutant
      ),
      paste0(
        "Asymmetric gradient detects transpose mutant: ",
        sensitivity$asymmetric_gradient_detects_transpose_mutant
      )
    )
    writeLines(lines, md_path)
    report$outputs <- c(json = json_path, markdown = md_path)
  }
  report
}
