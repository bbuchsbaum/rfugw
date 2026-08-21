gw_oracle_tensors <- function(C1, C2, G) {
  forward <- reverse <- matrix(0, nrow(G), ncol(G))
  for (a in seq_len(nrow(G))) {
    for (b in seq_len(ncol(G))) {
      for (j in seq_len(nrow(G))) {
        for (l in seq_len(ncol(G))) {
          forward[a, b] <- forward[a, b] +
            (C1[a, j] - C2[b, l])^2 * G[j, l]
        }
      }
      for (i in seq_len(nrow(G))) {
        for (k in seq_len(ncol(G))) {
          reverse[a, b] <- reverse[a, b] +
            (C1[i, a] - C2[k, b])^2 * G[i, k]
        }
      }
    }
  }
  list(
    forward = forward,
    reverse = reverse,
    loss = sum(forward * G),
    gradient = forward + reverse
  )
}

gw_partial_lp_direction <- function(cost, p, q, mass) {
  skip_if_not_installed("lpSolve")
  ns <- length(p)
  nt <- length(q)
  extended <- matrix(0, ns + 1L, nt + 1L)
  extended[seq_len(ns), seq_len(nt)] <- cost
  extended[ns + 1L, nt + 1L] <- max(cost) + max(abs(cost)) + 1
  fit <- lpSolve::lp.transport(
    cost.mat = extended,
    direction = "min",
    row.signs = rep("=", ns + 1L),
    row.rhs = c(p, sum(q) - mass),
    col.signs = rep("=", nt + 1L),
    col.rhs = c(q, sum(p) - mass),
    integers = integer(0)
  )
  expect_identical(fit$status, 0L)
  fit$solution[seq_len(ns), seq_len(nt), drop = FALSE]
}

gw_line_minimize <- function(objective, G0, Gc) {
  delta <- Gc - G0
  f0 <- objective(G0)
  half <- objective(G0 + 0.5 * delta)
  f1 <- objective(Gc)
  a <- 2 * (f1 - f0) - 4 * (half - f0)
  b <- (f1 - f0) - a
  steps <- c(0, 1)
  if (is.finite(a) && a > 0) {
    stationary <- -b / (2 * a)
    if (stationary > 0 && stationary < 1) steps <- c(steps, stationary)
  }
  values <- vapply(steps, function(step) objective(G0 + step * delta), numeric(1))
  step <- steps[[which.min(values)]]
  list(plan = G0 + step * delta, objective = min(values))
}

gw_entropic_partial_reference <- function(a, b, cost, reg, mass,
                                          max_iter = 2000L, tol = 1e-13) {
  plan <- exp(-cost / reg)
  plan <- plan * (mass / sum(plan))
  q1 <- q2 <- q3 <- matrix(1, nrow(cost), ncol(cost))
  error <- Inf
  iterations <- 0L

  while (iterations < max_iter && error > tol) {
    previous <- plan

    plan <- plan * q1
    projected_rows <- sweep(plan, 1L, pmin(a / pmax(rowSums(plan), 1e-300), 1), "*")
    q1 <- q1 * (previous / pmax(projected_rows, 1e-300))

    previous_rows <- projected_rows
    projected_rows <- projected_rows * q2
    projected_cols <- sweep(
      projected_rows, 2L,
      pmin(b / pmax(colSums(projected_rows), 1e-300), 1),
      "*"
    )
    q2 <- q2 * (previous_rows / pmax(projected_cols, 1e-300))

    previous_cols <- projected_cols
    projected_cols <- projected_cols * q3
    plan <- projected_cols * (mass / pmax(sum(projected_cols), 1e-300))
    q3 <- q3 * (previous_cols / pmax(plan, 1e-300))

    if ((iterations %% 10L) == 0L) {
      error <- sqrt(sum((previous - plan)^2))
    }
    iterations <- iterations + 1L
  }
  plan
}

test_that("canonical general tensors match rectangular O(n^4) oracles", {
  shapes <- list(c(2L, 3L), c(3L, 2L), c(3L, 4L), c(4L, 3L))
  for (seed in seq_len(12L)) {
    set.seed(9000L + seed)
    shape <- shapes[[1L + ((seed - 1L) %% length(shapes))]]
    ns <- shape[[1]]
    nt <- shape[[2]]
    C1 <- matrix(runif(ns * ns, -0.5, 1.5), ns, ns)
    C2 <- matrix(runif(nt * nt, -0.5, 1.5), nt, nt)
    G <- matrix(runif(ns * nt, 0, 0.4), ns, nt)

    oracle <- gw_oracle_tensors(C1, C2, G)
    native <- rfugw:::cpp_gw_square_terms_square(C1, C2, G, symmetric = FALSE)

    expect_equal(native$forward_tensor, oracle$forward, tolerance = 2e-12)
    expect_equal(native$reverse_tensor, oracle$reverse, tolerance = 2e-12)
    expect_equal(native$loss, oracle$loss, tolerance = 2e-12)
    expect_equal(native$grad, oracle$gradient, tolerance = 2e-12)
  }
})

test_that("general gradient obeys directional derivatives and permutations", {
  set.seed(20260820)
  C1 <- matrix(runif(16), 4, 4)
  C2 <- matrix(runif(9), 3, 3)
  G <- matrix(runif(12, 0.02, 0.2), 4, 3)
  direction <- matrix(rnorm(12), 4, 3)
  direction <- direction / sqrt(sum(direction^2))
  h <- 1e-6

  native <- rfugw:::cpp_gw_square_terms_square(C1, C2, G, symmetric = FALSE)
  finite_difference <- (
    gw_oracle_tensors(C1, C2, G + h * direction)$loss -
      gw_oracle_tensors(C1, C2, G - h * direction)$loss
  ) / (2 * h)
  expect_equal(sum(native$grad * direction), finite_difference, tolerance = 2e-7)

  source_perm <- c(4L, 2L, 1L, 3L)
  target_perm <- c(2L, 3L, 1L)
  permuted <- rfugw:::cpp_gw_square_terms_square(
    C1[source_perm, source_perm],
    C2[target_perm, target_perm],
    G[source_perm, target_perm],
    symmetric = FALSE
  )
  expect_equal(permuted$loss, native$loss, tolerance = 2e-12)
  expect_equal(
    permuted$grad,
    native$grad[source_perm, target_perm],
    tolerance = 2e-12
  )
})

test_that("symmetric and general canonical paths agree on symmetric inputs", {
  set.seed(44)
  A <- matrix(runif(16), 4, 4); C1 <- (A + t(A)) / 2
  B <- matrix(runif(9), 3, 3); C2 <- (B + t(B)) / 2
  G <- matrix(runif(12, 0.01, 0.3), 4, 3)

  symmetric <- rfugw:::cpp_gw_square_terms_square(C1, C2, G, symmetric = TRUE)
  general <- rfugw:::cpp_gw_square_terms_square(C1, C2, G, symmetric = FALSE)
  expect_equal(symmetric$loss, general$loss, tolerance = 2e-12)
  expect_equal(symmetric$grad, general$grad, tolerance = 2e-12)
  expect_equal(general$forward_tensor, general$reverse_tensor, tolerance = 2e-12)
})

test_that("one-step asymmetric exact partial FGW matches an independent LP oracle", {
  skip_if_not_installed("lpSolve")
  set.seed(91)
  ns <- 4L
  nt <- 3L
  C1 <- matrix(runif(ns * ns, 0.05, 1.2), ns, ns); diag(C1) <- 0
  C2 <- matrix(runif(nt * nt, 0.05, 1.3), nt, nt); diag(C2) <- 0
  M <- matrix(runif(ns * nt), ns, nt)
  p <- c(0.1, 0.2, 0.3, 0.4)
  q <- c(0.2, 0.3, 0.5)
  mass <- 0.7
  alpha <- 0.6
  G0 <- (p %o% q) * mass
  objective <- function(G) {
    (1 - alpha) * sum(M * G) + alpha * gw_oracle_tensors(C1, C2, G)$loss
  }
  gradient <- (1 - alpha) * M + alpha * gw_oracle_tensors(C1, C2, G0)$gradient
  direction <- gw_partial_lp_direction(gradient, p, q, mass)
  expected <- gw_line_minimize(objective, G0, direction)

  native <- rfugw:::cpp_partial_fgw_exact_square(
    M, C1, C2, p, q, mass, alpha, FALSE, G0,
    1L, 0, 1L, 20000L, 1e-12
  )
  expect_equal(native$plan, expected$plan, tolerance = 2e-9)
  expect_equal(native$objective, expected$objective, tolerance = 2e-10)
})

test_that("one-step asymmetric entropic partial FGW matches a Dykstra oracle", {
  set.seed(37)
  ns <- 3L
  nt <- 4L
  C1 <- matrix(runif(ns * ns, 0.05, 1.1), ns, ns); diag(C1) <- 0
  C2 <- matrix(runif(nt * nt, 0.05, 1.2), nt, nt); diag(C2) <- 0
  M <- matrix(runif(ns * nt), ns, nt)
  p <- c(0.2, 0.3, 0.5)
  q <- c(0.1, 0.2, 0.3, 0.4)
  mass <- 0.65
  alpha <- 0.55
  reg <- 0.2
  G0 <- (p %o% q) * mass
  cost <- (1 - alpha) * M + alpha * gw_oracle_tensors(C1, C2, G0)$gradient
  expected <- gw_entropic_partial_reference(p, q, cost, reg, mass)

  native <- rfugw:::cpp_partial_fgw_entropic_square(
    M, C1, C2, p, q, mass, reg, alpha, FALSE, G0,
    1L, 0, 2000L, 1e-13, 1L
  )
  expect_equal(native$plan, expected, tolerance = 5e-11)
  expect_equal(sum(native$plan), mass, tolerance = 2e-12)
  expect_true(all(rowSums(native$plan) <= p + 2e-10))
  expect_true(all(colSums(native$plan) <= q + 2e-10))
})

.gw_law_replay <- function(seed, ns, nt, label, observed, expected, atol, rtol) {
  sprintf(
    paste0(
      "case=%s seed=%d dims=%dx%d observed=%.17g expected=%.17g ",
      "atol=%.3g rtol=%.3g replay=RFUGW_REPLAY_SEED=%d ",
      "Rscript tools/numerical-trust/run-laws.R --family=gw"
    ),
    label, seed, ns, nt, observed, expected, atol, rtol, seed
  )
}

.expect_gw_directional_law <- function(C1, C2, G, M, alpha, label, seed) {
  set.seed(seed + 700000L)
  direction <- matrix(rnorm(length(G)), nrow(G), ncol(G))
  direction <- direction / sqrt(sum(direction^2))
  # The objective is quadratic in G, so a fixed plan-scale step avoids
  # catastrophic cancellation when structure costs are scaled by 1e5.
  h <- 1e-5 * max(1, max(abs(G)))
  native <- rfugw:::cpp_gw_square_terms_square(C1, C2, G, symmetric = FALSE)
  native_gradient <- (1 - alpha) * M + alpha * native$grad
  objective <- function(plan) {
    (1 - alpha) * sum(M * plan) + alpha * gw_oracle_tensors(C1, C2, plan)$loss
  }
  expected <- (objective(G + h * direction) - objective(G - h * direction)) / (2 * h)
  observed <- sum(native_gradient * direction)
  atol <- 3e-7
  rtol <- 5e-8
  scale <- max(abs(observed), abs(expected), 1)
  expect_true(
    abs(observed - expected) <= atol + rtol * scale,
    info = .gw_law_replay(
      seed, nrow(G), ncol(G), label, observed, expected, atol, rtol
    )
  )
}

test_that("seeded asymmetric GW laws span repeated, zero, and scaled structures", {
  n_cases <- if (identical(Sys.getenv("RFUGW_TRUST_SCOPE"), "nightly")) 120L else 24L
  for (case_id in seq_len(n_cases)) {
    seed <- 20260820L + case_id
    set.seed(seed)
    ns <- 2L + (case_id %% 3L)
    nt <- 2L + ((case_id + 1L) %% 3L)
    C1 <- matrix(runif(ns * ns, -1, 2), ns, ns)
    C2 <- matrix(runif(nt * nt, -1, 2), nt, nt)
    label <- c("dense", "repeated", "zero", "scaled")[[1L + ((case_id - 1L) %% 4L)]]
    if (identical(label, "repeated")) {
      C1[ns, ] <- C1[1L, ]
      C2[, nt] <- C2[, 1L]
    } else if (identical(label, "zero")) {
      C1[] <- 0
      C2[seq.int(1L, length(C2), by = nt + 1L)] <- 0
    } else if (identical(label, "scaled")) {
      C1 <- 1e5 * C1
      C2 <- 1e5 * C2
    }
    G <- matrix(runif(ns * nt, 0.01, 0.4), ns, nt)
    M <- matrix(runif(ns * nt, -0.5, 1.5), ns, nt)
    oracle <- gw_oracle_tensors(C1, C2, G)
    native <- rfugw:::cpp_gw_square_terms_square(C1, C2, G, symmetric = FALSE)
    scale <- max(1, max(abs(oracle$gradient)))
    tolerance <- 5e-12 * scale
    expect_true(
      max(abs(native$grad - oracle$gradient)) <= tolerance,
      info = .gw_law_replay(
        seed, ns, nt, label,
        max(abs(native$grad - oracle$gradient)), 0, tolerance, 0
      )
    )
    .expect_gw_directional_law(C1, C2, G, M, 0.63, label, seed)
  }
})

test_that("FGW alpha endpoints and objective decomposition are independent laws", {
  set.seed(121)
  C1 <- matrix(runif(9, -0.2, 1), 3, 3)
  C2 <- matrix(runif(16, -0.4, 1.3), 4, 4)
  M <- matrix(runif(12, -1, 2), 3, 4)
  G <- matrix(runif(12, 0, 0.3), 3, 4)
  linear <- sum(M * G)
  structural <- gw_oracle_tensors(C1, C2, G)$loss
  for (alpha in c(0, 0.2, 0.5, 1)) {
    observed <- ot_fgw_square(M, C1, C2, G, alpha, symmetric = FALSE)
    expect_equal(
      observed,
      (1 - alpha) * linear + alpha * structural,
      tolerance = 3e-12,
      info = sprintf("alpha=%g asymmetric rectangular decomposition", alpha)
    )
  }
  expect_equal(ot_fgw_square(M, C1, C2, G, 0, FALSE), linear, tolerance = 0)
  expect_equal(ot_fgw_square(M, C1, C2, G, 1, FALSE), structural, tolerance = 3e-12)
})

test_that("standard, partial, and semirelaxed solver plans obey the shared gradient law", {
  skip_if_not_installed("lpSolve")
  seed <- 8231L
  set.seed(seed)
  ns <- 3L
  nt <- 4L
  C1 <- matrix(runif(ns * ns, 0, 1), ns, ns); diag(C1) <- 0
  C2 <- matrix(runif(nt * nt, 0, 1), nt, nt); diag(C2) <- 0
  M <- matrix(runif(ns * nt), ns, nt)
  p <- c(0.2, 0.3, 0.5)
  q <- c(0.1, 0.2, 0.3, 0.4)
  alpha <- 0.57
  paths <- list(
    exact = fgw_exact_cg(
      M, C1, C2, p, q, alpha, symmetric = FALSE,
      max_iter = 2L, lp_solver = "cpp_transport"
    )$plan,
    entropic = fgw_entropic(
      M, C1, C2, p, q, alpha, epsilon = 0.2, symmetric = FALSE,
      max_iter = 2L, precision = "strict_double"
    )$plan,
    partial_exact = partial_fused_gromov_wasserstein(
      M, C1, C2, p, q, m = 0.7, alpha = alpha, symmetric = FALSE,
      numItermax = 2L, log = TRUE
    )$plan,
    partial_entropic = entropic_partial_fused_gromov_wasserstein(
      M, C1, C2, p, q, m = 0.7, reg = 0.2, alpha = alpha,
      symmetric = FALSE, numItermax = 2L, log = TRUE
    )$plan,
    semirelaxed_exact = semirelaxed_fused_gromov_wasserstein(
      M, C1, C2, p, alpha = alpha, symmetric = FALSE,
      max_iter = 2L, log = TRUE
    )$plan,
    semirelaxed_entropic = entropic_semirelaxed_fused_gromov_wasserstein(
      M, C1, C2, p, epsilon = 0.2, alpha = alpha,
      symmetric = FALSE, max_iter = 2L
    )$plan
  )
  for (label in names(paths)) {
    .expect_gw_directional_law(
      C1, C2, paths[[label]], M, alpha, label, seed + match(label, names(paths))
    )
  }
})

test_that("the asymmetric law kills the reviewed target-transpose mutation", {
  set.seed(991)
  C1 <- matrix(runif(16), 4, 4)
  C2 <- matrix(runif(9), 3, 3)
  G <- matrix(runif(12, 0.02, 0.2), 4, 3)
  oracle <- gw_oracle_tensors(C1, C2, G)
  wrong_reverse <- matrix(0, 4, 3)
  for (a in seq_len(4L)) for (b in seq_len(3L)) {
    for (i in seq_len(4L)) for (k in seq_len(3L)) {
      wrong_reverse[a, b] <- wrong_reverse[a, b] +
        (C1[i, a] - C2[b, k])^2 * G[i, k]
    }
  }
  mutant_gradient <- oracle$forward + wrong_reverse
  expect_gt(max(abs(mutant_gradient - oracle$gradient)), 1e-3)
  expect_equal(
    rfugw:::cpp_gw_square_terms_square(C1, C2, G, FALSE)$grad,
    oracle$gradient,
    tolerance = 2e-12
  )
})
