make_knn_similarity <- function(X, k = 5L) {
  D <- as.matrix(dist(X))
  sigma <- stats::median(D[D > 0])
  if (!is.finite(sigma) || sigma <= 0) {
    sigma <- 1
  }
  W <- exp(-(D^2) / (2 * sigma^2))
  diag(W) <- 0

  n <- nrow(W)
  out <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    ord <- order(W[i, ], decreasing = TRUE)
    ord <- ord[ord != i]
    keep <- ord[seq_len(min(k, length(ord)))]
    out[i, keep] <- W[i, keep]
  }
  out <- pmax(out, t(out))
  diag(out) <- 0
  out
}

test_that("graph diffusion coordinates are consistent for dense and sparse inputs", {
  skip_if_not_installed("Matrix")
  set.seed(901)

  X <- matrix(rnorm(30L * 3L), 30L, 3L)
  W_dense <- make_knn_similarity(X, k = 6L)
  W_sparse <- Matrix::Matrix(W_dense, sparse = TRUE)

  E_dense <- rfugw::graph_diffusion_coordinates(
    W_dense,
    n_components = 5L,
    diffusion_time = 1,
    self_loop = 1e-6
  )
  E_sparse <- rfugw::graph_diffusion_coordinates(
    W_sparse,
    n_components = 5L,
    diffusion_time = 1,
    self_loop = 1e-6
  )

  expect_equal(dim(E_dense), c(30L, 5L))
  expect_equal(dim(E_sparse), c(30L, 5L))
  expect_true(all(is.finite(E_dense)))
  expect_true(all(is.finite(E_sparse)))
  expect_equal(as.matrix(dist(E_dense)), as.matrix(dist(E_sparse)), tolerance = 1e-5)
})

test_that("graph diffusion coordinates handle disconnected and isolated nodes as expected", {
  W_disconnected <- matrix(0, nrow = 4L, ncol = 4L)
  W_disconnected[1L, 2L] <- 1
  W_disconnected[2L, 1L] <- 1
  W_disconnected[3L, 4L] <- 1
  W_disconnected[4L, 3L] <- 1

  E_ok <- rfugw::graph_diffusion_coordinates(
    W_disconnected,
    n_components = 2L,
    diffusion_time = 1
  )
  expect_equal(dim(E_ok), c(4L, 2L))
  expect_true(all(is.finite(E_ok)))

  W_isolated <- W_disconnected
  W_isolated[4L, ] <- 0
  W_isolated[, 4L] <- 0
  expect_error(
    rfugw::graph_diffusion_coordinates(W_isolated, n_components = 2L),
    "isolated"
  )

  E_fixed <- rfugw::graph_diffusion_coordinates(
    W_isolated,
    n_components = 2L,
    self_loop = 1e-3
  )
  expect_equal(dim(E_fixed), c(4L, 2L))
  expect_true(all(is.finite(E_fixed)))
})

test_that("coord-native sampled GW agrees with dense sampled GW objective on same fixtures", {
  set.seed(902)

  X1 <- matrix(rnorm(16L * 4L), 16L, 4L)
  X2 <- matrix(rnorm(14L * 4L), 14L, 4L)
  C1 <- as.matrix(dist(X1))
  C2 <- as.matrix(dist(X2))

  out_dense <- rfugw::sampled_gromov_wasserstein(
    C1 = C1,
    C2 = C2,
    nb_samples_grad = c(6L, 2L),
    epsilon = 0.05,
    max_iter = 35L,
    random_state = 902L,
    log = TRUE
  )
  out_coords <- rfugw::sampled_gromov_wasserstein_coords(
    X1 = X1,
    X2 = X2,
    nb_samples_grad = c(6L, 2L),
    epsilon = 0.05,
    max_iter = 35L,
    random_state = 902L,
    log = TRUE
  )

  p <- rep(1 / nrow(X1), nrow(X1))
  q <- rep(1 / nrow(X2), nrow(X2))
  val_dense <- rfugw:::.gw_square_value(C1, C2, out_dense$plan, p, q, symmetric = TRUE)
  val_coords <- rfugw:::.gw_square_value(C1, C2, out_coords$plan, p, q, symmetric = TRUE)
  rel_diff <- abs(val_coords - val_dense) / max(1, abs(val_dense))

  expect_lte(rel_diff, 0.08)
  expect_equal(rowSums(out_coords$plan), p, tolerance = 1e-3)
  expect_equal(colSums(out_coords$plan), q, tolerance = 1e-3)
  expect_true(is.finite(out_coords$gw_dist_estimated))
  expect_gte(out_coords$iterations, 1L)
})

test_that("coord-native sampled GW respects source permutation symmetry", {
  set.seed(903)

  ns <- 10L
  nt <- 9L
  X1 <- matrix(rnorm(ns * 3L), ns, 3L)
  X2 <- matrix(rnorm(nt * 3L), nt, 3L)
  perm <- sample.int(ns)

  out_ref <- rfugw::sampled_gromov_wasserstein_coords(
    X1 = X1,
    X2 = X2,
    nb_samples_grad = c(ns, nt),
    epsilon = 0.05,
    max_iter = 20L,
    random_state = 903L,
    log = TRUE
  )
  out_perm <- rfugw::sampled_gromov_wasserstein_coords(
    X1 = X1[perm, , drop = FALSE],
    X2 = X2,
    nb_samples_grad = c(ns, nt),
    epsilon = 0.05,
    max_iter = 20L,
    random_state = 903L,
    log = TRUE
  )

  plan_unperm <- out_perm$plan[order(perm), , drop = FALSE]
  expect_equal(plan_unperm, out_ref$plan, tolerance = 1e-8)

  C1 <- as.matrix(dist(X1))
  C2 <- as.matrix(dist(X2))
  p <- rep(1 / ns, ns)
  q <- rep(1 / nt, nt)
  val_ref <- rfugw:::.gw_square_value(C1, C2, out_ref$plan, p, q, symmetric = TRUE)
  val_perm <- rfugw:::.gw_square_value(C1, C2, plan_unperm, p, q, symmetric = TRUE)
  expect_equal(val_perm, val_ref, tolerance = 1e-10)
})

test_that("coord-native sampled GW C++ path remains aligned with R fallback", {
  set.seed(905)
  old_mixed <- Sys.getenv("RFUGW_SAMPLED_MIXED", unset = NA_character_)
  Sys.setenv(RFUGW_SAMPLED_MIXED = "0")
  on.exit({
    if (is.na(old_mixed)) {
      Sys.unsetenv("RFUGW_SAMPLED_MIXED")
    } else {
      Sys.setenv(RFUGW_SAMPLED_MIXED = old_mixed)
    }
  }, add = TRUE)

  ns <- 20L
  nt <- 19L
  X1 <- matrix(rnorm(ns * 4L), ns, 4L)
  X2 <- matrix(rnorm(nt * 4L), nt, 4L)
  p <- rep(1 / ns, ns)
  q <- rep(1 / nt, nt)
  C1 <- as.matrix(dist(X1))
  C2 <- as.matrix(dist(X2))

  out_cpp <- rfugw::sampled_gromov_wasserstein_coords(
    X1 = X1,
    X2 = X2,
    p = p,
    q = q,
    nb_samples_grad = c(8L, 2L),
    epsilon = 0.05,
    max_iter = 45L,
    random_state = 905L,
    use_cpp = TRUE,
    sampling = "deterministic",
    log = TRUE
  )
  out_r <- rfugw::sampled_gromov_wasserstein_coords(
    X1 = X1,
    X2 = X2,
    p = p,
    q = q,
    nb_samples_grad = c(8L, 2L),
    epsilon = 0.05,
    max_iter = 45L,
    random_state = 905L,
    use_cpp = FALSE,
    sampling = "deterministic",
    log = TRUE
  )

  val_cpp <- rfugw:::.gw_square_value(C1, C2, out_cpp$plan, p, q, symmetric = TRUE)
  val_r <- rfugw:::.gw_square_value(C1, C2, out_r$plan, p, q, symmetric = TRUE)
  rel_diff <- abs(val_cpp - val_r) / max(1, abs(val_r))

  expect_lte(rel_diff, 0.03)
  expect_equal(rowSums(out_cpp$plan), p, tolerance = 8e-4)
  expect_equal(colSums(out_cpp$plan), q, tolerance = 8e-4)
})

test_that("sparse graph to diffusion to sampled GW pipeline returns valid coupling", {
  skip_if_not_installed("Matrix")
  set.seed(904)

  X1 <- matrix(rnorm(18L * 3L), 18L, 3L)
  X2 <- matrix(rnorm(20L * 3L), 20L, 3L)
  W1 <- Matrix::Matrix(make_knn_similarity(X1, k = 5L), sparse = TRUE)
  W2 <- Matrix::Matrix(make_knn_similarity(X2, k = 5L), sparse = TRUE)

  E1 <- rfugw::graph_diffusion_coordinates(
    W = W1,
    n_components = 4L,
    diffusion_time = 1,
    self_loop = 1e-6
  )
  E2 <- rfugw::graph_diffusion_coordinates(
    W = W2,
    n_components = 4L,
    diffusion_time = 1,
    self_loop = 1e-6
  )

  out <- rfugw::sampled_gromov_wasserstein_coords(
    X1 = E1,
    X2 = E2,
    nb_samples_grad = c(7L, 2L),
    epsilon = 0.05,
    max_iter = 60L,
    random_state = 904L,
    log = TRUE
  )

  p <- rep(1 / nrow(E1), nrow(E1))
  q <- rep(1 / nrow(E2), nrow(E2))
  expect_equal(rowSums(out$plan), p, tolerance = 3e-4)
  expect_equal(colSums(out$plan), q, tolerance = 3e-4)
  expect_true(is.finite(out$gw_dist_estimated))
  expect_gte(out$iterations, 1L)
})

test_that("sampled_gw_from_graphs matches explicit two-step pipeline", {
  skip_if_not_installed("Matrix")
  set.seed(906)

  X1 <- matrix(rnorm(15L * 3L), 15L, 3L)
  X2 <- matrix(rnorm(17L * 3L), 17L, 3L)
  W1 <- Matrix::Matrix(make_knn_similarity(X1, k = 4L), sparse = TRUE)
  W2 <- Matrix::Matrix(make_knn_similarity(X2, k = 4L), sparse = TRUE)

  out_wrap <- rfugw::sampled_gw_from_graphs(
    W1 = W1,
    W2 = W2,
    n_components = 5L,
    diffusion_time = 1,
    self_loop = 1e-6,
    nb_samples_grad = c(6L, 2L),
    epsilon = 0.05,
    max_iter = 40L,
    random_state = 906L,
    use_cpp = FALSE,
    log = TRUE,
    return_embeddings = TRUE
  )

  E1 <- rfugw::graph_diffusion_coordinates(W1, n_components = 5L, diffusion_time = 1, self_loop = 1e-6)
  E2 <- rfugw::graph_diffusion_coordinates(W2, n_components = 5L, diffusion_time = 1, self_loop = 1e-6)
  out_manual <- rfugw::sampled_gromov_wasserstein_coords(
    X1 = E1,
    X2 = E2,
    nb_samples_grad = c(6L, 2L),
    epsilon = 0.05,
    max_iter = 40L,
    random_state = 906L,
    use_cpp = FALSE,
    log = TRUE
  )

  expect_equal(out_wrap$plan, out_manual$plan, tolerance = 1e-10)
  expect_equal(out_wrap$gw_dist_estimated, out_manual$gw_dist_estimated, tolerance = 1e-10)
  expect_equal(out_wrap$X1_embed, E1, tolerance = 1e-10)
  expect_equal(out_wrap$X2_embed, E2, tolerance = 1e-10)
})
