make_thread_set <- function(n, seed, d_struct = 3L, d_feat = 2L) {
  set.seed(seed)
  X <- matrix(rnorm(n * d_struct), n, d_struct)
  F <- matrix(rnorm(n * d_feat), n, d_feat)
  C <- as.matrix(dist(X))
  C <- C / max(C)
  list(C = C, F = F, w = rep(1 / n, n))
}

test_that("flagship single-subject solvers ignore OMP_NUM_THREADS", {
  d <- make_thread_set(8L, 20260816L)
  tpl <- make_thread_set(7L, 20260817L)
  M <- as.matrix(dist(rbind(d$F, tpl$F)))[seq_len(8L), 8L + seq_len(7L)]
  old <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
  on.exit({
    if (is.na(old)) Sys.unsetenv("OMP_NUM_THREADS") else Sys.setenv(OMP_NUM_THREADS = old)
  })

  Sys.setenv(OMP_NUM_THREADS = "1")
  a <- fgw_entropic(
    M, d$C, tpl$C, d$w, tpl$w,
    alpha = 0.5, epsilon = 0.08, max_iter = 40L, tol = 1e-8
  )
  Sys.setenv(OMP_NUM_THREADS = "4")
  b <- fgw_entropic(
    M, d$C, tpl$C, d$w, tpl$w,
    alpha = 0.5, epsilon = 0.08, max_iter = 40L, tol = 1e-8
  )
  expect_equal(a$fgw_dist, b$fgw_dist, tolerance = 1e-12)
  expect_equal(a$plan, b$plan, tolerance = 1e-12)
})

test_that("batched multialign 1-thread and 2-thread plans match", {
  s1 <- c(make_thread_set(12L, 301L), list(id = "s1"))
  s2 <- c(make_thread_set(11L, 302L), list(id = "s2"))
  s3 <- c(make_thread_set(10L, 303L), list(id = "s3"))
  tpl <- c(make_thread_set(9L, 304L), list(id = "tpl"))
  ctrl <- list(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    max_iter = 60L,
    sinkhorn_max_iter = 200L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE
  )
  one <- do.call(multialign_fit, c(ctrl, list(n_threads = 1L)))
  two <- do.call(multialign_fit, c(ctrl, list(n_threads = 2L)))
  expect_true(isTRUE(one$used_cpp_batch))
  expect_true(isTRUE(two$used_cpp_batch))
  expect_equal(one$used_threads, 1L)
  expect_gte(two$used_threads, 1L)
  expect_equal(one$objective_total, two$objective_total, tolerance = 1e-10)
  expect_equal(one$couplings$s1, two$couplings$s1, tolerance = 1e-10)
  expect_equal(one$couplings$s2, two$couplings$s2, tolerance = 1e-10)
  expect_equal(one$couplings$s3, two$couplings$s3, tolerance = 1e-10)
})

test_that("sampled GW random_state is reproducible and deterministic sampling matches", {
  set.seed(11L)
  n <- 10L
  X1 <- matrix(rnorm(n * 3L), n, 3L)
  X2 <- matrix(rnorm(n * 3L), n, 3L)
  C1 <- as.matrix(dist(X1)); C1 <- C1 / max(C1)
  C2 <- as.matrix(dist(X2)); C2 <- C2 / max(C2)
  p <- rep(1 / n, n)
  q <- p

  a <- sampled_gromov_wasserstein(
    C1, C2, p, q,
    nb_samples_grad = c(4L, 2L),
    epsilon = 0.1,
    max_iter = 20L,
    random_state = 99L
  )
  b <- sampled_gromov_wasserstein(
    C1, C2, p, q,
    nb_samples_grad = c(4L, 2L),
    epsilon = 0.1,
    max_iter = 20L,
    random_state = 99L
  )
  expect_equal(a, b, tolerance = 1e-12)

  det_a <- sampled_gromov_wasserstein_coords(
    X1, X2, p, q,
    nb_samples_grad = c(5L, 2L),
    epsilon = 0.08,
    max_iter = 25L,
    sampling = "deterministic",
    use_cpp = TRUE
  )
  det_b <- sampled_gromov_wasserstein_coords(
    X1, X2, p, q,
    nb_samples_grad = c(5L, 2L),
    epsilon = 0.08,
    max_iter = 25L,
    sampling = "deterministic",
    use_cpp = TRUE
  )
  expect_equal(det_a, det_b, tolerance = 1e-12)
})

test_that("structure_knn keeps dense storage", {
  C <- matrix(runif(64L), 8L, 8L)
  C <- (C + t(C)) / 2
  diag(C) <- 0
  sparse_metric <- rfugw:::.knn_sparsify_structure(C, k = 2L)
  expect_equal(dim(sparse_metric), dim(C))
  expect_false(inherits(sparse_metric, "sparseMatrix"))
  expect_gte(as.numeric(object.size(sparse_metric)), as.numeric(object.size(C)) * 0.8)
})
