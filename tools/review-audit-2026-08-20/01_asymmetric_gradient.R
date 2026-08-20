audit_lib <- Sys.getenv("RFUGW_AUDIT_LIB", unset = "")
if (nzchar(audit_lib)) {
  .libPaths(c(audit_lib, .libPaths()))
}

library(rfugw)

gw_brute <- function(C1, C2, G) {
  out <- 0
  for (i in seq_len(nrow(C1))) {
    for (j in seq_len(ncol(C1))) {
      for (k in seq_len(nrow(C2))) {
        for (l in seq_len(ncol(C2))) {
          out <- out + (C1[i, j] - C2[k, l])^2 * G[i, k] * G[j, l]
        }
      }
    }
  }
  out
}

gw_grad_brute <- function(C1, C2, G) {
  grad <- matrix(0, nrow(G), ncol(G))
  for (a in seq_len(nrow(G))) {
    for (b in seq_len(ncol(G))) {
      forward <- 0
      reverse <- 0
      for (j in seq_len(nrow(G))) {
        for (l in seq_len(ncol(G))) {
          forward <- forward + (C1[a, j] - C2[b, l])^2 * G[j, l]
        }
      }
      for (i in seq_len(nrow(G))) {
        for (k in seq_len(ncol(G))) {
          reverse <- reverse + (C1[i, a] - C2[k, b])^2 * G[i, k]
        }
      }
      grad[a, b] <- forward + reverse
    }
  }
  grad
}

set.seed(20260820)
C1 <- matrix(runif(16, 0.05, 1.4), 4, 4)
C2 <- matrix(runif(9, 0.1, 1.7), 3, 3)
diag(C1) <- 0
diag(C2) <- 0
G <- matrix(runif(12, 0.05, 0.25), 4, 3)
D <- matrix(rnorm(12), 4, 3)
D <- D / sqrt(sum(D^2))

cpp <- rfugw:::cpp_gw_square_terms_square(C1, C2, G, symmetric = FALSE)
brute_loss <- gw_brute(C1, C2, G)
brute_grad <- gw_grad_brute(C1, C2, G)
h <- 1e-6
finite_difference <- (gw_brute(C1, C2, G + h * D) -
  gw_brute(C1, C2, G - h * D)) / (2 * h)
cpp_directional <- sum(cpp$grad * D)
brute_directional <- sum(brute_grad * D)

C1s <- (C1 + t(C1)) / 2
C2s <- (C2 + t(C2)) / 2
cpp_symmetric <- rfugw:::cpp_gw_square_terms_square(C1s, C2s, G, symmetric = TRUE)
cpp_general_on_symmetric <- rfugw:::cpp_gw_square_terms_square(
  C1s, C2s, G, symmetric = FALSE
)

metrics <- c(
  loss_abs_error = abs(cpp$loss - brute_loss),
  gradient_max_abs_error = max(abs(cpp$grad - brute_grad)),
  finite_difference_vs_brute = abs(finite_difference - brute_directional),
  finite_difference_vs_cpp = abs(finite_difference - cpp_directional),
  symmetric_general_loss_error = abs(cpp_symmetric$loss - cpp_general_on_symmetric$loss),
  symmetric_general_gradient_error = max(abs(cpp_symmetric$grad - cpp_general_on_symmetric$grad))
)

print(metrics, digits = 16)

stopifnot(
  metrics[["loss_abs_error"]] < 1e-10,
  metrics[["gradient_max_abs_error"]] > 1e-4,
  metrics[["finite_difference_vs_brute"]] < 1e-7,
  metrics[["finite_difference_vs_cpp"]] > 1e-4,
  metrics[["symmetric_general_loss_error"]] < 1e-10,
  metrics[["symmetric_general_gradient_error"]] < 1e-10
)

sessionInfo()
