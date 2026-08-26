# Independent prototype used by the KL structural-loss decision record.
# This file is evidence tooling, not a package API.

kl_structural_loss <- function(a, b, log_floor = 0) {
  if (length(a) != length(b)) {
    stop("`a` and `b` must have the same length.", call. = FALSE)
  }
  if (any(!is.finite(a)) || any(!is.finite(b)) || any(a < 0) || any(b < 0)) {
    stop("KL structural inputs must be finite and nonnegative.", call. = FALSE)
  }
  if (length(log_floor) != 1L || !is.finite(log_floor) || log_floor < 0) {
    stop("`log_floor` must be one finite nonnegative scalar.", call. = FALSE)
  }

  if (log_floor > 0) {
    return(a * log(a + log_floor) - a + b - a * log(b + log_floor))
  }

  out <- numeric(length(a))
  source_zero <- a == 0
  target_zero <- b == 0
  out[source_zero] <- b[source_zero]
  out[!source_zero & target_zero] <- Inf
  regular <- !source_zero & !target_zero
  out[regular] <-
    a[regular] * log(a[regular] / b[regular]) - a[regular] + b[regular]
  out
}

gw_kl_enumerated <- function(C1, C2, T, log_floor = 0) {
  stopifnot(
    is.matrix(C1), is.matrix(C2), is.matrix(T),
    nrow(C1) == ncol(C1), nrow(C2) == ncol(C2),
    identical(dim(T), c(nrow(C1), nrow(C2)))
  )
  value <- 0
  for (i in seq_len(nrow(C1))) {
    for (j in seq_len(nrow(C2))) {
      for (k in seq_len(nrow(C1))) {
        for (l in seq_len(nrow(C2))) {
          loss <- kl_structural_loss(C1[i, k], C2[j, l], log_floor)
          value <- value + loss * T[i, j] * T[k, l]
        }
      }
    }
  }
  value
}

gw_kl_factorized <- function(C1, C2, T, log_floor = 0) {
  stopifnot(
    is.matrix(C1), is.matrix(C2), is.matrix(T),
    nrow(C1) == ncol(C1), nrow(C2) == ncol(C2),
    identical(dim(T), c(nrow(C1), nrow(C2)))
  )
  if (log_floor == 0 && any(C2 == 0)) {
    return(gw_kl_enumerated(C1, C2, T, log_floor = 0))
  }

  fC1 <- C1 * log(C1 + log_floor) - C1
  fC2 <- C2
  hC1 <- C1
  hC2 <- log(C2 + log_floor)
  p <- rowSums(T)
  q <- colSums(T)
  constC <- outer(as.vector(fC1 %*% p), rep(1, nrow(C2))) +
    outer(rep(1, nrow(C1)), as.vector(fC2 %*% q))
  tensor <- constC - hC1 %*% T %*% t(hC2)
  sum(tensor * T)
}
