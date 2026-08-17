#!/usr/bin/env Rscript
# Representative solver smoke used by sanitizer and serial-build CI.

library(rfugw)
set.seed(16)
ns <- 6
nt <- 7
C1 <- as.matrix(dist(matrix(rnorm(ns * 2), ns, 2))); C1 <- C1 / max(C1)
C2 <- as.matrix(dist(matrix(rnorm(nt * 2), nt, 2))); C2 <- C2 / max(C2)
M <- as.matrix(dist(rbind(matrix(rnorm(ns * 2), ns, 2), matrix(rnorm(nt * 2), nt, 2))))[1:ns, ns + 1:nt]
M <- M / max(M)

sink <- ot_sinkhorn(M, epsilon = 0.08, max_iter = 200L)
stopifnot(is.finite(sink$ot_dist), sink$status %in% c("converged", "max_iter"))
emd <- ot_emd(M)
stopifnot(is.finite(emd$ot_dist), isTRUE(emd$converged) || identical(emd$status, "max_iter"))
unb <- ot_sinkhorn_unbalanced(M, epsilon = 0.08, rho = 4, max_iter = 80L)
stopifnot(is.finite(unb$ot_dist))
ent <- fgw_entropic(M, C1, C2, epsilon = 0.08, max_iter = 40L)
stopifnot(is.finite(ent$fgw_dist), is.finite(ent$residual))
exact <- fgw_exact_cg(M, C1, C2, max_iter = 20L)
stopifnot(is.finite(exact$fgw_dist))
fugw <- fugw_kl(C1, C2, M = M, epsilon = 0.02, max_iter = 15L)
stopifnot(is.finite(fugw$fugw_cost))
cat("Flagship smoke passed.\n")
