library(rfugw)
library(bench)

set.seed(42)
ns <- 80
nt <- 90
X1 <- matrix(rnorm(ns * 3), ns, 3)
X2 <- matrix(rnorm(nt * 3), nt, 3)
Y1 <- matrix(rnorm(ns * 5), ns, 5)
Y2 <- matrix(rnorm(nt * 5), nt, 5)

C1 <- as.matrix(dist(X1)); C1 <- C1 / max(C1)
C2 <- as.matrix(dist(X2)); C2 <- C2 / max(C2)
M <- as.matrix(dist(rbind(Y1, Y2)))[1:ns, (ns + 1):(ns + nt)]
M <- M / max(M)

res <- bench::mark(
  fgw = fgw_entropic(
    M, C1, C2,
    alpha = 0.5,
    epsilon = 0.05,
    max_iter = 200,
    sinkhorn_max_iter = 500,
    sinkhorn_tol = 1e-9
  ),
  iterations = 5,
  check = FALSE
)

print(res)
