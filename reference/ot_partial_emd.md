# Exact partial linear optimal transport

Solves the nonnegative-cost partial transport problem `min <M, G>`
subject to `rowSums(G) <= p`, `colSums(G) <= q`, and `sum(G) = mass`.
Weights are normalized to probability vectors, so `mass` is in `[0, 1]`.
The implementation reduces the problem to the certified exact transport
primitive with one dummy source and target; it does not add a second
simplex implementation.

## Usage

``` r
ot_partial_emd(M, p = NULL, q = NULL, mass = 1, max_iter = 20000L, tol = 1e-12)
```

## Arguments

- M:

  Cost matrix (`ns x nt`).

- p:

  Source weights (default uniform). Renormalized to sum 1.

- q:

  Target weights (default uniform). Renormalized to sum 1.

- mass:

  Transported probability mass in `[0, 1]`.

- max_iter:

  Maximum simplex iterations.

- tol:

  Optimality tolerance.

## Value

An `rfugw_result` with `plan`, `ot_dist`, `partial_ot_dist`, partial
feasibility and mass certificates, and the exact certificate for the
equivalent augmented transport problem. Augmented dual potentials
include the final dummy coordinate.

## Examples

``` r
M <- matrix(c(0, 3, 2, 0), 2, 2)
out <- ot_partial_emd(M, mass = 0.5)
out$status
#> [1] "converged"
sum(rfugw_plan(out))
#> [1] 0.5
```
