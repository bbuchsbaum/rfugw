# Entropy-Regularized Fixed-Mass Partial Optimal Transport

Solves `sum(M * G) + epsilon * sum(G * (log(G) - 1))` over nonnegative
subcouplings with `rowSums(G) <= p`, `colSums(G) <= q`, and
`sum(G) = mass`. Zero plan entries contribute zero to the entropy. This
counting-measure, entropy-minus-one convention is explicit and differs
from product-reference KL conventions.

## Usage

``` r
ot_partial_sinkhorn(
  M,
  p = NULL,
  q = NULL,
  mass = 1,
  epsilon = 0.1,
  method = c("auto", "scaling", "log"),
  max_iter = 10000L,
  tol = 1e-09,
  check_every = 10L,
  init_state = NULL,
  verbose = FALSE
)
```

## Arguments

- M:

  Finite source-by-target transport cost.

- p, q:

  Nonnegative source and target weights. Each positive-mass vector is
  normalized to a probability measure, consistently with
  [`ot_partial_emd()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_emd.md).

- mass:

  Requested transported probability mass in `[0, 1]`.

- epsilon:

  Positive entropy coefficient.

- method:

  `"auto"`, `"scaling"`, or genuine `"log"` Dykstra.

- max_iter:

  Maximum Dykstra cycles.

- tol:

  Requested update and feasibility tolerance.

- check_every:

  Iteration interval for convergence diagnostics.

- init_state:

  Optional certified `warm_state` or prior converged result from the
  identical problem and effective backend. Raw plan starts are not
  accepted because they do not determine valid Dykstra corrections.

- verbose:

  If `TRUE`, print checked residuals.

## Value

An `rfugw_result` with the plan, exact entropy convention and objective
decomposition, fixed-mass feasibility, KKT/duality certificate,
dynamic-range dispatch, reusable warm state, and diagnostic trace.

## Details

A genuine log-domain Dykstra backend is available for large dynamic
range. The scaling backend is accepted only inside its declared range;
`"auto"` dispatches safely. Returned convergence additionally requires
primal feasibility, independent objective recomputation, dual
feasibility, complementarity, stationarity, and a primal-dual gap
certificate.

## Examples

``` r
M <- matrix(c(0, 2, 1, 0.3), 2, 2)
out <- ot_partial_sinkhorn(M, mass = 0.7, epsilon = 0.2)
out$status
#> [1] "converged"
sum(rfugw_plan(out))
#> [1] 0.7
```
