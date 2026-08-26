# Construct an explicit entropic fixed-mass partial transport problem

Construct an explicit entropic fixed-mass partial transport problem

## Usage

``` r
transport_problem_partial_sinkhorn(
  M,
  p = NULL,
  q = NULL,
  mass = 1,
  epsilon = 0.1,
  method = c("auto", "scaling", "log"),
  max_iter = 10000L,
  tol = 1e-09,
  check_every = 10L,
  mass_policy = NULL
)
```

## Arguments

- M:

  Finite source-by-target cost matrix.

- p, q:

  Nonnegative source and target weights.

- mass:

  Transported probability mass in `[0, 1]`.

- epsilon:

  Positive coefficient of `KL(plan || p %o% q)`, in the same units as
  `M`.

- method:

  Requested Sinkhorn method: `"scaling"`, `"log"`, or `"auto"`.

- max_iter:

  Maximum solver iterations.

- tol:

  Requested certificate tolerance.

- check_every:

  Iteration interval for diagnostics.

- mass_policy:

  For balanced and partial problems, either `"probability"` (verify that
  both sums are one) or `"normalize"` (separately normalize and record
  the original masses). For KL-unbalanced problems, see
  [`transport_problem_unbalanced()`](https://bbuchsbaum.github.io/rfugw/reference/transport_problem_unbalanced.md).

## Value

An `rfugw_transport_problem` for
[`transport_solve()`](https://bbuchsbaum.github.io/rfugw/reference/transport_solve.md).
