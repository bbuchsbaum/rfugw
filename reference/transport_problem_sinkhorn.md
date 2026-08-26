# Construct an explicit balanced entropic transport problem

These problem constructors are the versioned client boundary for
[`transport_solve()`](https://bbuchsbaum.github.io/rfugw/reference/transport_solve.md).
Rows of `M` are source support and columns are target support.
`mass_policy` is mandatory so a downstream caller cannot accidentally
confuse already-normalized probabilities with weights that it intends
rfugw to normalize.

## Usage

``` r
transport_problem_sinkhorn(
  M,
  p = NULL,
  q = NULL,
  epsilon = 0.05,
  method = c("auto", "scaling", "log"),
  max_iter = 1000L,
  tol = 1e-09,
  mass_policy = NULL
)
```

## Arguments

- M:

  Finite source-by-target cost matrix.

- p, q:

  Nonnegative source and target weights.

- epsilon:

  Positive coefficient of `KL(plan || p %o% q)`, in the same units as
  `M`.

- method:

  Requested Sinkhorn method: `"scaling"`, `"log"`, or `"auto"`.

- max_iter:

  Maximum solver iterations.

- tol:

  Requested certificate tolerance.

- mass_policy:

  For balanced and partial problems, either `"probability"` (verify that
  both sums are one) or `"normalize"` (separately normalize and record
  the original masses). For KL-unbalanced problems, see
  [`transport_problem_unbalanced()`](https://bbuchsbaum.github.io/rfugw/reference/transport_problem_unbalanced.md).

## Value

An `rfugw_transport_problem` for
[`transport_solve()`](https://bbuchsbaum.github.io/rfugw/reference/transport_solve.md).
