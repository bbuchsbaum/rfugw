# Construct an explicit exact balanced transport problem

Construct an explicit exact balanced transport problem

## Usage

``` r
transport_problem_emd(
  M,
  p = NULL,
  q = NULL,
  max_iter = 20000L,
  tol = 1e-12,
  mass_policy = NULL
)
```

## Arguments

- M:

  Finite source-by-target cost matrix.

- p, q:

  Nonnegative source and target weights.

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
