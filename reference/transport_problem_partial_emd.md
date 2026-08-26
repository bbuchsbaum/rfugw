# Construct an explicit exact fixed-mass partial transport problem

Construct an explicit exact fixed-mass partial transport problem

## Usage

``` r
transport_problem_partial_emd(
  M,
  p = NULL,
  q = NULL,
  mass = 1,
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

- mass:

  Transported probability mass in `[0, 1]`.

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
