# Construct an explicit KL-unbalanced transport problem

Unlike the probability-only problem constructors, this constructor can
preserve finite measures. The required `mass_policy` maps to the direct
solver without a hidden default: `"finite_measure"` preserves both
measures, `"joint"` applies one common scale, and
`"separate_probability"` normalizes them separately.

## Usage

``` r
transport_problem_unbalanced(
  M,
  p = NULL,
  q = NULL,
  epsilon = 0.05,
  rho = 10,
  method = c("auto", "scaling", "log"),
  max_iter = 500L,
  tol = 1e-07,
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

- rho:

  One or two positive marginal-KL penalties, in cost units.

- method:

  Requested Sinkhorn method: `"scaling"`, `"log"`, or `"auto"`.

- max_iter:

  Maximum solver iterations.

- tol:

  Requested certificate tolerance.

- mass_policy:

  One of `"finite_measure"`, `"joint"`, or `"separate_probability"`.

## Value

An `rfugw_transport_problem` for
[`transport_solve()`](https://bbuchsbaum.github.io/rfugw/reference/transport_solve.md).
