# Solve a versioned transport problem

`transport_solve()` validates the complete problem and any reusable
state before entering a numerical backend. When `state_policy = "cold"`,
an incompatible state is rejected with a structured reason and the
problem is solved cold; `"error"` fails closed instead. A state is never
inherited as evidence that the new result converged.

## Usage

``` r
transport_solve(problem, init_state = NULL, state_policy = c("error", "cold"))
```

## Arguments

- problem:

  An `rfugw_transport_problem`.

- init_state:

  Optional state returned by
  [`rfugw_state()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_state.md)
  or a certified prior `rfugw_result`.

- state_policy:

  What to do with invalid or unsupported state: error or explicitly
  record rejection and solve cold.

## Value

An `rfugw_result` with protocol, problem, initialization, mass, control,
plan-representation, status, and runtime provenance.
