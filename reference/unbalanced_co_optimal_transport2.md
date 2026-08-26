# Unbalanced Co-Optimal Transport Objective Value

Unbalanced Co-Optimal Transport Objective Value

## Usage

``` r
unbalanced_co_optimal_transport2(..., log = FALSE)
```

## Arguments

- ...:

  Additional arguments. Unused extras are rejected when the solver uses
  `.reject_unused_dots()`; otherwise they are forwarded to the primary
  solver.

- log:

  If `TRUE`, return diagnostics and objective decomposition.

## Value

Numeric UCOOT objective value. If `log = TRUE`, returns a list with
`ucoot` and detailed diagnostics.
