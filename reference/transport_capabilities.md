# Query public transport-solver capabilities

Query public transport-solver capabilities

## Usage

``` r
transport_capabilities(problem = NULL)
```

## Arguments

- problem:

  Optional `rfugw_transport_problem`; when supplied, return the row for
  its estimand.

## Value

A data frame with one row per distinct estimand family. Comma- separated
fields enumerate methods and plan representations; they are descriptive,
not accepted as solver-selection strings. `maturity` and
`coverage_family` link every row to the executable numerical-path matrix
used by release evidence.
