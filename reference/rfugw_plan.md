# Extract the coupling from an rfugw result

Extract the coupling from an rfugw result

## Usage

``` r
rfugw_plan(x, materialize = FALSE)
```

## Arguments

- x:

  An `rfugw_result`, transport plan, or a list with `plan` / `pi_samp`.

- materialize:

  If `TRUE`, explicitly convert a transport-plan representation to a
  dense matrix. Dense legacy results are unchanged.

## Value

The stored coupling matrix/representation, or a dense matrix when
explicitly requested.
