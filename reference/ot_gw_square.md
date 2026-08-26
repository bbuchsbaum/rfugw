# Square-loss Gromov-Wasserstein objective

Independent evaluator of `sum_ijkl (C1_ij - C2_kl)^2 plan_ik plan_jl`.
Uses the standard factored expansion; tests compare it to the O(n^4)
form.

## Usage

``` r
ot_gw_square(C1, C2, plan, symmetric = NULL)
```

## Arguments

- C1:

  Source structure matrix.

- C2:

  Target structure matrix.

- plan:

  Coupling or `rfugw_result`.

- symmetric:

  `NULL` auto-detects; `TRUE` requires symmetry.

## Value

Numeric scalar.
