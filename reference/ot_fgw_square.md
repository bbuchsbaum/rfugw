# Square-loss fused Gromov-Wasserstein objective

Unregularized FGW value:
`(1 - alpha) * <M, plan> + alpha * GW(C1, C2, plan)`.

## Usage

``` r
ot_fgw_square(M, C1, C2, plan, alpha = 0.5, symmetric = NULL)
```

## Arguments

- M:

  Feature cost matrix.

- C1:

  Source structure matrix.

- C2:

  Target structure matrix.

- plan:

  Coupling or `rfugw_result`.

- alpha:

  Feature/structure trade-off in `[0, 1]`.

- symmetric:

  `NULL` auto-detects; `TRUE` requires symmetry.

## Value

Numeric scalar.
