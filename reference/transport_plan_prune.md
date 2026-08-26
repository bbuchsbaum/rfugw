# Prune a Transport Plan and Invalidate Stale Certificates

Removes explicitly stored entries below `min_weight`, reports lost mass,
and returns a canonical edge-list plan. When pruning an `rfugw_result`,
any positive lost mass clears convergence/feasibility/objective
certificates and sets status to `"requires_revalidation"`.

## Usage

``` r
transport_plan_prune(plan, min_weight)
```

## Arguments

- plan:

  Explicit transport plan or `rfugw_result`.

- min_weight:

  Finite nonnegative retention threshold. Entries strictly below it are
  removed.

## Value

Pruned plan or result.
