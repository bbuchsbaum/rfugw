# Extract reusable solver state

Returned state is opaque, versioned, serializable with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html), and certified for
reuse only when the originating solve converged. Balanced dual-potential
state supports scaling/log interoperability and epsilon continuation.
Partial-Sinkhorn Dykstra state is intentionally bound to the identical
problem and effective method.

## Usage

``` r
rfugw_state(x)
```

## Arguments

- x:

  An `rfugw_result` or `rfugw_solver_state`.

## Value

An `rfugw_solver_state`.
