# Generalized KL divergence of a plan from a product reference

Generalized KL divergence of a plan from a product reference

## Usage

``` r
ot_kl(plan, p, q)
```

## Arguments

- plan:

  Coupling or `rfugw_result`.

- p:

  Finite nonnegative source reference weights.

- q:

  Finite nonnegative target reference weights.

## Value

Generalized `KL(plan || p %o% q)`, including the mass correction
`-sum(plan) + sum(p %o% q)`. Positive plan mass outside zero reference
support returns `Inf`; zero-over-zero contributes zero.
