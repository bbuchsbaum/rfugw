# Exact penalized variable-mass partial optimal transport

Minimizes
`<M,G> + discard_penalty * (sum(p) - sum(G)) + discard_penalty * (sum(q) - sum(G))`
over nonnegative plans with `rowSums(G) <= p` and `colSums(G) <= q`.
Thus a larger discard penalty weakly favors transporting more mass.
Unlike
[`ot_partial_emd()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_emd.md),
total transported mass is chosen by the optimization; unlike
KL-unbalanced OT, the marginal penalty is linear in discarded mass.

## Usage

``` r
ot_partial_penalized(
  M,
  p = NULL,
  q = NULL,
  discard_penalty,
  max_iter = 20000L,
  tol = 1e-12
)
```

## Arguments

- M:

  Finite nonnegative source-by-target transport cost.

- p:

  Source finite nonnegative measure (default uniform probability).

- q:

  Target finite nonnegative measure (default uniform probability).

- discard_penalty:

  Finite nonnegative cost per discarded unit on each marginal.
  Transporting one unit avoids two discard penalties.

- max_iter:

  Maximum augmented simplex iterations.

- tol:

  Augmented exact-solver tolerance.

## Value

An `rfugw_result` with the optimized `plan`, transported and discarded
masses, transport/penalty objective decomposition, and the scaled
primal-dual/reduced-cost certificate for the equivalent augmented
problem.

## Details

The implementation adds one dummy source and target. Real-to-dummy and
dummy-to-real edges cost exactly `discard_penalty`, while the
dummy-to-dummy edge costs zero. No hidden big-M or cost-derived penalty
is used.

## Examples

``` r
M <- matrix(c(0.2, 3, 2, 0.1), 2, 2)
out <- ot_partial_penalized(M, discard_penalty = 0.5)
out$transported_mass
#> [1] 1
rfugw_value(out)
#> [1] 0.15
```
