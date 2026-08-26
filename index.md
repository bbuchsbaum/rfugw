# rfugw

[![R CMD
check](https://github.com/bbuchsbaum/rfugw/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/bbuchsbaum/rfugw/actions/workflows/R-CMD-check.yml)
[![Numerical
trust](https://github.com/bbuchsbaum/rfugw/actions/workflows/numerical-trust.yml/badge.svg)](https://github.com/bbuchsbaum/rfugw/actions/workflows/numerical-trust.yml)
[![pkgdown](https://github.com/bbuchsbaum/rfugw/actions/workflows/pkgdown.yml/badge.svg)](https://github.com/bbuchsbaum/rfugw/actions/workflows/pkgdown.yml)

[Documentation](https://bbuchsbaum.github.io/rfugw/) · [Getting
started](https://bbuchsbaum.github.io/rfugw/articles/rfugw.html) ·
[Solver-client
protocol](https://bbuchsbaum.github.io/rfugw/articles/solver-contract.html)
· [Solver
guide](https://bbuchsbaum.github.io/rfugw/articles/solver-guide.html) ·
[Function reference](https://bbuchsbaum.github.io/rfugw/reference/) ·
[Changelog](https://bbuchsbaum.github.io/rfugw/NEWS.md)

`rfugw` is an R package for finding soft correspondences between
attributed collections when their coordinate systems are not directly
comparable. It provides linear optimal transport and Gromov–Wasserstein
solvers backed by C++ kernels, with explicit convergence and feasibility
diagnostics.

Use GW when only within-collection geometry is comparable, FGW when the
collections also share features, and partial or unbalanced formulations
when not all mass should be matched.

> **Status:** Version 0.1.0 defines the first supported solver and
> diagnostic contract. Sampled GW and dense-plan SVD remain explicitly
> experimental, and no end-to-end scalable relational-OT path is
> claimed.

## Install

Install the development version from GitHub:

``` r

install.packages("remotes")
remotes::install_github("bbuchsbaum/rfugw")
```

Installation requires R and a C++17 toolchain. The package has no Python
runtime dependency; OpenMP is optional.

## Quick start

This example aligns a point set with a rotated and reordered copy.
Within-set distances describe structure, while a shared scalar feature
disambiguates the objects.

``` r

library(rfugw)

source_xy <- rbind(A = c(0, 0), B = c(1, 0), C = c(1.4, 0.8),
                   D = c(0.7, 1.5), E = c(-0.3, 0.9), F = c(0.2, 0.4))
target_order <- c(4, 1, 6, 3, 5, 2)
theta <- pi / 3
rotation <- matrix(c(cos(theta), -sin(theta),
                     sin(theta),  cos(theta)), 2, 2)
target_xy <- source_xy[target_order, ] %*% rotation
rownames(target_xy) <- rownames(source_xy)[target_order]

C1 <- as.matrix(dist(source_xy)); C1 <- C1 / max(C1)
C2 <- as.matrix(dist(target_xy)); C2 <- C2 / max(C2)
feature <- c(0.05, 0.25, 0.50, 0.70, 0.90, 0.38)
M <- abs(outer(feature, feature[target_order], "-"))

fit <- fgw_entropic(
  M, C1, C2, alpha = 0.5,
  epsilon = 0.08, sinkhorn_tol = 1e-8
)
plan <- rfugw_plan(fit)
data.frame(
  source = rownames(source_xy),
  target = rownames(target_xy)[max.col(plan, ties.method = "first")]
)
#>   source target
#> 1      A      A
#> 2      B      B
#> 3      C      C
#> 4      D      D
#> 5      E      E
#> 6      F      F
```

The transport plan is a nonnegative matrix of soft correspondences, not
just a set of hard matches. Before interpreting it, check the solver
status and its diagnostics:

``` r

diagnostics <- rfugw_residuals(fit)
identical(rfugw_status(fit), "converged") &&
  isTRUE(diagnostics$feasible) &&
  isTRUE(diagnostics$inner_converged)
#> [1] TRUE
```

## What it covers

- Solve balanced, exact, partial, and KL-unbalanced linear transport
  problems.
- Align metric spaces with GW, FGW, partial, semirelaxed, and FUGW
  models.
- Build fixed-size GW/FGW barycenters with learned structure and
  features.
- Optimize ordinary Wasserstein barycenter weights on a fixed common
  support, with exact or explicitly regularized objectives and KKT
  diagnostics.
- Align several attributed metric spaces to a shared template.
- Inspect plans, objectives, stopping status, and residuals through
  stable result accessors.
- Build explicit, versioned transport problems and reuse certified
  iterative state through the downstream solver-client protocol.

## Choosing a solver

| Data and matching assumption | Start with |
|----|----|
| A direct source-to-target cost is meaningful | [`ot_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn.md) |
| Only within-domain structure is comparable | [`entropic_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_gromov_wasserstein.md) |
| Structure and shared features both matter | [`fgw_entropic()`](https://bbuchsbaum.github.io/rfugw/reference/fgw_entropic.md) |
| Some mass should remain unmatched | a partial solver |
| Marginals may change with a soft penalty | [`fugw_kl()`](https://bbuchsbaum.github.io/rfugw/reference/fugw_kl.md) or an unbalanced solver |
| Support locations are fixed; only barycenter weights are unknown | [`ot_barycenter_weights()`](https://bbuchsbaum.github.io/rfugw/reference/ot_barycenter_weights.md) |

The [solver
guide](https://bbuchsbaum.github.io/rfugw/articles/solver-guide.html)
explains the modeling assumptions, regularization, and trade-offs in
detail.

## Fit and boundaries

`rfugw` focuses on square-loss optimal-transport formulations rather
than the full Python Optimal Transport (POT) API. Unregularized GW and
FGW use conditional-gradient procedures with exact linear-transport
directions; the outer problems remain non-convex, so a converged result
is not a certificate of the global minimum.

Sampled GW functions and
[`dense_gromov_wasserstein_plan_svd()`](https://bbuchsbaum.github.io/rfugw/reference/dense_gromov_wasserstein_plan_svd.md)
are experimental. The latter compresses an already-materialized dense GW
plan; it is not an end-to-end low-rank or low-memory solver. Its former
name,
[`lowrank_gromov_wasserstein_samples()`](https://bbuchsbaum.github.io/rfugw/reference/lowrank_gromov_wasserstein_samples.md),
is deprecated. Coordinate and sparse- graph sampled inputs also return
and optimize a dense coupling, so they are not supported end-to-end
scalable paths. The measured [admission
decision](https://bbuchsbaum.github.io/rfugw/inst/scalable-relational-ot-decision.md)
records the quadratic boundary and the quality, memory, and operator
gates required for future promotion. The default exact FGW/GW and
partial linear-transport paths do not require `lpSolve`; only the
optional `lp_transport` and `lp_matrix` backends do. See the [solver
contract](https://bbuchsbaum.github.io/rfugw/inst/solver-contract.md)
for supported formulations, diagnostic meanings, and compatibility
commitments.

## Documentation

- [Getting
  started](https://bbuchsbaum.github.io/rfugw/articles/rfugw.html) —
  build and interpret a first FGW alignment.
- [Linear optimal
  transport](https://bbuchsbaum.github.io/rfugw/articles/linear-ot.html)
  — exact, entropic, partial, and unbalanced transport.
- [Solver-client
  protocol](https://bbuchsbaum.github.io/rfugw/articles/solver-contract.html)
  — explicit estimands, mass policy, reusable state, capabilities, and
  provenance for downstream libraries.
- [Choosing a
  solver](https://bbuchsbaum.github.io/rfugw/articles/solver-guide.html)
  — select a formulation and interpret its diagnostics.
- [Barycenters](https://bbuchsbaum.github.io/rfugw/articles/barycenters.html)
  — construct representative metric spaces.
- [Multiset
  alignment](https://bbuchsbaum.github.io/rfugw/articles/multiset-alignment.html)
  — align several collections to a template.
- [Function reference](https://bbuchsbaum.github.io/rfugw/reference/) —
  generated help pages for every exported function.
- [Numerical trust
  charter](https://bbuchsbaum.github.io/rfugw/inst/numerical-trust-charter.md)
  and [benchmark
  protocol](https://bbuchsbaum.github.io/rfugw/inst/bench/PROTOCOL.md) —
  correctness and performance evidence requirements.

## Contributing

See
[CONTRIBUTING.md](https://bbuchsbaum.github.io/rfugw/CONTRIBUTING.md)
for development checks, numerical-test requirements, and benchmark
policy.

## Citation

``` r

citation("rfugw")
```

## License

MIT © 2026 Bradley Buchsbaum.
