# rfugw

[![R CMD check](https://github.com/bbuchsbaum/rfugw/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/bbuchsbaum/rfugw/actions/workflows/R-CMD-check.yml)
[![Numerical trust](https://github.com/bbuchsbaum/rfugw/actions/workflows/numerical-trust.yml/badge.svg)](https://github.com/bbuchsbaum/rfugw/actions/workflows/numerical-trust.yml)
[![pkgdown](https://github.com/bbuchsbaum/rfugw/actions/workflows/pkgdown.yml/badge.svg)](https://github.com/bbuchsbaum/rfugw/actions/workflows/pkgdown.yml)

[Getting started](vignettes/rfugw.Rmd) ·
[Solver guide](vignettes/solver-guide.Rmd) ·
[Function reference](man/) ·
[Changelog](NEWS.md)

`rfugw` is an R package for finding soft correspondences between attributed
collections when their coordinate systems are not directly comparable. It
provides linear optimal transport and Gromov--Wasserstein solvers backed by C++
kernels, with explicit convergence and feasibility diagnostics.

Use GW when only within-collection geometry is comparable, FGW when the
collections also share features, and partial or unbalanced formulations when
not all mass should be matched.

> **Status:** Pre-release and currently installed from source. The project is
> working toward 0.1, and APIs may still change before that release.

## Install

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("bbuchsbaum/rfugw")
```

Installation requires R and a C++17 toolchain. The package has no Python
runtime dependency; OpenMP is optional.

## Quick start

This example aligns a point set with a rotated and reordered copy. Within-set
distances describe structure, while a shared scalar feature disambiguates the
objects.

```r
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

The transport plan is a nonnegative matrix of soft correspondences, not just a
set of hard matches. Before interpreting it, check the solver status and its
diagnostics:

```r
diagnostics <- rfugw_residuals(fit)
identical(rfugw_status(fit), "converged") &&
  isTRUE(diagnostics$feasible) &&
  isTRUE(diagnostics$inner_converged)
#> [1] TRUE
```

## What it covers

- Solve balanced, exact, partial, and KL-unbalanced linear transport problems.
- Align metric spaces with GW, FGW, partial, semirelaxed, and FUGW models.
- Build fixed-size GW/FGW barycenters with learned structure and features.
- Align several attributed metric spaces to a shared template.
- Inspect plans, objectives, stopping status, and residuals through stable
  result accessors.

## Choosing a solver

| Data and matching assumption | Start with |
|---|---|
| A direct source-to-target cost is meaningful | `ot_sinkhorn()` |
| Only within-domain structure is comparable | `entropic_gromov_wasserstein()` |
| Structure and shared features both matter | `fgw_entropic()` |
| Some mass should remain unmatched | a partial solver |
| Marginals may change with a soft penalty | `fugw_kl()` or an unbalanced solver |

The [solver guide](vignettes/solver-guide.Rmd) explains the modeling
assumptions, regularization, and trade-offs in detail.

## Fit and boundaries

`rfugw` focuses on square-loss optimal-transport formulations rather than the
full Python Optimal Transport (POT) API. Unregularized GW and FGW use
conditional-gradient procedures with exact linear-transport directions; the
outer problems remain non-convex, so a converged result is not a certificate of
the global minimum.

Sampled and low-rank GW functions are experimental. The default exact FGW/GW
and partial linear-transport paths do not require `lpSolve`; only the optional
`lp_transport` and `lp_matrix` backends do. See the
[solver contract](inst/solver-contract.md) for supported formulations,
diagnostic meanings, and compatibility commitments.

## Documentation

- [Getting started](vignettes/rfugw.Rmd) — build and interpret a first FGW
  alignment.
- [Linear optimal transport](vignettes/linear-ot.Rmd) — exact, entropic,
  partial, and unbalanced transport.
- [Barycenters](vignettes/barycenters.Rmd) — construct representative metric
  spaces.
- [Multiset alignment](vignettes/multiset-alignment.Rmd) — align several
  collections to a template.
- [Function reference](man/) — generated help pages for every exported
  function.
- [Numerical trust charter](inst/numerical-trust-charter.md) and
  [benchmark protocol](inst/bench/PROTOCOL.md) — correctness and performance
  evidence requirements.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development checks, numerical-test
requirements, and benchmark policy.

## Citation

```r
citation("rfugw")
```

## License

MIT © 2026 Bradley Buchsbaum.
