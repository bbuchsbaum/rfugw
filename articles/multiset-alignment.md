# Multiset alignment

Multiset alignment puts several attributed structures into
correspondence with one shared template. A subject might be a brain
parcellation, graph, or point cloud; subjects may have different numbers
of nodes.
[`multialign_fit()`](https://bbuchsbaum.github.io/rfugw/reference/multialign_fit.md)
returns one subject-to-template coupling per subject, plus the template
and solver diagnostics.

This vignette aligns three related point clouds, checks their couplings,
then uses those couplings to move features and square matrices between
subject and template coordinates.

``` r

library(rfugw)
```

## What does a subject look like?

Each subject is a list containing a square structure matrix `C`. It may
also contain node features `F`, probability weights `w`, and a character
`id`. Features must use the same columns across subjects.

``` r

make_subject <- function(n, phase, seed) {
  set.seed(seed)
  angle <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)] + phase
  xy <- cbind(cos(angle), 0.65 * sin(angle))
  xy <- xy + matrix(rnorm(2 * n, sd = 0.04), ncol = 2)
  C <- as.matrix(dist(xy))
  F <- cbind(horizontal = xy[, 1], vertical = xy[, 2])
  list(C = C / max(C), F = F, id = paste0("subject_", seed))
}

subjects <- list(
  make_subject(12, 0.00, 1),
  make_subject(15, 0.12, 2),
  make_subject(10, -0.10, 3)
)
```

``` r

data.frame(
  id = vapply(subjects, `[[`, character(1), "id"),
  nodes = vapply(subjects, function(x) nrow(x$C), integer(1)),
  features = vapply(subjects, function(x) ncol(x$F), integer(1))
)
#>          id nodes features
#> 1 subject_1    12        2
#> 2 subject_2    15        2
#> 3 subject_3    10        2
```

Missing `w` means uniform node weights. A missing `F` is allowed only
when a compatible template is supplied; automatic template construction
pools the subject features.

## How do you align to a shared template?

With `template = NULL`,
[`multialign_fit()`](https://bbuchsbaum.github.io/rfugw/reference/multialign_fit.md)
builds a template by clustering the pooled features. `k_template`
chooses its number of nodes. The default `fgw_entropic` backend balances
feature similarity and structural similarity.

``` r

fit <- multialign_fit(
  subjects,
  k_template = 8,
  alpha = 0.5,
  epsilon = 0.08,
  max_iter = 50L,
  sinkhorn_max_iter = 250L,
  tol = 1e-9
)
```

``` r

data.frame(
  id = names(fit$couplings),
  subject_nodes = vapply(fit$couplings, nrow, integer(1)),
  template_nodes = vapply(fit$couplings, ncol, integer(1)),
  status = fit$diagnostics$status,
  converged = fit$diagnostics$converged,
  inner_status = fit$diagnostics$inner_status,
  objective = fit$objectives
)
#>                  id subject_nodes template_nodes
#> subject_1 subject_1            12              8
#> subject_2 subject_2            15              8
#> subject_3 subject_3            10              8
#>                                      status converged inner_status  objective
#> subject_1 outer_converged_inner_unavailable        NA  unavailable 0.09128723
#> subject_2 outer_converged_inner_unavailable        NA  unavailable 0.09563146
#> subject_3 outer_converged_inner_unavailable        NA  unavailable 0.07962000
```

Each coupling has one row per subject node and one column per template
node. For balanced FGW, its row and column sums reproduce the subject
and template weights. Different node counts do not by themselves require
an unbalanced solver. This run uses the native batch path: the outer
stopping residual and balanced feasibility checks pass, but the batch
kernel does not return nested Sinkhorn certificates. Consequently
`converged` is `NA` and `inner_status` is `"unavailable"`; the result
does not claim fully certified convergence.

## How do you reuse a template?

Build a template explicitly when several analyses should use the same
target, or when you want to inspect and save it before fitting.

``` r

template <- multialign_make_template(subjects, k = 8, seed = 42)
fit_fixed <- multialign_fit(
  subjects, template = template,
  alpha = 0.5, epsilon = 0.08,
  max_iter = 50L, sinkhorn_max_iter = 250L
)
```

``` r

fit_fixed$diagnostics[, c(
  "id", "status", "converged", "inner_status", "objective"
)]
#>                  id                            status converged inner_status
#> subject_1 subject_1 outer_converged_inner_unavailable        NA  unavailable
#> subject_2 subject_2 outer_converged_inner_unavailable        NA  unavailable
#> subject_3 subject_3 outer_converged_inner_unavailable        NA  unavailable
#>            objective
#> subject_1 0.09241252
#> subject_2 0.08997615
#> subject_3 0.08276241
```

An explicit template can also be any valid list with `C`, `F`, and `w`.
Reusing one makes coordinates comparable across separate runs; an
automatically built template is specific to the subjects supplied in
that call.

## How do you project template features to a subject?

[`multialign_project_features()`](https://bbuchsbaum.github.io/rfugw/reference/multialign_project_features.md)
maps a template-level matrix back to subject resolution. Row
normalization turns each subject row of the coupling into weights over
template nodes.

``` r

P1 <- fit$couplings$subject_1
subject_1_features <- multialign_project_features(P1, fit$template$F)
data.frame(
  feature = colnames(subject_1_features),
  correlation = diag(cor(subjects[[1]]$F, subject_1_features))
)
#>               feature correlation
#> horizontal horizontal   0.9971606
#> vertical     vertical   0.9922839
```

The output has one row per subject node and the same columns as the
template features. The high correlations show that this synthetic
alignment preserves the known feature pattern. Because the template was
built from these same features, this is an illustrative sanity check,
not independent accuracy evidence.

## How do you move a square matrix to template space?

[`multialign_project_matrix()`](https://bbuchsbaum.github.io/rfugw/reference/multialign_project_matrix.md)
pushes a subject-level square matrix through both sides of the coupling.
This is useful for connectivity or similarity matrices.

``` r

template_structure <- multialign_project_matrix(
  P1, subjects[[1]]$C, normalize = "col"
)
dim(template_structure)
#> [1] 8 8
```

Normalization changes the interpretation. Use `"row"` when interpolating
template values at subject nodes and `"col"` when aggregating subject
values at template nodes; use
[`multialign_normalize_plan()`](https://bbuchsbaum.github.io/rfugw/reference/multialign_normalize_plan.md)
directly when you need to inspect those weights.

## When should you use unbalanced FUGW?

The default `method = "fgw_entropic"` fixes both marginals. Choose
`method = "fugw_kl"` when unmatched regions, contamination, or genuine
mass variation make exact marginal preservation inappropriate.

``` r

fit_unbalanced <- multialign_fit(
  subjects,
  template = template,
  method = "fugw_kl",
  alpha = 0.5,
  epsilon = 0.5,
  reg_marginals = c(1, 1),
  max_iter = 200L,
  tol = 1e-6,
  max_iter_ot = 1000L,
  tol_ot = 1e-6,
  use_cpp_batch = FALSE
)
data.frame(
  id = names(fit_unbalanced$couplings),
  status = fit_unbalanced$diagnostics$status,
  inner_status = fit_unbalanced$diagnostics$inner_status,
  transported_mass = fit_unbalanced$diagnostics$plan_mass
)
#>          id    status inner_status transported_mass
#> 1 subject_1 converged    converged        0.9401245
#> 2 subject_2 converged    converged        0.9425568
#> 3 subject_3 converged    converged        0.9403135
```

Unbalanced row and column sums are intentionally relaxed. The
`reg_marginals` penalties control that relaxation, so transported mass
and solver diagnostics should be checked rather than judged by balanced
marginal equalities.

## Should the template be fixed or learned?

The default `template_mode = "fixed"` estimates couplings to an
unchanged template. `template_mode = "learned"` alternates alignment
with barycenter-style template updates.

``` r

learned <- multialign_fit(
  subjects,
  template = template,
  template_mode = "learned",
  alpha = 0.5,
  epsilon = 0.08,
  max_iter = 40L,
  template_max_iter = 3L
)
last_template_step <- learned$template_history[learned$outer_iterations, ]
data.frame(
  outer_iterations = learned$outer_iterations,
  stopped_at_cap = learned$outer_iterations == 3L,
  template_tolerance_met =
    last_template_step$delta_template < 1e-6 &&
    last_template_step$delta_objective < 1e-6
)
#>   outer_iterations stopped_at_cap template_tolerance_met
#> 1                3           TRUE                  FALSE
```

A learned template is a fitted summary of the supplied collection rather
than a fixed external reference. It costs additional alignment passes
and makes the result specific to those inputs. In this example the run
reaches its three-iteration cap; it must not be described as converged
unless both reported deltas meet `template_tol`. Even when they do, the
non-convex optimization does not establish a unique or globally optimal
template.

## What should you tune first?

- `alpha` controls structure (`1`) versus features (`0`). Choose it from
  the scientific question or downstream validation.
- `epsilon` controls entropic smoothing. Smaller values can sharpen
  couplings but are harder to optimize.
- `precision`, batching, low-rank structure, and `autotune` are
  performance controls. Establish a trustworthy full-precision result
  before relying on approximations, and compare both objectives and
  coupling feasibility.

## Where should you go next?

- [`vignette("barycenters")`](https://bbuchsbaum.github.io/rfugw/articles/barycenters.md)
  explains direct fixed-support GW and FGW barycenters.
- [`vignette("solver-guide")`](https://bbuchsbaum.github.io/rfugw/articles/solver-guide.md)
  covers solver formulations, convergence, and diagnostics.
- See
  [`?multialign_fit`](https://bbuchsbaum.github.io/rfugw/reference/multialign_fit.md)
  and
  [`?multialign_project_features`](https://bbuchsbaum.github.io/rfugw/reference/multialign_project_features.md)
  for the complete interfaces.
