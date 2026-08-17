# rfugw 0.0.1.9000

Development version toward the focused 0.1 release.

## Compatibility

- Public solver names, aliases, and diagnostic field meanings are specified
  in `inst/solver-contract.md`.
- Unregularized GW/FGW names mean conditional gradient with an exact
  linear-OT subproblem. They do not claim a global minimizer of the
  non-convex GW/FGW objective.

## Solver honesty

- UCOOT and across-spaces defaults now use a positive `epsilon` and succeed
  on documented minimal examples. Unsupported `l2` / `mm` / `lbfgsb` choices
  are rejected before computation.
- `G0` warm starts are consumed by the C++ exact-CG backend, not only by
  lpSolve.
- Flagship solvers auto-detect cost symmetry. `symmetric = TRUE` is validated
  rather than treated as a silent fast-path override.
- Results expose `status`, `converged`, `residual`, and marginal residuals.
  Iteration-limit and numerical-failure exits are distinct from convergence.
- Alpha, regularization, iteration limits, and non-finite inputs are
  validated consistently on the flagship solvers.

## Linear optimal transport

- Public primitives `ot_sinkhorn()`, `ot_emd()`, and
  `ot_sinkhorn_unbalanced()` wrap the C++ Sinkhorn and simplex backends and
  return `rfugw_result` objects. They are certified against analytic
  identities and POT 0.9.6.post1 fixtures (`inst/extdata/fixtures/linear_ot_fixture.json`).
- Independent helpers `ot_validate_plan()`, `ot_linear_cost()`,
  `ot_entropy()`, `ot_kl()`, `ot_gw_square()`, `ot_fgw_square()`, and
  `ot_barycentric_project()` check plans and recompute objectives.

## Results

- Flagship solvers return class `rfugw_result` with `print()`, `summary()`,
  and accessors `rfugw_plan()`, `rfugw_value()`, `rfugw_status()`, and
  `rfugw_residuals()`. Legacy list fields remain.

## Documentation

- `_pkgdown.yml` groups the reference by family and labels sampled / low-rank
  APIs as experimental. Flagship solvers have `\examples{}`. The solver
  guide compares formulation, regularization, maturity, and performance.
  `vignette("linear-ot")` is the task-oriented linear-OT guide.

## Build and optional dependencies

- Default compilation is portable (`-O2`, no `-march=native` / `-ffast-math`).
  Aggressive flags are opt-in via `RFUGW_FAST_FLAGS`. Windows `Makevars.win`
  is provided. OpenMP remains optional. Sanitizer and `-Werror` flags are
  opt-in via `RFUGW_EXTRA_CXXFLAGS` / `RFUGW_EXTRA_LIBS`.
- Hosted CI runs multi-OS `R CMD check`, serial/ASan smokes, generated-file
  consistency, and a pkgdown rebuild. The local one-command gate is
  `Rscript tools/release-gate.R`.
- Speed evidence uses `inst/bench/PROTOCOL.md`. The entry point
  `Rscript inst/bench/run_protocol.R` records environment metadata and
  rejects invalid-quality runs.
- Unregularized partial solvers document and test the `lpSolve` Suggests
  route. Entropic partial and default `cpp_transport` exact CG do not need it.
  Exact and entropic partial GW/FGW now run their outer loops in C++. The
  inner entropic projection keeps a finite plan at high transported mass.
- UCOOT and across-spaces KL Sinkhorn now run their BCD loop in C++ with
  reusable sample/feature workspaces and inner warm-start telemetry.
- Sampled and low-rank GW stay experimental. The certified envelope is
  tiny-versus-full budget quality, rank reconstruction, and input-memory
  scaling only; unusable budgets and ranks warn or error instead of
  silently clamping. See `inst/bench/sampled-budget-curves.md`.
- Threading is certified only on the batched multiset kernels:
  1-thread and N-thread plans match, BLAS is pinned to avoid
  oversubscription, and `structure_knn` is documented as a metric
  change rather than a memory saving. See
  `inst/bench/threading-memory.md`.

## Metadata

- Package authorship, license holder, repository URL, bug reports, and
  citation metadata now identify the real maintainer.
- `CONTRIBUTING.md` records the compatibility, test, benchmark, generated
  documentation, and evidence policy for 0.1.
