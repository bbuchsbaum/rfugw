# rfugw 0.0.1.9000

Development version toward the focused 0.1 release.

## Compatibility

- Public solver names, aliases, and diagnostic field meanings are specified
  in `inst/solver-contract.md`.
- Unregularized GW/FGW names mean conditional gradient with an exact
  linear-OT subproblem. They do not claim a global minimizer of the
  non-convex GW/FGW objective.

## Solver honesty

- An executable numerical path matrix now maps every advertised symmetry,
  regularization, Sinkhorn, precision, warm-start, adapter, threading, and
  approximation path to its comparison law and CI scope.
- A bounded mutation-proof harness recreates and kills the reviewed transpose,
  false-simplex-success, tolerance-floor, ignored-inner-status, missing-KL-mass,
  and zero-support-KL defect classes without modifying shipping source.
- Distinct PR, nightly, and release numerical-trust workflows now publish
  replayable ledgers. The release matrix checks and installs the exact tarball
  on Linux, macOS, and Windows, with a separate ASan/UBSan gate; evidence files
  keep hosted success separate from publication status.
- Native batch kernels now copy R inputs into owned C++ storage before OpenMP,
  use only read-only shared caches and per-job writes, capture worker exceptions
  for main-thread reporting, suppress nested OpenMP, and report requested,
  effective, and maximum thread counts. A pure-C++ ThreadSanitizer harness and
  1/2/4-thread differential tests guard the contract.
- Canonical square-loss GW algebra and batch-thread adapters now have dedicated
  translation units; transport certificates and approximation caches have
  R-independent headers, with float/double cache invariants implemented once as
  templates. `inst/native-architecture.md` records the remaining solver-engine
  boundary and the measured compile/object-size cost of the split.
- Benchmark evidence now distinguishes certified comparisons from fixed-budget
  performance sentinels. `max_iter` and `experimental` rows can never become
  convergence claims; quality, setup/solve/end-to-end time, allocations,
  precision, threads, commit, and environment are separate fields/artifacts.
  Sampled GW records exact-reference error at three budgets, and threshold
  changes require a digest-matched retained evidence-history entry.
- Public convergence is now certificate-based across balanced, partial, and
  semirelaxed, and unbalanced families: returned plans must be finite and nonnegative,
  formulation feasibility must pass, required inner solves must be certified,
  and reported objectives and named components must agree with independent
  recomputation. Failures report `infeasible`, `objective_mismatch`,
  `inner_failure`, or `numerical_failure`; iteration count alone cannot imply
  success.
- Precision-selecting solvers report requested/effective tolerances,
  requested/effective/actual compute precision, backend transitions, and a
  termination reason. Tight mixed requests are promoted to strict double
  instead of silently flooring `1e-9` requests to `1e-6`; explicit
  `strict_double` blocks the large-problem float acceleration path.
  Balanced entropic FGW now defaults to an explicit dynamic-range `auto`
  Sinkhorn policy, selecting scaling only within the documented exponent
  regime and otherwise selecting and reporting the robust log-domain backend.
  Linear KL-unbalanced OT has a genuine log-domain implementation; the former
  UCOOT `sinkhorn_log` scaling alias now errors as deprecated. Entropic partial
  GW/FGW explicitly reject unsafe scaling and log requests until a genuine
  log-domain Dykstra backend exists.
- Nested exact and entropic solvers now report final and maximum inner
  residuals, total inner iterations, and `inner_converged`. Exact partial
  directions preserve simplex certificates instead of extracting only the
  plan, and an uncertified required inner solve can no longer imply outer
  convergence. Deterministic nested fault injection covers iteration-limit,
  numerical, infeasible, and non-finite failures.
- Asymmetric exact and entropic partial GW/FGW were first quarantined, then
  re-enabled after the target-marginal transpose defect was fixed in a shared
  square-loss tensor module. O(n^4) gradient, directional-derivative,
  permutation, exact-LP, and entropic-projection oracles now guard the path.
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
- Fractional iteration/rank/budget counts are rejected instead of truncated;
  partial solvers record whether transported mass was explicit or defaulted;
  symmetry detection now uses a scale-aware absolute-plus-relative tolerance.
- `ot_kl()` now implements generalized KL, including unequal-mass correction,
  zero-reference support rules, and finite nonnegative reference validation.

## Linear optimal transport

- Exact transport now returns explicit termination reasons, dual potentials,
  primal and dual objectives, marginal residuals, minimum nonbasic reduced
  cost, duality gap, and the tolerances used to certify them. Only a feasible,
  dual-feasible, zero-gap result is reported as optimal; deterministic fault
  injection covers every internal termination branch.
- Public primitives `ot_sinkhorn()`, `ot_emd()`, and
  `ot_sinkhorn_unbalanced()` wrap the C++ Sinkhorn and simplex backends and
  return `rfugw_result` objects. They are certified against analytic
  identities and POT 0.9.6.post1 fixtures (`inst/extdata/fixtures/linear_ot_fixture.json`).
- `ot_partial_emd()` now exposes exact nonnegative-cost partial linear OT by
  reducing it to the existing certified simplex with one dummy source/target.
  It reports partial mass/feasibility plus the augmented primal/dual
  certificate and is guarded by analytic, independent-LP, metamorphic,
  adversarial, regression, and local baseline evidence. Entropic partial OT,
  Sinkhorn divergence, and fixed-support Wasserstein weight barycenters remain
  explicitly deferred with admission criteria in
  `inst/linear-ot-foundations.md`.
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
- Hosted flagship gates run the canonical protocol with quality
  checks first. The PR gate covers FGW, FUGW, and semirelaxed. Nightly
  adds partial, UCOOT, sampled, larger sizes, and a 1-vs-2 thread
  smoke. Time caps are median `solve_ms` slack, not tight speed
  claims. Exit 1 is a solver regression; exit 2 is infrastructure.

## Metadata

- Package authorship, license holder, repository URL, bug reports, and
  citation metadata now identify the real maintainer.
- `CONTRIBUTING.md` records the compatibility, test, benchmark, generated
  documentation, and evidence policy for 0.1.
