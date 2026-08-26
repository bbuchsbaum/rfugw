# Changelog

## rfugw 0.1.0

Initial focused release.

- The scalable relational-OT admission review closes with a measured
  deferral: no sampled, graph, operator, or low-rank candidate currently
  avoids both dense structure and dense coupling materialization while
  meeting a certified quality envelope. Coordinate inputs reduce input
  storage, but the native sampled kernel retains at least five dense
  coupling-sized matrices and its measured 400-to-800 runtime slope is
  1.936. The experimental labels remain;
  `inst/scalable-relational-ot-decision.md` records target sizes,
  rejected alternatives, exact workspace lower bounds, and promotion
  thresholds.

- [`ot_barycenter_weights()`](https://bbuchsbaum.github.io/rfugw/reference/ot_barycenter_weights.md)
  now optimizes probability weights on a fixed user-supplied support
  without relabeling a GW/FGW support-learning update. Exact mode solves
  one joint LP and re-certifies every component with exact EMD.
  Regularized mode minimizes a semi-debiased product-reference-KL
  Sinkhorn objective with analytic dual gradients, safeguarded descent,
  projected-simplex KKT, a declared numerical weight floor, and reusable
  component state. Independent LP/convex optimizers, finite differences,
  identity/permutation laws, zero weights, duplicate support, imbalanced
  coefficients, tiny entropy, warm reuse, and retained solve/allocation
  evidence guard the two distinct objective contracts.

- A versioned solver-client protocol now lets downstream packages
  construct explicit balanced entropic, exact balanced, exact/entropic
  fixed-mass partial, and finite-measure KL-unbalanced problems, solve
  them through one validated boundary, and consume stable problem, mass,
  provenance, plan, value, status, residual, and state accessors. Mass
  policy is mandatory and formulation-specific. Balanced dual state
  supports scaling/log/auto and epsilon continuation; partial Dykstra
  state remains bound to the identical problem/backend. Structured
  fail-closed rejection, cold fallback, RDS/fresh- session reuse, a
  capability matrix, a public-only generic client, live DKGE
  differential tests, and a matched-quality continuation benchmark
  define v1.

- [`ot_sinkhorn_unbalanced_ti()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn_unbalanced_ti.md)
  now provides translation-invariant KL-UOT on dense costs and validated
  sparse edge support. The default result is an implicit transport
  operator; explicit edge or dense plans are opt-in. Canonical
  translated potentials, independent fixed-point and KKT residuals,
  generalized-KL primal/dual components, gauge invariance, transported
  mass, exact balanced-support max-flow diagnostics, and runtime/support
  provenance are reported. Analytic, direct-optimization, dense/sparse,
  stiff, unequal- mass, zero-entry, dynamic-range, Hall-deficiency,
  operator, serialization, installed-client, and manifoldalign
  differential evidence guard the path.

- [`as_transport_plan()`](https://bbuchsbaum.github.io/rfugw/reference/as_transport_plan.md)
  and
  [`transport_operator()`](https://bbuchsbaum.github.io/rfugw/reference/transport_operator.md)
  now provide dense, Matrix-CSC sparse, canonical edge-list, and
  implicit coupling contracts. Sparse-safe shape, total/row/column mass,
  forward and adjoint application, barycentric projection, validation,
  linear objective, entropy, and KL avoid hidden dense plan allocation;
  [`transport_plan_materialize()`](https://bbuchsbaum.github.io/rfugw/reference/transport_plan_materialize.md)
  is explicit. Zero transported rows/columns abstain rather than
  becoming uniform matches, and adjoint/reverse barycentric actions are
  documented as non-inverses. Pruning reports lost mass and invalidates
  every stale result certificate. Dense/sparse/operator parity,
  malformed support, representative allocation, saveRDS/readRDS,
  installed-package, and legacy dense-result tests define the boundary.

- [`ot_partial_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_sinkhorn.md)
  now provides certified entropy-regularized fixed-mass partial linear
  OT under the explicit counting-measure entropy-minus-one convention.
  Moderate scaling and genuine log-domain Dykstra paths agree; unsafe
  scaling fails and auto dispatch selects log. Results expose the full
  objective, capacity/mass feasibility, capacity and mass duals,
  stationarity, complementarity, primal-dual gap, trace, and a
  backend-bound warm state. Independent constrained-dual,
  zero/full-mass, tiny-entropy, zero-weight, rectangular,
  constant/duplicate-cost, warm-state, and dispatch tests certify the
  boundary. Entropic partial GW and FGW now reuse this public primitive
  and propagate every required inner certificate.

- [`penalized_partial_fused_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/penalized_partial_fused_gromov_wasserstein.md)
  now optimizes square-loss FGW over the variable-mass subcoupling
  domain with an explicit linear penalty on unmatched source and target
  mass. The penalty is present in the full derivative, certified
  [`ot_partial_penalized()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_penalized.md)
  direction solves, exact line-search polynomial, and reported unrooted
  objective. Results expose transported/discarded mass, every objective
  term, trace feasibility, nested solver status, and a Frank-Wolfe
  stationarity gap. Independent continuous, derivative, line-search,
  fixed-mass-limit, permutation, scaling, asymmetric, and downstream
  fail-closed fixtures distinguish it from fixed-mass partial FGW and
  generalized-KL FUGW.

- [`ot_partial_penalized()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_penalized.md)
  now provides certified exact variable-mass partial linear OT for
  finite, potentially unequal measures. Its `discard_penalty` has
  explicit cost-per-discarded-unit semantics on both marginals, so
  larger values weakly favor more transport. A symmetric dummy reduction
  uses the declared penalty directly—no hidden big-M—and returns
  transported/discarded masses, objective components, and scaled
  augmented primal/dual certificates. Independent original-variable LP
  and direction-mutation tests guard against the historical
  reversed-penalty and accidentally fixed-mass defects.

- [`ot_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn.md)
  now exposes a regularized primal and dual under the explicit
  convention `<M,T> + epsilon * KL(T || p %o% q)`, their gap
  certificate, and the exact constant offset from the entropy-minus-one
  convention, while preserving legacy `ot_dist` as the linear plan cost.
  The new
  [`ot_sinkhorn_divergence()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn_divergence.md)
  debiases one cross and two self objectives, retains every component
  result/status/residual/runtime, supports reusable component dual
  states, and fails closed on component or nonnegativity failure.

- [`ot_wasserstein_cost()`](https://bbuchsbaum.github.io/rfugw/reference/ot_wasserstein_cost.md)
  and
  [`ot_wasserstein_distance()`](https://bbuchsbaum.github.io/rfugw/reference/ot_wasserstein_distance.md)
  now make the ground power and root convention explicit. They accept
  raw supports with a declared metric or a nonnegative cost matrix with
  an explicit existing power, retain scale and probability-normalization
  provenance, distinguish exact values from entropic-plan estimates, and
  fail closed on uncertified solver output. Analytic and independent
  one-dimensional p=1/p=2 oracles, permutation, symmetry, scaling,
  duplicates, zero weights, and rectangular fixtures guard the semantic
  boundary.

- KL structural loss remains deliberately unsupported after an explicit
  scope evaluation. An independent four-index prototype matches the
  factorized formula on positive inputs, while ordinary zero-diagonal
  distance costs make the unfloored directed-KL objective infinite and
  common logarithm floors materially change it. Public GW-family APIs
  now name this decision in their rejection;
  `inst/kl-structural-loss-decision.md` records the evidence and the
  gate for reconsideration.

- The misleading experimental name
  [`lowrank_gromov_wasserstein_samples()`](https://bbuchsbaum.github.io/rfugw/reference/lowrank_gromov_wasserstein_samples.md)
  is deprecated in favor of
  [`dense_gromov_wasserstein_plan_svd()`](https://bbuchsbaum.github.io/rfugw/reference/dense_gromov_wasserstein_plan_svd.md).
  The replacement makes its lifecycle explicit: it materializes two
  dense structure costs and a dense GW plan before truncated SVD,
  reports those memory semantics in its result, and makes no end-to-end
  low-rank scaling claim. Legacy POT-shaped factorized-cost, Dykstra,
  seed, and warning parameters were never operational and are now
  rejected when supplied.

- [`ot_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn.md)
  now accepts reusable `init_duals` and `init_plan` warm starts and
  returns canonical source/target potentials plus serializable
  `dual_state`. Scaling and log backends share the
  source-weighted-mean-zero gauge, duals take documented precedence over
  a simultaneously supplied plan, and every warm solve must satisfy a
  fresh marginal certificate. A deterministic continuation baseline and
  a public-API-only manifoldalign OT-Procrustes fixture guard the
  downstream contract.

- [`ot_sinkhorn_unbalanced()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn_unbalanced.md)
  now has explicit finite-measure semantics through
  `normalization = "none"`, `"joint"`, or `"separate"`. Results retain
  original and effective measures/masses, transported mass,
  generalized-KL objective terms, and a fixed-point certificate. The
  backward-compatible 0.1 default is still separate probability
  normalization; absolute-mass callers should opt into `"none"`
  explicitly.

### Compatibility

- Public solver names, aliases, and diagnostic field meanings are
  specified in `inst/solver-contract.md`.
- Unregularized GW/FGW names mean conditional gradient with an exact
  linear-OT subproblem. They do not claim a global minimizer of the
  non-convex GW/FGW objective.

### Solver honesty

- An executable numerical path matrix now maps every advertised
  symmetry, regularization, Sinkhorn, precision, warm-start, adapter,
  threading, and approximation path to its comparison law and CI scope.
- A bounded mutation-proof harness recreates and kills the reviewed
  transpose, false-simplex-success, tolerance-floor,
  ignored-inner-status, missing-KL-mass, and zero-support-KL defect
  classes without modifying shipping source.
- Distinct PR, nightly, and release numerical-trust workflows now
  publish replayable ledgers. The release matrix checks and installs the
  exact tarball on Linux, macOS, and Windows, with a separate ASan/UBSan
  gate; evidence files keep hosted success separate from publication
  status.
- Native batch kernels now copy R inputs into owned C++ storage before
  OpenMP, use only read-only shared caches and per-job writes, capture
  worker exceptions for main-thread reporting, suppress nested OpenMP,
  and report requested, effective, and maximum thread counts. A pure-C++
  ThreadSanitizer harness and 1/2/4-thread differential tests guard the
  contract.
- Canonical square-loss GW algebra and batch-thread adapters now have
  dedicated translation units; transport certificates and approximation
  caches have R-independent headers, with float/double cache invariants
  implemented once as templates. `inst/native-architecture.md` records
  the remaining solver-engine boundary and the measured
  compile/object-size cost of the split.
- Benchmark evidence now distinguishes certified comparisons from
  fixed-budget performance sentinels. `max_iter` and `experimental` rows
  can never become convergence claims; quality, setup/solve/end-to-end
  time, allocations, precision, threads, commit, and environment are
  separate fields/artifacts. Sampled GW records exact-reference error at
  three budgets, and threshold changes require a digest-matched retained
  evidence-history entry.
- Public convergence is now certificate-based across balanced, partial,
  and semirelaxed, and unbalanced families: returned plans must be
  finite and nonnegative, formulation feasibility must pass, required
  inner solves must be certified, and reported objectives and named
  components must agree with independent recomputation. Failures report
  `infeasible`, `objective_mismatch`, `inner_failure`, or
  `numerical_failure`; iteration count alone cannot imply success.
- Precision-selecting solvers report requested/effective tolerances,
  requested/effective/actual compute precision, backend transitions, and
  a termination reason. Tight mixed requests are promoted to strict
  double instead of silently flooring `1e-9` requests to `1e-6`;
  explicit `strict_double` blocks the large-problem float acceleration
  path. Balanced entropic FGW now defaults to an explicit dynamic-range
  `auto` Sinkhorn policy, selecting scaling only within the documented
  exponent regime and otherwise selecting and reporting the robust
  log-domain backend. Linear KL-unbalanced OT has a genuine log-domain
  implementation; the former UCOOT `sinkhorn_log` scaling alias now
  errors as deprecated. Entropic partial GW/FGW explicitly reject unsafe
  scaling and log requests until a genuine log-domain Dykstra backend
  exists.
- Nested exact and entropic solvers now report final and maximum inner
  residuals, total inner iterations, and `inner_converged`. Exact
  partial directions preserve simplex certificates instead of extracting
  only the plan, and an uncertified required inner solve can no longer
  imply outer convergence. Deterministic nested fault injection covers
  iteration-limit, numerical, infeasible, and non-finite failures.
- Asymmetric exact and entropic partial GW/FGW were first quarantined,
  then re-enabled after the target-marginal transpose defect was fixed
  in a shared square-loss tensor module. O(n^4) gradient,
  directional-derivative, permutation, exact-LP, and entropic-projection
  oracles now guard the path.
- UCOOT and across-spaces defaults now use a positive `epsilon` and
  succeed on documented minimal examples. Unsupported `l2` / `mm` /
  `lbfgsb` choices are rejected before computation.
- `G0` warm starts are consumed by the C++ exact-CG backend, not only by
  lpSolve.
- Flagship solvers auto-detect cost symmetry. `symmetric = TRUE` is
  validated rather than treated as a silent fast-path override.
- Results expose `status`, `converged`, `residual`, and marginal
  residuals. Iteration-limit and numerical-failure exits are distinct
  from convergence.
- Alpha, regularization, iteration limits, and non-finite inputs are
  validated consistently on the flagship solvers.
- Fractional iteration/rank/budget counts are rejected instead of
  truncated; partial solvers record whether transported mass was
  explicit or defaulted; symmetry detection now uses a scale-aware
  absolute-plus-relative tolerance.
- [`ot_kl()`](https://bbuchsbaum.github.io/rfugw/reference/ot_kl.md) now
  implements generalized KL, including unequal-mass correction,
  zero-reference support rules, and finite nonnegative reference
  validation.

### Linear optimal transport

- Exact transport now returns explicit termination reasons, dual
  potentials, primal and dual objectives, marginal residuals, minimum
  nonbasic reduced cost, duality gap, and the tolerances used to certify
  them. Only a feasible, dual-feasible, zero-gap result is reported as
  optimal; deterministic fault injection covers every internal
  termination branch.
- Public primitives
  [`ot_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn.md),
  [`ot_emd()`](https://bbuchsbaum.github.io/rfugw/reference/ot_emd.md),
  and
  [`ot_sinkhorn_unbalanced()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn_unbalanced.md)
  wrap the C++ Sinkhorn and simplex backends and return `rfugw_result`
  objects. They are certified against analytic identities and POT
  0.9.6.post1 fixtures (`inst/extdata/fixtures/linear_ot_fixture.json`).
- [`ot_partial_emd()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_emd.md)
  now exposes exact nonnegative-cost partial linear OT by reducing it to
  the existing certified simplex with one dummy source/target. It
  reports partial mass/feasibility plus the augmented primal/dual
  certificate and is guarded by analytic, independent-LP, metamorphic,
  adversarial, regression, and local baseline evidence. Entropic partial
  OT, Sinkhorn divergence, and fixed-support Wasserstein weight
  barycenters remain explicitly deferred with admission criteria in
  `inst/linear-ot-foundations.md`.
- Independent helpers
  [`ot_validate_plan()`](https://bbuchsbaum.github.io/rfugw/reference/ot_validate_plan.md),
  [`ot_linear_cost()`](https://bbuchsbaum.github.io/rfugw/reference/ot_linear_cost.md),
  [`ot_entropy()`](https://bbuchsbaum.github.io/rfugw/reference/ot_entropy.md),
  [`ot_kl()`](https://bbuchsbaum.github.io/rfugw/reference/ot_kl.md),
  [`ot_gw_square()`](https://bbuchsbaum.github.io/rfugw/reference/ot_gw_square.md),
  [`ot_fgw_square()`](https://bbuchsbaum.github.io/rfugw/reference/ot_fgw_square.md),
  and
  [`ot_barycentric_project()`](https://bbuchsbaum.github.io/rfugw/reference/ot_barycentric_project.md)
  check plans and recompute objectives.

### Results

- Flagship solvers return class `rfugw_result` with
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), and accessors
  [`rfugw_plan()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_plan.md),
  [`rfugw_value()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_value.md),
  [`rfugw_status()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_status.md),
  and
  [`rfugw_residuals()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_residuals.md).
  Legacy list fields remain.

### Documentation

- `_pkgdown.yml` groups the reference by family and labels sampled /
  low-rank APIs as experimental. Flagship solvers have `\examples{}`.
  The solver guide compares formulation, regularization, maturity, and
  performance.
  [`vignette("linear-ot")`](https://bbuchsbaum.github.io/rfugw/articles/linear-ot.md)
  is the task-oriented linear-OT guide.

### Build and optional dependencies

- Default compilation is portable (`-O2`, no `-march=native` /
  `-ffast-math`). Aggressive flags are opt-in via `RFUGW_FAST_FLAGS`.
  Windows `Makevars.win` is provided. OpenMP remains optional. Sanitizer
  and `-Werror` flags are opt-in via `RFUGW_EXTRA_CXXFLAGS` /
  `RFUGW_EXTRA_LIBS`.
- Hosted CI runs multi-OS `R CMD check`, serial/ASan smokes,
  generated-file consistency, and a pkgdown rebuild. The local
  one-command gate is `Rscript tools/release-gate.R`.
- Speed evidence uses `inst/bench/PROTOCOL.md`. The entry point
  `Rscript inst/bench/run_protocol.R` records environment metadata and
  rejects invalid-quality runs.
- Unregularized partial solvers document and test the `lpSolve` Suggests
  route. Entropic partial and default `cpp_transport` exact CG do not
  need it. Exact and entropic partial GW/FGW now run their outer loops
  in C++. The inner entropic projection keeps a finite plan at high
  transported mass.
- UCOOT and across-spaces KL Sinkhorn now run their BCD loop in C++ with
  reusable sample/feature workspaces and inner warm-start telemetry.
- Sampled GW and dense-plan SVD compression stay experimental. The
  certified envelope is tiny-versus-full budget quality, rank
  reconstruction, and the explicitly separated input/solve-memory
  evidence only; unusable budgets and ranks warn or error instead of
  silently clamping. See `inst/bench/sampled-budget-curves.md`.
- Threading is certified only on the batched multiset kernels: 1-thread
  and N-thread plans match, BLAS is pinned to avoid oversubscription,
  and `structure_knn` is documented as a metric change rather than a
  memory saving. See `inst/bench/threading-memory.md`.
- Hosted flagship gates run the canonical protocol with quality checks
  first. The PR gate covers FGW, FUGW, and semirelaxed. Nightly adds
  partial, UCOOT, sampled, larger sizes, and a 1-vs-2 thread smoke. Time
  caps are median `solve_ms` slack, not tight speed claims. Exit 1 is a
  solver regression; exit 2 is infrastructure.

### Metadata

- Package authorship, license holder, repository URL, bug reports, and
  citation metadata now identify the real maintainer.
- `CONTRIBUTING.md` records the compatibility, test, benchmark,
  generated documentation, and evidence policy for 0.1.
