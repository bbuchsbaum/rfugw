# rfugw numerical trust charter

This charter is the executable-evidence policy for numerical changes. It is
stricter than the pre-audit suite that passed while asymmetric partial
gradients, tolerance reporting, nested status, convergence certification, and
generalized KL semantics were wrong. A passing example or objective snapshot is
not sufficient evidence for a solver contract.

## Contracts and supported paths

Every supported solver must declare its input domain, returned formulation,
stopping rule, feasibility certificate, nested-solver obligations, numerical
precision, and implementation backend. The matrix below is the minimum test
surface; aliases must share the primary implementation and diagnostics.

| Family | Public paths and backends | Preconditions | Required postconditions and certificates |
|---|---|---|---|
| Balanced linear OT | `ot_sinkhorn` scaling/log; `ot_emd` native simplex | finite cost; finite nonnegative positive-mass weights | finite nonnegative plan; both marginals; recomputed linear objective; Sinkhorn final marginal residual or simplex primal feasibility, reduced cost, and duality gap |
| KL-unbalanced linear OT | `ot_sinkhorn_unbalanced` | finite cost/weights; positive epsilon and relaxation | finite nonnegative plan; recomputed full objective; fixed-point residual; certified inner termination |
| Balanced GW/FGW | `fgw_entropic` PGD/PPA, scaling/log, mixed/double/strict-double, symmetric/general/low-rank; `fgw_exact_cg` native/lpSolve | square finite structures; feasible weights/start; declared square loss | finite nonnegative balanced plan; independently recomputed objective; final inner certificate; truthful local-CG termination |
| Partial GW/FGW | exact native/lpSolve and entropic native; GW/FGW aliases | finite structures; feasible upper-bound weights/start; declared transported mass | nonnegative plan; marginal upper bounds and transported mass; independently recomputed objective; every final LP/projection certified |
| Semirelaxed GW/FGW | exact/entropic, symmetric/general, R/native | finite structures; fixed source weights | nonnegative plan; source marginal; learned target mass; independently recomputed objective and stopping rule |
| FUGW | `fugw_kl` double/mixed, warm/cold inner starts | finite structures/cost/weights; positive epsilon and relaxation | both couplings finite/nonnegative; independently recomputed components; final and maximum inner residuals; all required inner solves certified |
| UCOOT/across spaces | joint/independent KL paths | finite data/weights; supported solver and positive epsilon | both couplings finite/nonnegative; recomputed components; final/max inner residuals; truthful BCD termination |
| Barycenters and multialign | structure/fused barycenters; serial/threaded batches | every constituent solver contract plus compatible shapes | constituent certificates, finite template, deterministic serial/thread equivalence, no hidden failed subject |
| Approximate/experimental | sampled dense/coordinate/graph and post-hoc low-rank | explicit budget/rank and experimental status | finite feasible output; quality versus exact baseline; declared approximation; monotone envelope only where tested; memory/performance evidence separate from correctness |

## Failure modes and mandatory test families

| Failure mode | Required test family | Independent evidence |
|---|---|---|
| Wrong objective, tensor, or gradient algebra | contract, differential, metamorphic, regression | O(n^4) enumeration, finite differences, directional derivatives, symmetric/general equivalence |
| Wrong feasible set, mass, or support | contract, property/fuzz, adversarial | analytic marginal/mass laws, zero support/weights, singletons and rectangular plans |
| False convergence or hidden nested failure | contract, fault injection, regression | `converged => certificates`; starved inner budgets; injected max-iter, nonfinite, infeasible, and branch failures |
| Exact-OT false optimality | differential, adversarial, native branch tests | analytic cases, independent LP, primal/dual objectives, reduced costs, gap and degenerate bases |
| Precision, tolerance, or dispatch drift | boundary, metamorphic, regression | requested/effective fields, immediately-above/below thresholds, double/mixed/strict-double comparisons |
| Scale, shift, and conditioning instability | metamorphic, adversarial, property | additive-shift and scaling laws, tiny epsilon, large dynamic range, near-degenerate inputs |
| Nondeterminism or unreplayable fuzz failure | property/fuzz | fixed RNG kind and seed, serialized minimal counterexample, one-command replay |
| Approximation quality regression | differential, metamorphic, performance | exact baseline, budget/rank curve, correctness threshold before timing threshold |
| Thread/backend divergence | differential, concurrency, sanitizer | 1-versus-N worker comparison, native thread harness, ASan/UBSan/TSan where supported |
| Performance erosion | benchmark after quality gate | fixed problem generator, seed, environment, solve/e2e/memory measures, quality-valid rows only |

## Tolerances

Floating-point comparisons use `abs(x - y) <= atol + rtol * scale`, where
`scale = max(abs(x), abs(y), problem_scale)`. Each test records `atol`, `rtol`,
precision, conditioning or scale rationale, and the expected error-growth path.
Float/mixed tests may use wider justified tolerances than strict double tests,
but a tolerance must remain at least ten times tighter than the defect it is
designed to catch. Exact equality is reserved for discrete status, shape,
seed, dispatch, and serialized-provenance contracts.

Random tests use an explicit seed and RNG kind. A failure prints the seed,
solver arguments, backend, precision, environment controls, and a replay
command. Nightly fuzzers save the smallest failing fixture as RDS plus a text
manifest; no user data or absolute machine paths belong in the fixture.

## CI scopes

- PR: fast contract tests, known regressions, small analytic/oracle cases,
  gradient laws, critical metamorphic properties, and mutation sentinels.
- Nightly: seeded randomized differentials, adversarial precision/scale grids,
  optional backends, threading, sanitizer jobs, approximation curves, and
  answer-quality-controlled performance trends.
- Release: clean source tarball and installed package; Linux, macOS, and
  Windows; required hosted checks; exact release commit; full oracle/fuzz
  corpus; sanitizers; generated artifacts; benchmark evidence only after all
  numerical gates; artifact digest and replayable evidence ledger.

The executable workflows are `.github/workflows/numerical-trust.yml`
(PR), `numerical-trust-nightly.yml`, and `numerical-trust-release.yml`.
They invoke `tools/numerical-trust/run-laws.R`; the release workflow tests the
installed tarball on Linux, macOS, and Windows and runs a separate ASan/UBSan
job. `collect-evidence.R` records the exact commit, artifact digest, seed and
replay contract, selected environment controls, representative certificates,
toolchain session, and evidence channel. Hosted workflow completion is hosted
evidence only; publication remains explicitly `not_evaluated`.

Local, hosted, benchmark, and publication evidence are distinct. A focused
test pass proves only that focused scope. A performance row is inadmissible if
its solver result lacks the required numerical certificate.

`Rscript tools/numerical-trust/run-mutation-proof.R` is the bounded mutation
sentinel. It recreates the six reviewed defect classes without editing shipping
source, requires each targeted invariant to reject its mutant, records runtime
and assertion evidence, and demonstrates why scalar-objective-only and
symmetric-only suites cannot detect the asymmetric transpose defect.

## Admission rule

A new numerical feature or backend is not supported until it has an explicit
contract, a metamorphic or property test, an adversarial edge case, an
independent oracle where feasible, mutation-sensitive regression evidence, and
a representative performance baseline. Unsupported paths must fail clearly or
remain explicitly experimental.
