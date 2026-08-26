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
| Balanced regularized OT and Sinkhorn divergence | product-reference KL primal/dual; cross/source-self/target-self solves | one declared entropy/reference convention; symmetric zero-diagonal self costs or raw metric supports | independent primal recomputation and dual gap; explicit constant offset; all three divergence components certified; symmetry, debiasing, nonnegativity tolerance, and joint cost-epsilon scale law |
| Wasserstein summaries | `ot_wasserstein_cost`, `ot_wasserstein_distance`; exact or Sinkhorn-plan path | raw supports with declared metric, or nonnegative cost with explicit existing power; positive-mass weights | independent 1-D p=1/p=2 oracle; explicit p-cost/root convention; exact versus entropic-plan label; fail-closed underlying certificate; scale and permutation laws |
| Fixed-support Wasserstein barycenter weights | `ot_barycenter_weights` exact joint LP or semi-debiased Sinkhorn mode | nonnegative cross costs with common target column count; probability measures; normalized nonnegative coefficients; symmetric common-support self cost for regularized mode | independent joint LP or convex optimizer; certified components; exact LP primal/dual/reduced costs or analytic dual gradient; simplex feasibility and KKT; monotone safeguarded descent; identity/permutation/zero/duplicate/tiny-entropy/warm-state laws |
| KL-unbalanced linear OT | `ot_sinkhorn_unbalanced` | finite cost/weights; positive epsilon and relaxation | finite nonnegative plan; recomputed full objective; fixed-point residual; certified inner termination |
| Translation-invariant sparse KL-UOT | `ot_sinkhorn_unbalanced_ti`; dense, edge/CSR sparse, implicit operator, explicit edge plan | finite measures with positive totals; positive epsilon/rho; validated duplicate-free support; every active node covered | dense/complete-sparse parity; direct finite-measure optimizer; translated-potential KKT and independently recomputed fixed point; generalized-KL primal-dual gap; gauge invariance; exact max-flow support diagnostic; coupling-free apply and installed downstream parity |
| Entropic fixed-mass partial OT | `ot_partial_sinkhorn` scaling/log Dykstra and auto dispatch | finite cost; probability-normalized weights; mass in `[0,1]`; positive epsilon; backend-matched certified warm state | row/column capacities and exact mass; explicit counting-measure entropy-minus-one objective; independent constrained-dual oracle; capacity/mass dual KKT and primal-dual gap; safe dynamic-range dispatch; scaling/log and partial-GW/FGW parity |
| Transport plan/operator representations | dense matrix, Matrix CSC sparse, canonical edge list, implicit forward/adjoint operator | declared positive shape; finite nonnegative stored weights or masses; duplicate policy; dimension-valid callbacks | dense/sparse parity for mass, actions, barycentric projection, linear objective, entropy/KL, and validation; explicit zero-mass abstention; no hidden dense allocation; pruning loss and certificate invalidation; serialization/install parity |
| Penalized variable-mass partial OT | `ot_partial_penalized` exact symmetric-dummy reduction | finite nonnegative cost and finite measures; nonnegative discard penalty | row/column capacities; optimized transported/discarded mass; independent original-variable LP; objective decomposition; scaled augmented primal/dual/reduced-cost certificate; monotone penalty direction |
| Penalized variable-mass partial FGW | certified Frank-Wolfe with `ot_partial_penalized` directions | finite feature/structure costs and finite measures; nonnegative discard penalty; feasible free-mass start | every trace iterate satisfies capacities; finite-difference gradient; direct line-search polynomial; independent tiny optimizer; full objective decomposition; nested direction certificates; Frank-Wolfe stationarity gap; lambda, fixed-mass-limit, permutation, and compatible-scale laws |
| Balanced GW/FGW | `fgw_entropic` PGD/PPA, scaling/log, mixed/double/strict-double, symmetric/general/low-rank; `fgw_exact_cg` native/lpSolve | square finite structures; feasible weights/start; declared square loss | finite nonnegative balanced plan; independently recomputed objective; final inner certificate; truthful local-CG termination |
| Partial GW/FGW | exact native/lpSolve and entropic native; GW/FGW aliases | finite structures; feasible upper-bound weights/start; declared transported mass | nonnegative plan; marginal upper bounds and transported mass; independently recomputed objective; every final LP/projection certified |
| Semirelaxed GW/FGW | exact/entropic, symmetric/general, R/native | finite structures; fixed source weights | nonnegative plan; source marginal; learned target mass; independently recomputed objective and stopping rule |
| FUGW | `fugw_kl` double/mixed, warm/cold inner starts | finite structures/cost/weights; positive epsilon and relaxation | both couplings finite/nonnegative; independently recomputed components; final and maximum inner residuals; all required inner solves certified |
| UCOOT/across spaces | joint/independent KL paths | finite data/weights; supported solver and positive epsilon | both couplings finite/nonnegative; recomputed components; final/max inner residuals; truthful BCD termination |
| GW/FGW barycenters and multialign | structure/fused support-learning barycenters; serial/threaded batches | every constituent solver contract plus compatible shapes | constituent certificates, finite template, deterministic serial/thread equivalence, no hidden failed subject; never conflated with fixed-support linear-Wasserstein weights |
| Approximate/experimental | sampled dense/coordinate/graph and `dense_gromov_wasserstein_plan_svd`; no scalable path admitted | explicit budget/rank and experimental status | finite feasible output; quality versus exact baseline; declared approximation; monotone envelope only where tested; explicit dense-plan materialization and solve-memory order; measured admission deferral with target-size workspace lower bounds; performance evidence separate from correctness |

Structural GW/FGW support is deliberately square-loss only. Directed KL
structural loss is deferred because zero-diagonal distance costs make the exact
objective infinite and an implicit logarithm floor changes the estimand; see
`inst/kl-structural-loss-decision.md` for the independent prototype and gate.

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
