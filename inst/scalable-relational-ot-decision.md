# Scalable relational OT admission decision

## Decision

No relational-OT approximation is promoted to supported status in this
release. `fugw_factorized()` and `fugw_multiscale()` now provide an experimental
full-support Moment-FUGW path whose implementation avoids dense structure,
feature-cost, and coupling matrices. It passes the frozen candidate
exact-factor, operator-action, four-size scaling, multiscale-efficacy,
metamorphic/adversarial robustness, partial-overlap image, and planted cortical held-out-map gates. It remains experimental because those measurements are
from a dirty source-tree candidate, not a reviewed clean same-commit installed
artifact with hosted cross-platform evidence.

The coordinate-native sampled-GW path also remains deferred: it avoids dense
square *structure inputs* while still materializing the dense rectangular
coupling and multiple equally sized work matrices. The graph wrapper inherits
that boundary. Both sampled paths remain experimental.

This is a measured deferral, not a rejection of sampled GW as a research
method. It prevents an input-memory improvement from being presented as an
end-to-end scalable solver.

## Intended workload and admission bar

The original target is balanced, square-loss relational transport between
5,000 to 10,000 objects per domain, initially for coordinate or sparse-graph
geometry. Moment-FUGW extends the intended workload to fused, KL-unbalanced
alignment and an initial 10,000-to-32,000-node surface target. At those sizes a
supported path must:

- accept a provenance-carrying block, sparse, or factorized structure-cost
  representation and never densify it silently;
- avoid every dense `ns x ns`, `nt x nt`, and `ns x nt` allocation during
  solve and normal operator application;
- return an implicit or sparse `transport_operator()` representation, with
  dense materialization available only through an explicit request;
- show empirical incremental peak-memory exponent below 1.5 and a fixed-work
  full-support runtime exponent in 1.8--2.2 on four sizes reaching 5,000 by
  5,000; the primary runtime fit uses the largest three sizes and retains the
  all-size fit as an overhead diagnostic;
- on deterministic dense-oracle cases up to size 128, meet marginal residual
  `<= 1e-8`, relative square-loss GW objective error `<= 1%`, and normalized
  plan-action error `<= 5%` for held-out signals in at least 95% of declared
  fixtures and seeds;
- at larger sizes, independently certify feasibility, recompute the declared
  approximate objective without materializing the full coupling, and fail the
  gate when any answer-quality threshold is missed.

These thresholds define an admission experiment. They are not claims about
the current experimental functions.

## Candidate audit

### Moment-FUGW: implemented, matrix-free, and still experimental

`factorized_cost()` represents a cost as source-additive, target-additive, and
bilinear terms. Exact squared-Euclidean, cosine, and correlation constructors
never form their full matrices. `ot_sinkhorn_unbalanced_ti()` evaluates these
costs with native blocked log-sum-exp reductions and retains the full-support
coupling through translated potentials. `fugw_factorized()` closes each FUGW
block through marginal Gram matrices and the cross moment `H = Bx' P By`, and
recomputes its objective using implicit plan actions.

The test contract compares factorized TI-UOT to the existing certified
dense solver and compares `fugw_factorized()` to dense `fugw_kl()` on exact
fixtures. It checks both FUGW couplings, every objective component, transported
mass, moments, apply/adjoint actions, zero-weight nodes, asymmetric factors,
permutation equivariance, serialization, and the explicit dense-materialize
boundary. An independent pure-R four-index oracle covers the FUGW algebra. A
large shape test replaces the cost callback with a sentinel that would fail on
a full-block request and verifies tile-sized rather than coupling-sized stored
state. A separate robustness contract now exercises both couplings under
independent permutations, coordinate isometries, factor-gauge changes,
zero-factor padding, split-weight duplicates, cost composition, block-size
changes, serialization, and thread-setting changes. Small epsilon, extreme
rho-to-epsilon ratios, constant correlation profiles, extreme coordinate
scales, and forced underflow must either certify honestly or fail closed. The
local replay currently retains 12 cases and 150 successful expectations across
15 named properties with no failures, warnings, errors, or skips.

The separate core-certification contract turns AC1--AC5 into a release
artifact rather than relying on ordinary test execution. It regenerates 12
dense/oracle comparisons through 128 by 127 (six exact-zero and three positive
near-zero weight cases) from a frozen case-design CSV, the ten-seed multiscale
curve, and 27 named oracle,
post-rescale KKT, geometry, hierarchy, and transfer cases. The local replay has
281 successful expectations with no failures, errors, warnings, or skips.

`fugw_multiscale()` accepts arbitrary-depth coarse-to-fine `fugw_domain()`
lists. It validates parent coverage and weight aggregation, prolongs both pairs
of TI potentials, and lifts each coarse implicit coupling by within-parent
conditional weights. The lift preserves transported mass and recomputes the
fine marginal, entropy, cross-moment, and Gram state. A separate explicit tiny
matrix oracle checks those identities. Per-level epsilon schedules and traces
are retained.

The current support is complete and implicit. It removes coupling-sized
storage but retains full-support `O(ns * nt * r_eff)` rank-dependent arithmetic
per affine-bilinear sweep; it is matrix-free, not subquadratic. The amended
checked-tarball installed profile at 625, 1,250, 2,500, and 5,000 has
incremental peak-RSS exponent 1.324 and primary largest-three runtime exponent
1.943. The all-point runtime exponent, 1.895, is retained as the required
overhead-sensitive diagnostic rather than substituted for the primary
asymptotic fit. All 60 rows and 244 raw receipts are retained, including
separate cumulative `Rprofmem` curves for solve, objective, apply, and adjoint;
native RSS remains the peak-memory measure. The largest completed end-to-end
certified endpoint in that candidate artifact is 512.

The ten-seed planted hierarchy reduced fine-level inner work in all cases with
a 34.6 percent median reduction, objective excess below numerical noise, and
held-out score noninferiority. The frozen image evaluation retained all 40
requested rows and found lower false occluded mass for unbalanced FUGW at each
partial-overlap level. Five fsaverage6-derived cortical fixtures reached
held-out correlations of 0.99942--0.99949; the maps are planted synthetic
geodesic radial-basis channels, not observed neurobiology. A 16-landmark dense
parity case has maximum component/action error 2.27e-11. Geometry-embedding
error is reported separately from solver error.

`support = "adaptive"` still fails closed because no omitted-kernel-mass or
reduced-cost certificate has landed. Domain-specific image/mesh constructors
and automatic adaptive support remain future scope. The implementation makes
no 100k-node or global-optimality claim. Promotion still requires a clean
same-commit installed full artifact, hosted cross-platform replay, and
independent review of the retained evidence. The release workflow now produces
a deterministic technical candidate receipt only after recomputing both RSS
slopes, all frozen scientific gates, the three-platform core and robustness
replays, three successful package checks, and exact artifact provenance; that
workflow has not yet run for this dirty candidate.

### Coordinate-native sampled GW: deferred

`sampled_gromov_wasserstein_coords()` stores coordinate inputs in `O(nd)`, and
the 400/800-node benchmark shows close objective estimates to the dense-input
implementation at the same sampled budget. That comparison is useful
implementation parity, but the reference is itself experimental and is not a
certified dense-GW accuracy oracle.

The native kernel declares full `ns x nt` matrices `T`, `Lik`, and `Lik_eff`,
then creates `new_T` and `diff` while the first three remain live. Thus the
source-audited lower bound is five double matrices, or `40 * ns * nt` bytes,
before allocator, Sinkhorn, conversion, and returned-object overhead. For
equal domains the lower bound is 1.0 GB at `n = 5,000` and 4.0 GB at
`n = 10,000`; the returned plan alone is 200 MB and 800 MB respectively.

The retained one-thread, three-run fixture used budget `(16, 2)`, 120 sampled
iterations, seed 123, and sizes 400 and 800. Median coordinate-solver time rose
from 2.314 s to 8.851 s, a doubling slope of 1.936. R-visible allocations rose
from 3.68 MB to 10.27 MB, but those values exclude native Armadillo allocation
and therefore are not peak-memory evidence. The exact source-audited workspace
lower bound grows with slope 2.0. See
`inst/bench/scalable-relational-admission.csv`.

The existing `n = 16` quality curve establishes only that a full `(ns, nt)`
budget is closer to dense entropic GW than `(2, 1)`. Intermediate budgets were
not monotone, and a full budget does not solve the practical complexity
problem. Logged dense-input sampled results explicitly report
`converged = FALSE` and `experimental_no_convergence_certificate`; the
coordinate path does not add an independent feasibility or objective
certificate.

### Sparse graph wrapper: deferred with its dense-coupling dependency

`sampled_gw_from_graphs()` preserves a sparse graph long enough to construct
diffusion coordinates, then calls the coordinate-native solver. It therefore
inherits the dense `ns x nt` plan and workspaces. Sparse graph provenance does
not make the relational solve or returned coupling sparse, and the current
wrapper cannot satisfy the transport-operator admission requirement.

### Coupling operators and sparse plans

The supported `transport_operator()` and sparse-plan contracts avoid dense
materialization. Moment-FUGW now produces a native implicit coupling rather
than wrapping an already dense plan. Sampled GW still allocates before the
operator boundary, so wrapping its result would not make that path admissible.
Automatic sparse support discovery remains unimplemented.

### Existing low-rank paths: deferred or explicitly ineligible

`dense_gromov_wasserstein_plan_svd()` is post-hoc compression: it constructs
dense structure costs and a dense plan before SVD. It is ineligible by design.
The `structure_rank` optimization used by dense FGW reduces tensor-product
arithmetic but still accepts dense structures and maintains the full
coupling. Neither is a genuine low-rank relational solver with rank-bounded
solve memory.

### Block/operator structure costs: factorized candidate only

The public cost protocol now validates generic block costs and affine-bilinear
factors. The Moment-FUGW solver requires algebraic factors and fails closed for
a callback-only operator structure cost because its exact moment closure is not
available. Generic block reduction remains useful infrastructure, but it is
not advertised as a relational solver path.

## What remains supported and experimental

Dense certified GW/FGW solvers remain the supported reference for problems
that fit their declared memory boundary. Moment-FUGW is an experimental
matrix-free joint-KL FUGW implementation with full implicit support. Sampled
dense, coordinate, and graph functions remain experimental square-loss
approximations. Dense-plan SVD remains experimental post-processing. None is
described as globally optimal or as a supported substitute for dense GW/FUGW.

Reconsider promotion only when a new solver—not a wrapper around a dense
plan—passes the workload, complexity, quality, provenance, and operator gates
above.
