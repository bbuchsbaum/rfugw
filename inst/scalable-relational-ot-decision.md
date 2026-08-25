# Scalable relational OT admission decision

## Decision

No relational-OT approximation is promoted to supported status in this
release. The coordinate-native sampled-GW path is the closest candidate, but
it fails the defining complexity requirement: it avoids dense square
*structure inputs* while still materializing the dense rectangular coupling
and multiple equally sized work matrices. The graph wrapper inherits that
boundary. Both remain experimental.

This is a measured deferral, not a rejection of sampled GW as a research
method. It prevents an input-memory improvement from being presented as an
end-to-end scalable solver.

## Intended workload and admission bar

The target is balanced, square-loss relational transport between 5,000 to
10,000 objects per domain, initially for coordinate or sparse-graph geometry.
At those sizes a supported path must:

- accept a provenance-carrying block, sparse, or factorized structure-cost
  representation and never densify it silently;
- avoid every dense `ns x ns`, `nt x nt`, and `ns x nt` allocation during
  solve and normal operator application;
- return an implicit or sparse `transport_operator()` representation, with
  dense materialization available only through an explicit request;
- show empirical peak-memory slope below 1.5 under domain-size doubling and a
  runtime slope below 2.0 on at least four sizes reaching 5,000 by 5,000;
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

### Coordinate-native sampled GW: closest, but deferred

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

### Coupling operators and sparse plans: valuable, but downstream of solve

The supported `transport_operator()` and sparse-plan contracts avoid dense
materialization when a solver actually produces an implicit or sparse
coupling. They do not currently provide block structure costs or a relational
optimization algorithm. Wrapping an already dense sampled plan would move the
allocation boundary after the prohibited allocation, so it is not an
admissible scalable path.

### Existing low-rank paths: deferred or explicitly ineligible

`dense_gromov_wasserstein_plan_svd()` is post-hoc compression: it constructs
dense structure costs and a dense plan before SVD. It is ineligible by design.
The `structure_rank` optimization used by dense FGW reduces tensor-product
arithmetic but still accepts dense structures and maintains the full
coupling. Neither is a genuine low-rank relational solver with rank-bounded
solve memory.

### Block/operator structure costs: no candidate implementation

No public GW/FGW solver accepts a validated block or operator structure cost.
Promoting a name without an implementation, provenance checks, independent
objectives, or failure semantics would create a parallel uncaught contract.
This remains the most direct future design direction, potentially paired with
a factorized or sparse coupling representation.

## What remains supported and experimental

Dense certified GW/FGW solvers remain the supported reference for problems
that fit their declared memory boundary. Sampled dense, coordinate, and graph
functions remain experimental square-loss approximations. Dense-plan SVD
remains experimental post-processing. None is described as globally optimal,
end-to-end low-memory, or a supported substitute for dense GW.

Reconsider promotion only when a new solver—not a wrapper around a dense
plan—passes the workload, complexity, quality, provenance, and operator gates
above.
