# Package index

## Linear optimal transport

Public balanced, exact, and KL-unbalanced linear-OT primitives.

- [`ot_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn.md)
  : Balanced entropic optimal transport
- [`ot_sinkhorn_divergence()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn_divergence.md)
  : Debiased Sinkhorn divergence with certified component solves
- [`ot_emd()`](https://bbuchsbaum.github.io/rfugw/reference/ot_emd.md) :
  Exact balanced linear transport
- [`ot_partial_emd()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_emd.md)
  : Exact partial linear optimal transport
- [`ot_partial_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_sinkhorn.md)
  : Entropy-Regularized Fixed-Mass Partial Optimal Transport
- [`ot_partial_penalized()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_penalized.md)
  : Exact penalized variable-mass partial optimal transport
- [`ot_sinkhorn_unbalanced()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn_unbalanced.md)
  : KL-unbalanced entropic optimal transport
- [`ot_sinkhorn_unbalanced_ti()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn_unbalanced_ti.md)
  : Translation-Invariant Sparse KL-UOT
- [`ot_wasserstein_cost()`](https://bbuchsbaum.github.io/rfugw/reference/ot_wasserstein_cost.md)
  : Wasserstein p-cost with explicit cost-power semantics
- [`ot_wasserstein_distance()`](https://bbuchsbaum.github.io/rfugw/reference/ot_wasserstein_distance.md)
  : Wasserstein p-distance with explicit root semantics
- [`ot_validate_plan()`](https://bbuchsbaum.github.io/rfugw/reference/ot_validate_plan.md)
  : Validate a transport plan
- [`ot_linear_cost()`](https://bbuchsbaum.github.io/rfugw/reference/ot_linear_cost.md)
  : Linear transport cost
- [`ot_entropy()`](https://bbuchsbaum.github.io/rfugw/reference/ot_entropy.md)
  : Entropic term of a plan
- [`ot_kl()`](https://bbuchsbaum.github.io/rfugw/reference/ot_kl.md) :
  Generalized KL divergence of a plan from a product reference
- [`ot_gw_square()`](https://bbuchsbaum.github.io/rfugw/reference/ot_gw_square.md)
  : Square-loss Gromov-Wasserstein objective
- [`ot_fgw_square()`](https://bbuchsbaum.github.io/rfugw/reference/ot_fgw_square.md)
  : Square-loss fused Gromov-Wasserstein objective
- [`ot_barycentric_project()`](https://bbuchsbaum.github.io/rfugw/reference/ot_barycentric_project.md)
  : Barycentric projection of a coupling
- [`transport_problem_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/transport_problem_sinkhorn.md)
  : Construct an explicit balanced entropic transport problem
- [`transport_problem_emd()`](https://bbuchsbaum.github.io/rfugw/reference/transport_problem_emd.md)
  : Construct an explicit exact balanced transport problem
- [`transport_problem_partial_emd()`](https://bbuchsbaum.github.io/rfugw/reference/transport_problem_partial_emd.md)
  : Construct an explicit exact fixed-mass partial transport problem
- [`transport_problem_partial_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/transport_problem_partial_sinkhorn.md)
  : Construct an explicit entropic fixed-mass partial transport problem
- [`transport_problem_unbalanced()`](https://bbuchsbaum.github.io/rfugw/reference/transport_problem_unbalanced.md)
  : Construct an explicit KL-unbalanced transport problem
- [`transport_solve()`](https://bbuchsbaum.github.io/rfugw/reference/transport_solve.md)
  : Solve a versioned transport problem
- [`transport_capabilities()`](https://bbuchsbaum.github.io/rfugw/reference/transport_capabilities.md)
  : Query public transport-solver capabilities
- [`print(`*`<rfugw_transport_problem>`*`)`](https://bbuchsbaum.github.io/rfugw/reference/print.rfugw_transport_problem.md)
  : Print a transport-problem contract
- [`print(`*`<rfugw_solver_state>`*`)`](https://bbuchsbaum.github.io/rfugw/reference/print.rfugw_solver_state.md)
  : Print reusable solver state
- [`rfugw_plan()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_plan.md)
  : Extract the coupling from an rfugw result
- [`rfugw_value()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_value.md)
  : Extract the documented objective value
- [`rfugw_status()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_status.md)
  : Extract solver status
- [`rfugw_residuals()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_residuals.md)
  : Extract residual diagnostics
- [`rfugw_state()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_state.md)
  : Extract reusable solver state
- [`rfugw_problem()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_problem.md)
  : Extract the explicit transport problem contract
- [`rfugw_masses()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_masses.md)
  : Extract source, target, and transported-mass provenance
- [`rfugw_provenance()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_provenance.md)
  : Extract stable solver and runtime provenance
- [`print(`*`<rfugw_result>`*`)`](https://bbuchsbaum.github.io/rfugw/reference/print.rfugw_result.md)
  : Print an rfugw result
- [`summary(`*`<rfugw_result>`*`)`](https://bbuchsbaum.github.io/rfugw/reference/summary.rfugw_result.md)
  : Summarize an rfugw result
- [`print(`*`<summary.rfugw_result>`*`)`](https://bbuchsbaum.github.io/rfugw/reference/print.summary.rfugw_result.md)
  : Print an rfugw result summary

## Transport plan representations

Dense, sparse, edge-supported, and implicit coupling operations.

- [`as_transport_plan()`](https://bbuchsbaum.github.io/rfugw/reference/as_transport_plan.md)
  [`dim(`*`<rfugw_transport_plan>`*`)`](https://bbuchsbaum.github.io/rfugw/reference/as_transport_plan.md)
  [`as.matrix(`*`<rfugw_transport_plan>`*`)`](https://bbuchsbaum.github.io/rfugw/reference/as_transport_plan.md)
  [`print(`*`<rfugw_transport_plan>`*`)`](https://bbuchsbaum.github.io/rfugw/reference/as_transport_plan.md)
  : Construct a Transport Plan Representation
- [`transport_operator()`](https://bbuchsbaum.github.io/rfugw/reference/transport_operator.md)
  : Construct an Implicit Transport Operator
- [`transport_plan_shape()`](https://bbuchsbaum.github.io/rfugw/reference/transport_plan_shape.md)
  : Inspect Transport Plan Shape
- [`transport_plan_mass()`](https://bbuchsbaum.github.io/rfugw/reference/transport_plan_mass.md)
  : Inspect Transported Mass
- [`transport_plan_apply()`](https://bbuchsbaum.github.io/rfugw/reference/transport_plan_apply.md)
  : Apply a Transport Plan as a Linear Operator
- [`transport_plan_adjoint()`](https://bbuchsbaum.github.io/rfugw/reference/transport_plan_adjoint.md)
  : Apply the Transport Adjoint
- [`transport_plan_barycentric()`](https://bbuchsbaum.github.io/rfugw/reference/transport_plan_barycentric.md)
  : Barycentric Projection Through a Transport Plan
- [`transport_plan_materialize()`](https://bbuchsbaum.github.io/rfugw/reference/transport_plan_materialize.md)
  : Explicitly Materialize a Transport Plan
- [`transport_plan_representation()`](https://bbuchsbaum.github.io/rfugw/reference/transport_plan_representation.md)
  : Inspect Transport Representation
- [`transport_plan_prune()`](https://bbuchsbaum.github.io/rfugw/reference/transport_plan_prune.md)
  : Prune a Transport Plan and Invalidate Stale Certificates

## Flagship GW and FGW

Primary unregularized and entropic solvers. Unregularized names are
conditional-gradient procedures, not global optimizers.

- [`fgw_entropic()`](https://bbuchsbaum.github.io/rfugw/reference/fgw_entropic.md)
  : Entropic Fused Gromov-Wasserstein (square loss)
- [`fgw_entropic2()`](https://bbuchsbaum.github.io/rfugw/reference/fgw_entropic2.md)
  : Entropic Fused Gromov-Wasserstein Objective Value
- [`fgw_exact_cg()`](https://bbuchsbaum.github.io/rfugw/reference/fgw_exact_cg.md)
  : Unregularized FGW via Conditional Gradient + LP Direction (square
  loss)
- [`fugw_kl()`](https://bbuchsbaum.github.io/rfugw/reference/fugw_kl.md)
  : Fused Unbalanced Gromov-Wasserstein (KL divergence, Sinkhorn inner
  solver)
- [`fugw_kl2()`](https://bbuchsbaum.github.io/rfugw/reference/fugw_kl2.md)
  : Fused Unbalanced Gromov-Wasserstein Objective Value
- [`gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/gromov_wasserstein.md)
  : Unregularized Gromov-Wasserstein via Conditional Gradient (square
  loss)
- [`gromov_wasserstein2()`](https://bbuchsbaum.github.io/rfugw/reference/gromov_wasserstein2.md)
  : Unregularized Gromov-Wasserstein Objective Value
- [`fused_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/fused_gromov_wasserstein.md)
  : POT Alias for Unregularized FGW
- [`fused_gromov_wasserstein2()`](https://bbuchsbaum.github.io/rfugw/reference/fused_gromov_wasserstein2.md)
  : POT Alias for Unregularized FGW Value
- [`entropic_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_gromov_wasserstein.md)
  : Entropic Gromov-Wasserstein (square loss)
- [`entropic_gromov_wasserstein2()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_gromov_wasserstein2.md)
  : Entropic Gromov-Wasserstein Objective Value
- [`entropic_fused_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_fused_gromov_wasserstein.md)
  : POT Alias for Entropic FGW
- [`entropic_fused_gromov_wasserstein2()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_fused_gromov_wasserstein2.md)
  : POT Alias for Entropic FGW Value
- [`fused_unbalanced_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/fused_unbalanced_gromov_wasserstein.md)
  : POT Alias for Fused Unbalanced GW
- [`fused_unbalanced_gromov_wasserstein2()`](https://bbuchsbaum.github.io/rfugw/reference/fused_unbalanced_gromov_wasserstein2.md)
  : POT Alias for Fused Unbalanced GW Value

## Partial, semirelaxed, and unbalanced variants

- [`partial_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/partial_gromov_wasserstein.md)
  : Partial Gromov-Wasserstein (square loss)
- [`partial_gromov_wasserstein2()`](https://bbuchsbaum.github.io/rfugw/reference/partial_gromov_wasserstein2.md)
  : Partial Gromov-Wasserstein Objective Value
- [`partial_fused_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/partial_fused_gromov_wasserstein.md)
  : Partial Fused Gromov-Wasserstein (square loss)
- [`partial_fused_gromov_wasserstein2()`](https://bbuchsbaum.github.io/rfugw/reference/partial_fused_gromov_wasserstein2.md)
  : Partial Fused Gromov-Wasserstein Objective Value
- [`penalized_partial_fused_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/penalized_partial_fused_gromov_wasserstein.md)
  : Penalized Variable-Mass Partial Fused Gromov-Wasserstein
- [`entropic_partial_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_partial_gromov_wasserstein.md)
  : Entropic Partial Gromov-Wasserstein (square loss)
- [`entropic_partial_gromov_wasserstein2()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_partial_gromov_wasserstein2.md)
  : Entropic Partial Gromov-Wasserstein Objective Value
- [`entropic_partial_fused_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_partial_fused_gromov_wasserstein.md)
  : Entropic Partial Fused Gromov-Wasserstein (square loss)
- [`entropic_partial_fused_gromov_wasserstein2()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_partial_fused_gromov_wasserstein2.md)
  : Entropic Partial Fused Gromov-Wasserstein Objective Value
- [`semirelaxed_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/semirelaxed_gromov_wasserstein.md)
  : Semi-Relaxed Gromov-Wasserstein (non-entropic, square loss)
- [`semirelaxed_gromov_wasserstein2()`](https://bbuchsbaum.github.io/rfugw/reference/semirelaxed_gromov_wasserstein2.md)
  : Semi-Relaxed Gromov-Wasserstein Objective Value
- [`semirelaxed_fused_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/semirelaxed_fused_gromov_wasserstein.md)
  : Semi-Relaxed Fused Gromov-Wasserstein (non-entropic, square loss)
- [`semirelaxed_fused_gromov_wasserstein2()`](https://bbuchsbaum.github.io/rfugw/reference/semirelaxed_fused_gromov_wasserstein2.md)
  : Semi-Relaxed Fused Gromov-Wasserstein Objective Value
- [`entropic_semirelaxed_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_semirelaxed_gromov_wasserstein.md)
  : Entropic Semirelaxed Gromov-Wasserstein (square loss)
- [`entropic_semirelaxed_gromov_wasserstein2()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_semirelaxed_gromov_wasserstein2.md)
  : Entropic Semirelaxed Gromov-Wasserstein Objective Value
- [`entropic_semirelaxed_fused_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_semirelaxed_fused_gromov_wasserstein.md)
  : Entropic Semirelaxed Fused Gromov-Wasserstein (square loss)
- [`entropic_semirelaxed_fused_gromov_wasserstein2()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_semirelaxed_fused_gromov_wasserstein2.md)
  : Entropic Semirelaxed Fused Gromov-Wasserstein Objective Value
- [`fused_unbalanced_across_spaces_divergence()`](https://bbuchsbaum.github.io/rfugw/reference/fused_unbalanced_across_spaces_divergence.md)
  : Fused Unbalanced Across-Spaces Divergence (KL Sinkhorn)
- [`unbalanced_co_optimal_transport()`](https://bbuchsbaum.github.io/rfugw/reference/unbalanced_co_optimal_transport.md)
  : Unbalanced Co-Optimal Transport
- [`unbalanced_co_optimal_transport2()`](https://bbuchsbaum.github.io/rfugw/reference/unbalanced_co_optimal_transport2.md)
  : Unbalanced Co-Optimal Transport Objective Value

## Barycenters and multiset alignment

- [`ot_barycenter_weights()`](https://bbuchsbaum.github.io/rfugw/reference/ot_barycenter_weights.md)
  : Fixed-support Wasserstein barycenter weights
- [`print(`*`<rfugw_barycenter_result>`*`)`](https://bbuchsbaum.github.io/rfugw/reference/print.rfugw_barycenter_result.md)
  : Print a fixed-support Wasserstein barycenter result
- [`fgw_barycenters()`](https://bbuchsbaum.github.io/rfugw/reference/fgw_barycenters.md)
  : Fixed-Support FGW Barycenters
- [`gromov_barycenters()`](https://bbuchsbaum.github.io/rfugw/reference/gromov_barycenters.md)
  : Gromov-Wasserstein Barycenters
- [`entropic_gromov_barycenters()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_gromov_barycenters.md)
  : Entropic Gromov-Wasserstein Barycenters
- [`entropic_fused_gromov_barycenters()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_fused_gromov_barycenters.md)
  : Entropic Fused Gromov-Wasserstein Barycenters
- [`fused_gromov_barycenters()`](https://bbuchsbaum.github.io/rfugw/reference/fused_gromov_barycenters.md)
  : Fused Gromov-Wasserstein Barycenters
- [`multialign_fit()`](https://bbuchsbaum.github.io/rfugw/reference/multialign_fit.md)
  : Fit Multiset Alignment to a Fixed Template
- [`multialign_make_template()`](https://bbuchsbaum.github.io/rfugw/reference/multialign_make_template.md)
  : Build a Fixed Template Support for Multiset Alignment
- [`multialign_normalize_plan()`](https://bbuchsbaum.github.io/rfugw/reference/multialign_normalize_plan.md)
  : Normalize a Coupling Matrix
- [`multialign_project_features()`](https://bbuchsbaum.github.io/rfugw/reference/multialign_project_features.md)
  : Project Template Features to Subject Nodes
- [`multialign_project_matrix()`](https://bbuchsbaum.github.io/rfugw/reference/multialign_project_matrix.md)
  : Project a Subject Matrix into Template Space

## Experimental

Approximate APIs. See inst/solver-contract.md. Do not treat these as
certified 0.1 flagship solvers.

- [`sampled_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/sampled_gromov_wasserstein.md)
  : Sampled Gromov-Wasserstein (square loss)
- [`sampled_gromov_wasserstein_coords()`](https://bbuchsbaum.github.io/rfugw/reference/sampled_gromov_wasserstein_coords.md)
  : Sampled Gromov-Wasserstein Directly from Coordinates
- [`sampled_gw_from_graphs()`](https://bbuchsbaum.github.io/rfugw/reference/sampled_gw_from_graphs.md)
  : Sampled GW from Sparse or Dense Similarity Graphs
- [`graph_diffusion_coordinates()`](https://bbuchsbaum.github.io/rfugw/reference/graph_diffusion_coordinates.md)
  : Diffusion Coordinates from a Similarity Graph
- [`dense_gromov_wasserstein_plan_svd()`](https://bbuchsbaum.github.io/rfugw/reference/dense_gromov_wasserstein_plan_svd.md)
  : Dense GW Plan Followed by Truncated SVD (Experimental)

## Deprecated

Compatibility names retained only for migration.

- [`lowrank_gromov_wasserstein_samples()`](https://bbuchsbaum.github.io/rfugw/reference/lowrank_gromov_wasserstein_samples.md)
  : Deprecated pseudo-low-rank GW compatibility wrapper
