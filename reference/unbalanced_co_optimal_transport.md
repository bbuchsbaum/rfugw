# Unbalanced Co-Optimal Transport

POT-compatible UCOOT wrapper (`reg_type = "independent"`).

## Usage

``` r
unbalanced_co_optimal_transport(
  X,
  Y,
  wx_samp = NULL,
  wx_feat = NULL,
  wy_samp = NULL,
  wy_feat = NULL,
  reg_marginals = 10,
  epsilon = 0.01,
  divergence = c("kl"),
  unbalanced_solver = c("sinkhorn", "sinkhorn_log"),
  alpha = 0,
  M_samp = NULL,
  M_feat = NULL,
  rescale_plan = TRUE,
  init_pi = NULL,
  init_duals = NULL,
  max_iter = 100L,
  tol = 1e-07,
  max_iter_ot = 500L,
  tol_ot = 1e-07,
  log = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- X:

  Source matrix (`n_sample_x x n_feature_x`).

- Y:

  Target matrix (`n_sample_y x n_feature_y`).

- wx_samp:

  Source sample weights.

- wx_feat:

  Source feature weights.

- wy_samp:

  Target sample weights.

- wy_feat:

  Target feature weights.

- reg_marginals:

  Marginal relaxation(s), length 1 or 2.

- epsilon:

  Entropic regularization(s), scalar or length 2. Must be positive;
  default `1e-2`.

- divergence:

  Only `"kl"` is supported. `"l2"` is rejected.

- unbalanced_solver:

  `"sinkhorn"` is the supported scaling-domain implementation. The
  former `"sinkhorn_log"` scaling alias is deprecated and errors because
  it was not a genuine log-domain solver. `"mm"` and `"lbfgsb"` are also
  rejected.

- alpha:

  Linear-term coefficient(s), scalar or length 2.

- M_samp:

  Optional sample linear cost matrix.

- M_feat:

  Optional feature linear cost matrix.

- rescale_plan:

  Rescale sample/feature plans to equal mass each BCD step.

- init_pi:

  Optional list with `pi_samp` and `pi_feat` initial couplings.

- init_duals:

  Accepted for POT-shaped signatures and ignored.

- max_iter:

  Max BCD iterations.

- tol:

  BCD stopping tolerance on sample coupling change.

- max_iter_ot:

  Max iterations for inner unbalanced Sinkhorn solves.

- tol_ot:

  Inner unbalanced Sinkhorn tolerance.

- log:

  If `TRUE`, return diagnostics and objective decomposition.

- verbose:

  If `TRUE`, print BCD diagnostics.

- ...:

  Additional arguments. Unused extras are rejected when the solver uses
  `.reject_unused_dots()`; otherwise they are forwarded to the primary
  solver.

## Value

A list with sample and feature couplings; with `log = TRUE`, includes
objective diagnostics.
