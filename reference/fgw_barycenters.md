# Fixed-Support FGW Barycenters

Learns a fixed-size barycenter support by alternating between (entropic)
FGW couplings and closed-form barycenter updates of structure and
features.

## Usage

``` r
fgw_barycenters(
  N,
  Ys,
  Cs,
  ps = NULL,
  lambdas = NULL,
  alpha = 0.5,
  fixed_structure = FALSE,
  fixed_features = FALSE,
  p = NULL,
  loss_fun = "square_loss",
  epsilon = 0.05,
  symmetric = TRUE,
  max_iter = 100L,
  tol = 1e-09,
  warmstartT = TRUE,
  solver = c("PGD", "PPA"),
  sinkhorn_max_iter = 500L,
  sinkhorn_tol = 1e-09,
  sinkhorn_method = c("auto", "scaling", "log"),
  precision = c("mixed", "double", "strict_double"),
  check_every = 10L,
  init_C = NULL,
  init_X = NULL,
  random_state = 1L,
  feature_cost_metric = c("euclidean", "sqeuclidean"),
  feature_cost_normalize = TRUE,
  barycenter_eps = 1e-12,
  verbose = FALSE
)
```

## Arguments

- N:

  Number of barycenter nodes.

- Ys:

  List of feature matrices, one per sample (`ns x d`).

- Cs:

  List of structure matrices, one per sample (`ns x ns`).

- ps:

  Optional list of sample weights for each sample.

- lambdas:

  Optional sample weights across sets (default uniform).

- alpha:

  FGW feature/structure trade-off.

- fixed_structure:

  If `TRUE`, keep barycenter structure fixed at `init_C`.

- fixed_features:

  If `TRUE`, keep barycenter features fixed at `init_X`.

- p:

  Optional barycenter weights (`N`) (default uniform).

- loss_fun:

  Inner loss, currently only `"square_loss"` is supported.

- epsilon:

  Entropic regularization used for inner FGW solves.

- symmetric:

  Assume symmetric structure matrices in inner FGW solves.

- max_iter:

  Maximum outer barycenter iterations.

- tol:

  Stopping tolerance on barycenter updates.

- warmstartT:

  If `TRUE`, warm-start inner FGW with previous couplings.

- solver:

  Inner FGW solver (`"PGD"` or `"PPA"`).

- sinkhorn_max_iter:

  Inner Sinkhorn iterations for FGW solves.

- sinkhorn_tol:

  Inner Sinkhorn tolerance for FGW solves.

- sinkhorn_method:

  Sinkhorn variant (`"scaling"` or `"log"`).

- precision:

  Numeric precision mode for inner FGW solves.

- check_every:

  Outer FGW check interval for inner solves.

- init_C:

  Optional initial barycenter structure (`N x N`).

- init_X:

  Optional initial barycenter feature matrix (`N x d`).

- random_state:

  Optional seed used only when `init_C` is `NULL`.

- feature_cost_metric:

  Feature cost metric (`"euclidean"` or `"sqeuclidean"`).

- feature_cost_normalize:

  If `TRUE`, normalize feature costs to `[0, 1]`.

- barycenter_eps:

  Stabilizer for divisions in barycenter updates.

- verbose:

  If `TRUE`, print outer-iteration diagnostics.

## Value

A list with `X`, `C`, `p`, `couplings`, `history`, `iterations`,
`objective`, outer `status`/`converged`, and compact final
`solver_diagnostics` for each inner FGW problem. `history` includes
aggregate inner-solver diagnostics for every barycenter iteration.
