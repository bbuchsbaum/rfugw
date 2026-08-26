# Fit Multiset Alignment to a Fixed Template

Aligns a list of attributed metric-measure sets to a fixed template
using FGW (`method = "fgw_entropic"`) or unbalanced FUGW
(`method = "fugw_kl"`). This is the first stage of general multiset
alignment and returns subject-to-template couplings.

## Usage

``` r
multialign_fit(
  subjects,
  template_mode = c("fixed", "learned"),
  template = NULL,
  k_template = NULL,
  method = c("fgw_entropic", "fugw_kl"),
  alpha = 0.5,
  epsilon = 0.05,
  feature_normalization = c("none", "zscore"),
  feature_cost_metric = c("euclidean", "sqeuclidean"),
  structure_normalize = TRUE,
  structure_knn = 0L,
  structure_knn_min_n = 600L,
  feature_cost_normalize = TRUE,
  solver = c("PGD", "PPA"),
  precision = c("mixed", "double", "strict_double"),
  sinkhorn_method = c("auto", "scaling", "log"),
  max_iter = 200L,
  tol = 1e-09,
  sinkhorn_max_iter = 500L,
  sinkhorn_tol = 1e-09,
  reg_marginals = c(10, 10),
  max_iter_ot = 500L,
  tol_ot = 1e-07,
  rescale_plan = TRUE,
  check_every = 10L,
  coarse_init = c("none", "sampled_graph"),
  coarse_init_components = 10L,
  coarse_init_diffusion_time = 1,
  coarse_init_self_loop = 1e-06,
  coarse_init_nb_samples = c(8L, 2L),
  coarse_init_epsilon = 0.1,
  coarse_init_max_iter = 40L,
  coarse_init_use_cpp = TRUE,
  coarse_init_sampling = c("stochastic", "deterministic"),
  structure_rank = 0L,
  autotune = FALSE,
  autotune_level = c("moderate", "aggressive"),
  seed = 1L,
  use_cpp_batch = TRUE,
  use_cpp_feature_fused = TRUE,
  use_warm_start = TRUE,
  n_threads = 1L,
  batch_min_subjects = 2L,
  template_max_iter = 8L,
  template_tol = 1e-06,
  template_relax = 1,
  barycenter_eps = 1e-12,
  final_refit = FALSE
)
```

## Arguments

- subjects:

  A non-empty list of subject sets. Each subject is a list with `C`
  (required), optional `F`, optional `w`, optional `id`.

- template_mode:

  Template strategy: `"fixed"` (fit couplings to a fixed template) or
  `"learned"` (alternating align-and-update barycenter-style template
  refinement).

- template:

  Fixed template set with the same fields as a subject list. If `NULL`,
  a template is built from pooled subject features via
  [`multialign_make_template()`](https://bbuchsbaum.github.io/rfugw/reference/multialign_make_template.md).

- k_template:

  Number of template nodes used only when `template = NULL`.

- method:

  Alignment backend: `"fgw_entropic"` or `"fugw_kl"`.

- alpha:

  Feature/structure trade-off parameter.

- epsilon:

  Entropic regularization for the selected backend.

- feature_normalization:

  Feature normalization mode before cost assembly.

- feature_cost_metric:

  Feature cost metric used to form cross-set `M`.

- structure_normalize:

  If `TRUE`, each structure matrix is scaled to have max value 1 before
  solving.

- structure_knn:

  Optional kNN sparsification level for structure matrices. If positive,
  rows keep `k+1` smallest entries (including diagonal) and other
  entries are set to a large fill value before solving.

- structure_knn_min_n:

  Minimum matrix size for applying `structure_knn`.

- feature_cost_normalize:

  If `TRUE`, each cross-feature cost matrix is scaled to have max value
  1.

- solver:

  `fgw_entropic` solver (`"PGD"` or `"PPA"`).

- precision:

  Numeric precision for `fgw_entropic` (defaults to `"mixed"` for speed
  in multiset alignment).

- sinkhorn_method:

  Sinkhorn variant for `fgw_entropic`.

- max_iter:

  Outer iterations (`fgw_entropic`) or BCD iterations (`fugw_kl`).

- tol:

  Outer stopping tolerance for selected backend.

- sinkhorn_max_iter:

  Inner Sinkhorn iterations for `fgw_entropic`.

- sinkhorn_tol:

  Inner Sinkhorn tolerance for `fgw_entropic`.

- reg_marginals:

  Marginal relaxation for `fugw_kl`.

- max_iter_ot:

  Inner Sinkhorn iterations for `fugw_kl`.

- tol_ot:

  Inner Sinkhorn tolerance for `fugw_kl`.

- rescale_plan:

  Plan rescaling option for `fugw_kl`.

- check_every:

  Iteration check interval for selected backend.

- coarse_init:

  Optional coarse coupling initialization mode. `"none"` disables coarse
  initialization. `"sampled_graph"` computes a fast sampled-GW plan on
  graph-derived diffusion coordinates to initialize FGW.

- coarse_init_components:

  Diffusion components for `"sampled_graph"`.

- coarse_init_diffusion_time:

  Diffusion-time exponent for coarse init.

- coarse_init_self_loop:

  Diagonal stabilizer for coarse graph init.

- coarse_init_nb_samples:

  Sample counts for coarse sampled-GW gradients.

- coarse_init_epsilon:

  Entropic regularization for coarse init.

- coarse_init_max_iter:

  Iterations for coarse sampled-GW init.

- coarse_init_use_cpp:

  Use C++ fast path for coarse sampled-GW init.

- coarse_init_sampling:

  Sampling mode for coarse init.

- structure_rank:

  Optional low-rank rank used in the FGW tensor-product update. Use
  `0`/`"off"` to disable or `"auto"` for size-aware rank selection.

- autotune:

  If `TRUE`, apply size-aware caps to FGW iterations for faster
  large-`n` runs.

- autotune_level:

  Autotune aggressiveness.

- seed:

  Optional seed used only when `template = NULL`.

- use_cpp_batch:

  If `TRUE`, use the C++ batched FGW path for fixed-template
  `method = "fgw_entropic"`.

- use_cpp_feature_fused:

  If `TRUE` and conditions allow, use a fused C++ path that computes
  feature costs and solves FGW in one batched call without materializing
  `M_list` in R.

- use_warm_start:

  If `TRUE`, reuse previous couplings as initialization across
  learned-template outer iterations (and for `fugw_kl` fixed solves).

- n_threads:

  Requested subject-parallel threads for the C++ batched path.

- batch_min_subjects:

  Minimum number of subjects before enabling batch path.

- template_max_iter:

  Maximum alternating updates for `template_mode = "learned"`.

- template_tol:

  Relative stopping tolerance for learned-template updates.

- template_relax:

  Relaxation in `(0, 1]` when updating learned templates.

- barycenter_eps:

  Stabilizer used in learned-template denominator updates.

- final_refit:

  If `TRUE`, run one final alignment against the last updated learned
  template. If `FALSE` (default), skip this extra pass for speed.

## Value

A list with `couplings`, `objectives`, `template`, `diagnostics`, and
run settings. Per-subject `diagnostics` include outer and nested-solver
status and residual fields when the selected backend exposes them; batch
paths mark unavailable nested certificates explicitly.
