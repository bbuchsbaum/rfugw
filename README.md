# rfugw

[![Sparse Sampled Perf Gate (Fast)](https://github.com/bbuchsbaum/rfugw/actions/workflows/sparse-sampled-perf-gate.yml/badge.svg)](https://github.com/bbuchsbaum/rfugw/actions/workflows/sparse-sampled-perf-gate.yml)
[![Sparse Sampled Perf Gate (Nightly)](https://github.com/bbuchsbaum/rfugw/actions/workflows/sparse-sampled-perf-gate-nightly.yml/badge.svg)](https://github.com/bbuchsbaum/rfugw/actions/workflows/sparse-sampled-perf-gate-nightly.yml)
[![Flagship protocol gate (PR)](https://github.com/bbuchsbaum/rfugw/actions/workflows/flagship-gate.yml/badge.svg)](https://github.com/bbuchsbaum/rfugw/actions/workflows/flagship-gate.yml)
[![Flagship protocol gate (Nightly)](https://github.com/bbuchsbaum/rfugw/actions/workflows/flagship-gate-nightly.yml/badge.svg)](https://github.com/bbuchsbaum/rfugw/actions/workflows/flagship-gate-nightly.yml)

Fast R + C++ implementations of:

- Balanced and KL-unbalanced linear optimal transport (`ot_sinkhorn`,
  `ot_emd`, `ot_sinkhorn_unbalanced`)
- Entropic Fused Gromov-Wasserstein (FGW, square loss)
- Unregularized FGW via conditional gradient + LP direction
- Fused Unbalanced Gromov-Wasserstein (FUGW, KL divergence with Sinkhorn inner solver)

The implementation follows the POT reference formulations in:

- `ot.gromov.entropic_fused_gromov_wasserstein`
- `ot.gromov.fused_unbalanced_gromov_wasserstein`

## Install (local)

```r
install.packages("rfugw", repos = NULL, type = "source")
```

or from the repo root:

```r
devtools::install("rfugw")
```

Default builds are conservative and portable (R's `-O2`, no `-march=native`,
no `-ffast-math`). OpenMP uses `$(SHLIB_OPENMP_CXXFLAGS)` when the
toolchain provides it. On macOS Homebrew, set `RFUGW_OPENMP_FLAGS` /
`RFUGW_OPENMP_LIBS` if those defaults are not enough; see
`inst/build-profiles.md`.

To opt in to aggressive, machine-local flags (not for CRAN, sanitizers, or
reproducible comparisons):

```bash
RFUGW_FAST_FLAGS="-march=native -ffast-math -fno-math-errno -fno-trapping-math" \
  R CMD INSTALL .
```

The default exact FGW/GW and unregularized partial GW/FGW paths
(`lp_solver = "cpp_transport"`) do not need `lpSolve`. The
`"lp_transport"` and `"lp_matrix"` backends still do. At `n >= 24`,
`lp_solver = "lp_matrix"` is typically faster when `lpSolve` is installed:

```r
install.packages("lpSolve")
```

Entropic partial solvers also do not need `lpSolve`.

## Example

```r
library(rfugw)

set.seed(1)
ns <- 20
nt <- 25
C1 <- as.matrix(dist(matrix(rnorm(ns * 2), ns, 2)))
C2 <- as.matrix(dist(matrix(rnorm(nt * 2), nt, 2)))
M  <- as.matrix(dist(rbind(matrix(rnorm(ns * 2), ns, 2), matrix(rnorm(nt * 2), nt, 2))))[1:ns, (ns + 1):(ns + nt)]

C1 <- C1 / max(C1)
C2 <- C2 / max(C2)
M  <- M / max(M)

out <- fgw_entropic(M, C1, C2, alpha = 0.5, epsilon = 0.05)
str(out)

out_exact <- fgw_exact_cg(M, C1, C2, alpha = 0.5, max_iter = 120)
str(out_exact)
```

## Guides

- `vignette("linear-ot")` — balanced and unbalanced linear OT
- `vignette("solver-guide")` — which GW/FGW solver to call
- `vignette("barycenters")` — fixed-support barycenters
- `vignette("multiset-alignment")` — many subjects to one template

## Multiset Alignment (Fixed Template)

`rfugw` now includes a general multiset interface that aligns multiple
attributed metric-measure sets to a shared template support:

```r
subjects <- list(
  list(C = C1, F = matrix(rnorm(ns * 3), ns, 3), id = "s1"),
  list(C = C2, F = matrix(rnorm(nt * 3), nt, 3), id = "s2")
)

fit <- multialign_fit(
  subjects = subjects,
  template_mode = "fixed",
  template = NULL,      # auto-build fixed template from pooled subject features
  k_template = 30,
  method = "fgw_entropic",
  precision = "mixed",  # default in multiset mode (speed-first)
  alpha = 0.5,
  epsilon = 0.05
)

# Subject -> template coupling
P1 <- fit$couplings$s1
```

Performance knobs for large multiset runs:

- `precision = "mixed"`: fastest default in `multialign_fit()`
- `autotune = TRUE`: size-aware caps for `max_iter`, `sinkhorn_max_iter`, `check_every`
- `structure_rank = r` or `"auto"`: low-rank FGW tensor-product approximation (`0`/`"off"` disables)
- `structure_knn = k`: optional kNN structure sparsification (approximate) for large `n`

## Current Scope

- `square_loss` only
- No Python runtime dependency
- Differential tests against frozen POT fixtures shipped in `inst/extdata/fixtures`
- Supported names, aliases, parameters, and diagnostic meanings are specified
  in [`inst/solver-contract.md`](inst/solver-contract.md)
- Unregularized GW/FGW solvers are conditional-gradient procedures with an
  exact linear-OT subproblem. They do not claim a global minimizer.

## Citation

```r
citation("rfugw")
```

See `NEWS.md` and `CONTRIBUTING.md` for the 0.1 compatibility and evidence
policy.

## Benchmarking

Quality-controlled speed evidence uses `inst/bench/PROTOCOL.md`. The
hosted PR gate is FGW, FUGW, and semirelaxed; nightly also covers
partial, UCOOT, sampled, larger sizes, and 1-vs-2 threads.

```bash
inst/bench/run_flagship_gate.sh inst/bench/results/current \
  fgw,fugw,semirelaxed pr 1 20260816 1
```

Run rfugw benchmarks:

```bash
RFUGW_RLIB=/tmp/Rlib_test Rscript inst/bench/benchmark_suite.R 3 inst/bench/results/benchmark_latest.csv 42
```

By default this runs a strict pre-benchmark accuracy gate and aborts on drift.
To skip the gate (not recommended):

```bash
RFUGW_SKIP_ACCURACY_GATE=1 RFUGW_RLIB=/tmp/Rlib_test Rscript inst/bench/benchmark_suite.R 3 inst/bench/results/benchmark_latest.csv 42
```

Run multiset alignment benchmarks:

```bash
RFUGW_RLIB=/tmp/Rlib_test OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
Rscript inst/bench/benchmark_multiset.R 2 inst/bench/results/benchmark_multiset_latest.csv 42
```

Run large-scale multiset profiling (`S=10`, `n=500/1000/2000`):

```bash
RFUGW_RLIB=/tmp/Rlib_test RFUGW_SKIP_ACCURACY_GATE=1 \
Rscript inst/bench/benchmark_multiset_large.R 1 inst/bench/results/benchmark_multiset_large_latest.csv 42 "1 2 4 8"
```

Profile mixed-path throughput and isolate C++ batch kernel cost:

```bash
RFUGW_RLIB=/tmp/Rlib_test OMP_NUM_THREADS=2 \
Rscript inst/bench/profile_mixed_path.R 2 inst/bench/results/profile_mixed_path_latest.csv 20260222 8 600 2
```

Run the nightly guard bundle (accuracy + large-scale + POT parity checks):

```bash
inst/bench/run_nightly_guard.sh inst/bench/results/nightly 20260222 2
```

You can limit sizes to profile faster, e.g. `n = 500,1000`:

```bash
RFUGW_RLIB=/tmp/Rlib_test RFUGW_SKIP_ACCURACY_GATE=1 \
Rscript inst/bench/benchmark_multiset_large.R 1 inst/bench/results/benchmark_multiset_large_latest.csv 42 "1 4" "500 1000"
```

Run POT reference benchmarks:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python inst/bench/benchmark_pot_reference.py 3 inst/bench/results/pot_benchmark_latest.csv 42
```

Run sparse sampled POT-vs-rfugw performance gate (`n >= 400` by default):

```bash
RFUGW_PYTHON=python3 inst/bench/run_sparse_sampled_perf_gate.sh \
  inst/bench/results/benchmark_sparse_sampled_latest.csv \
  123 \
  1 \
  "400" \
  60 \
  1e-6 \
  400 \
  1.0 \
  1.0 \
  0.005
```

This runs:

- `inst/bench/benchmark_sparse_sampled.R` (with tuned sampled Sinkhorn controls)
- `inst/bench/compare_sparse_sampled_benchmark.R` (crossover markdown report)
- `inst/bench/gate_sparse_sampled_perf.R` (strict ratio/objective gate)

CI automation for this path lives in:

- Fast PR/push gate: `.github/workflows/sparse-sampled-perf-gate.yml`
- Nightly stability gate: `.github/workflows/sparse-sampled-perf-gate-nightly.yml`

Nightly-scale local run example:

```bash
RFUGW_PYTHON=python3 inst/bench/run_sparse_sampled_perf_gate.sh \
  inst/bench/results/benchmark_sparse_sampled_nightly_latest.csv \
  123 \
  3 \
  "400,800" \
  60 \
  1e-6 \
  400 \
  1.0 \
  1.0 \
  0.005
```

Run thread scaling:

```bash
RFUGW_RLIB=/tmp/Rlib_test inst/bench/run_thread_scaling.sh inst/bench/results/thread_scaling_latest.csv 2 42 "1 2 4 8"
```

Generate a markdown report:

```bash
python inst/bench/make_benchmark_report.py \
  inst/bench/results/benchmark_latest.csv \
  inst/bench/results/pot_benchmark_latest.csv \
  inst/bench/results/thread_scaling_latest.csv \
  inst/bench/results/benchmark_report.md
```

## Next Performance Steps

- Adaptive rank selection and sparse-kernel crossover (choose dense vs low-rank automatically).
- Precompute/cache template-side low-rank factors across batched subject solves.
- Add optional block-sparse GW tensor kernels for very large, sparse structure matrices.
