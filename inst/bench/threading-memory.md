# Threading, determinism, and memory

Evidence for `bd-01M05QY8W4132GQQX6VYRPBRCR`. The certified envelope is
in `inst/solver-contract.md`. This note records the measurements.

Regenerate the memory table with:

```sh
RFUGW_RLIB=/path/to/.tmp-lib Rscript inst/bench/profile_threading_memory.R \
  inst/bench/results/threading_memory.csv 20260816 96 4
```

## What is parallel

Only the C++ batched kernels used by `multialign_fit()` (`n_threads`)
and `cpp_feature_cost_batch()` are parallel. They farm subjects, not
inner Sinkhorn or FGW iterations.

Single-subject flagship solvers are serial. Changing `OMP_NUM_THREADS`
does not change their plans. Tests: `tests/testthat/test-threading.R`
and `tests/testthat/test-linear-ot-cert.R`.

`n_threads = 1` is the serial fallback. On a 3-subject FGW batch it
matches `n_threads = 2` to `1e-10` in the objective and couplings.

## Oversubscription

When `n_threads > 1`, `multialign_fit()` pins OpenBLAS / MKL / vecLib /
BLIS to one thread unless `RFUGW_PIN_BLAS_THREADS=0`. Thread-scaling
speed claims apply only to the batch kernels with BLAS pinned. See
`inst/bench/benchmark_thread_scaling.R`. Mixing subject OpenMP with a
threaded BLAS is not a certified speed path.

## RNG and deterministic sampling

Flagship solvers do not draw random numbers. Sampled GW is stochastic;
`random_state` reproduces a run. `sampling = "deterministic"` on the
coordinate path uses top-k indices and is bit-stable across repeated
C++ calls. C++ versus R deterministic plans can differ at mixed
precision; that parity is checked on the GW value in
`tests/testthat/test-sparse-graph-sampled-gw.R`, not as a bit-exact
claim.

## Sanitizers

Hosted ASan/UBSan (`.github/workflows/sanitizer.yml`) builds with
OpenMP flags cleared. They certify serial memory safety. They do not
certify data races. There is no TSan job. Threaded correctness is the
1-versus-N equivalence test.

## Peak memory and densification

Bytes are `object.size()` of inputs, not allocator peak. Conservative
`-O2` install, one thread, seed `20260816`, 4 subjects plus a template.

| n | structure C | kNN-filled C | features | R `M_list` | plans |
|---|---|---|---|---|---|
| 48 | 125320 | 125320 | 6840 | 100256 | 74592 |
| 96 | 432520 | 432520 | 12600 | 346016 | 295776 |

- `structure_knn` writes a still-dense `n x n` matrix, filling far
  entries with `max(C)`. Storage does not drop.
- `use_cpp_feature_fused = TRUE` avoids the R-side `M_list`
  (`n x n` per subject). At `n = 96` that is about 346 KB of extra R
  matrices.
- Returned plans stay dense and already rival the structure storage.

Times at these sizes are at timer resolution and are not a speed claim.
