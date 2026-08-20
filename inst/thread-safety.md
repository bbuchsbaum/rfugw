# Native batch thread-safety contract

`rfugw` parallelizes across independent batch items. It does not call R, use an
Rcpp proxy, allocate an R object, emit an R warning, or raise an R exception
inside an OpenMP worker. Each batch entry crosses the ownership boundary on the
R main thread: matrices and vectors from R lists are copied into owned
`arma::mat` / `arma::vec` storage before a parallel region begins. Shared
template matrices and low-rank caches are read-only. Every worker writes only
to result, timing, diagnostic, and failure slots indexed by its own job.

Native exceptions are caught by `run_worker_guarded()` inside each worker and
stored in a per-job string. After the region completes, the R main thread turns
the first captured failure into an R error. No C++ exception is permitted to
cross an OpenMP boundary.

Nested native parallelism is suppressed. Inner matrix kernels check
`omp_in_parallel()` before starting OpenMP work, and batch results report
`nested_parallel = FALSE`. Batch outputs also report `requested_threads`,
`used_threads`, and `max_threads`. `runtime_provenance` records the OpenMP and
BLAS environment, while `RFUGW_PIN_BLAS_THREADS=1` (the default) pins supported
BLAS libraries during `multialign_fit()` batches to avoid OpenMP-by-BLAS
oversubscription.

## Executable assurance

- `tests/testthat/test-native-concurrency.R` compares deterministic outputs for
  1, 2, and 4 requested threads, checks isolated exception capture, verifies
  thread provenance, and statically rejects R/Rcpp calls in batch regions.
- `tests/testthat/test-threading.R` compares real FGW batch plans and objectives
  at 1, 2, and 4 requested threads.
- `tools/numerical-trust/thread_kernel_harness.cpp` is a pure-C++ OpenMP
  harness for shared read-only caches, disjoint writes, nested suppression, and
  an injected worker exception. The nightly workflow builds and runs it under
  ThreadSanitizer.

The regular ASan/UBSan package check is deliberately serial: sanitizer evidence
for memory and undefined behavior is kept separate from concurrency evidence.
The ThreadSanitizer harness is the authoritative hosted race check. On builds
without OpenMP, the same tests run the serial fallback and honestly report
`used_threads = max_threads = 1`; they do not claim local threaded execution.
