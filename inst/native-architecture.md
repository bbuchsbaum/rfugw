# Native architecture and ownership boundaries

The native layer is split by contract, while `fgw_core.cpp` remains the
consolidated translation unit for tightly coupled numerical engines. Public R
API policy and result certification remain in R, where they can be inspected
without reading optimized kernels.

| Boundary | Canonical location | Responsibility |
|---|---|---|
| Linear OT certificate | `src/transport_simplex.h` | R-independent termination enum and primal/dual/feasibility certificate shape. The pivot engine and its R adapter remain in the consolidated core. |
| Square-loss GW algebra | `src/gw_square.h`, `src/gw_square.cpp` | One compiled implementation of forward/reverse tensors, loss, and gradient for symmetric and general costs. Exact and partial solvers call this contract. |
| Approximation caches | `src/approximation_cache.h` | Precision-templated low-rank cache dimensions and square-cache initialization shared by float and double engines. |
| Batching and threading | `src/batch_worker.h`, `src/thread_probe.cpp` | Pure-C++ guarded workers, shared-cache probe, nested suppression, and the thin R adapter. |
| Solver engines | `src/fgw_core.cpp` | Balanced entropic, exact CG, semirelaxed, partial, sampled, FUGW, and UCOOT arithmetic. Exported functions are adapters around internal core results. |
| Diagnostics and provenance | `R/result.R`, `R/solver_helpers.R` | Status, convergence implications, independent objective/feasibility checks, precision transitions, runtime environment, and stable accessors. |
| Solver-family APIs | `R/fgw.R`, `R/linear_ot.R`, `R/pot_api.R`, `R/pot_api_ports.R`, `R/multiset.R` | Validation, formulation-specific policy, dispatch, and user-facing compatibility. |

Optimized entropic paths retain BLAS-factorized symmetric/general tensor
products in `fgw_core.cpp`; this is the one deliberate second implementation.
It is performance-critical and is justified by the brute-force O(n^4),
directional-derivative, permutation, exact-LP, and entropic differential laws in
`test-gw-gradient-laws.R`, `test-partial-asymmetric.R`, and
`test-advertised-path-laws.R`. Exact and partial CG paths do not reconstruct
forward/reverse tensors: they call `rfugw::gw_square_terms()`.

The headers containing algebra, certificates, approximation contracts, and
worker contracts do not depend on Rcpp. Their R adapters live in translation
units or on the R main thread, so the contracts can be compiled and tested
without invoking a complete solver wrapper. `thread_kernel_harness.cpp` does
exactly this for the worker contract; focused R adapters expose the algebra and
certificate contracts to the independent oracle suites.

## Compile and object evidence

Measurements used the same macOS arm64 debug `devtools::load_all()` toolchain
(Homebrew clang 20.1.8, C++17). The baseline was captured immediately before
the module extraction; the post-split rebuild was captured after moving GW
implementations out of headers.

| Measure | Before split | After split | Interpretation |
|---|---:|---:|---|
| `fgw_core.cpp` bytes | 270,408 | 267,752 | Canonical adapters/algebra left the monolith. |
| `fgw_core.o` bytes | 6,989,832 | 6,634,808 | Core debug object decreased 5.1%. |
| All new native debug objects | 6,989,832 | 8,474,840 | Separate debug objects add symbol/debug overhead; this is not presented as a size optimization. |
| Linked `rfugw.so` bytes | 2,668,744 | 2,678,664 | Linked size increased 0.4%. |
| Timed debug rebuild wall time | 21.4 s | 21.92 s | Within 0.52 s on a single local observation; a benchmark claim would require repeated clean builds. |

These numbers are architecture evidence, not a performance certification. The
post-split runtime and numerical behavior is guarded by the full Gate 2 law
suite and the release gate.
