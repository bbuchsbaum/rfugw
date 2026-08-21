# Runtime controls and provenance

Every certified `rfugw_result` contains `runtime_provenance`. The record keeps
the raw environment value, parsed value, default, source, validity, and impact
for each environment control. It also records native effective values (including
the OpenMP thread count and initialized matvec thresholds) and requested-to-
effective automatic transitions for precision, stopping tolerances, Sinkhorn
method, approximation rank, and sampled budgets when those fields apply.

The behavior-changing runtime environment is:

| Control | Effect | Default |
|---|---|---:|
| `RFUGW_ENABLE_WARM_START` | FUGW inner warm starts | `0` |
| `RFUGW_ADAPTIVE_INNER_TOL` | staged FUGW inner tolerances | `0` |
| `RFUGW_ADAPTIVE_INNER_TOL_STAGE1` | early-stage tolerance | `1e-6` |
| `RFUGW_ADAPTIVE_INNER_TOL_STAGE2` | middle-stage tolerance | `5e-7` |
| `RFUGW_SINKHORN_TARGET_TOL_ABS_MULT` | absolute target multiplier | `100` |
| `RFUGW_SINKHORN_TARGET_TOL_MULT` | target multiplier | `52` |
| `RFUGW_SINKHORN_TARGET_TOL_CAP_D` | double target cap | `1.5e-4` |
| `RFUGW_SINKHORN_TARGET_TOL_CAP_F` | float target cap | `2e-4` |
| `RFUGW_SINKHORN_REL_TOL_MULT` | relative update multiplier | `200` |
| `RFUGW_SINKHORN_REL_TOL_FLOOR_D` | double relative floor | `4e-5` |
| `RFUGW_SINKHORN_REL_TOL_FLOOR_F` | float relative floor | `2e-4` |
| `RFUGW_SINKHORN_COL_REL_TOL_MULT` | column-relative multiplier | `1000` |
| `RFUGW_SINKHORN_COL_REL_TOL_FLOOR_D` | double column floor | `1.5e-4` |
| `RFUGW_SINKHORN_COL_REL_TOL_FLOOR_F` | float column floor | `2e-4` |
| `RFUGW_MATVEC_BLOCKED_MIN_WORK_D` | double blocked-matvec threshold | `25000` |
| `RFUGW_MATVEC_GEMV_MIN_WORK_D` | double BLAS threshold | `120000` |
| `RFUGW_MATVEC_BLOCKED_MIN_WORK_F` | float blocked-matvec threshold | `160000` |
| `RFUGW_MATVEC_GEMV_MIN_WORK_F` | float BLAS threshold | `25000` |
| `RFUGW_AUTOTUNE_MATVEC` | session matvec autotuning | `1` |
| `RFUGW_SAMPLED_MIXED` | sampled-solver mixed precision | `0` |
| `RFUGW_SEMIRELAXED_MIXED` | semirelaxed mixed precision | `1` |
| `RFUGW_PIN_BLAS_THREADS` | batch BLAS thread pinning | `1` |
| `OMP_NUM_THREADS` | OpenMP thread request | runtime default |
| `OPENBLAS_NUM_THREADS` | OpenBLAS thread request | library default |
| `MKL_NUM_THREADS` | MKL thread request | library default |
| `VECLIB_MAXIMUM_THREADS` | Accelerate thread request | library default |
| `BLIS_NUM_THREADS` | BLIS thread request | library default |

Boolean controls accept only `0` or `1`; numeric controls accept finite positive
values, except matvec work thresholds which may be zero. Invalid values use the
documented fallback and are named in `runtime_provenance$warnings`. Matvec
thresholds are initialized once per native library session; the provenance says
whether initialization has occurred and reports the actual initialized values.

Build-only controls (`RFUGW_FAST_FLAGS`, `RFUGW_EXTRA_CXXFLAGS`,
`RFUGW_EXTRA_LIBS`, `RFUGW_OPENMP_FLAGS`, and `RFUGW_OPENMP_LIBS`) do not alter
an already installed library. The release gate records them separately alongside
the source identity and artifact digest. See `inst/thread-safety.md` for the
native ownership boundary, nested-parallelism policy, and sanitizer evidence.
