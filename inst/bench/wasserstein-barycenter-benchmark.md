# Fixed-support Wasserstein barycenter benchmark

Run `Rscript inst/bench/benchmark_wasserstein_barycenter.R OUTPUT.csv` against
an installed rfugw package. The deterministic fixture holds a five-point
support, three probability measures, barycentric coefficients, costs, and
answer-quality requirements fixed. It compares the exact joint linear program,
a cold regularized solve, and continuation from the cold solve's certified
state.

The retained run used rfugw 0.0.1.9000 from the working tree based on
`4b701b6`, R 4.5.1, Apple arm64, and macOS 14.3. The machine-readable
measurements are in `wasserstein-barycenter-baseline.csv`.

| mode | objective | KKT residual / tolerance | outer iterations | component solves | component iterations | warm reuses | elapsed (s) | allocated bytes |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| exact joint LP | 0.0125000 | 0 / 1e-8 | 1 | 3 | 0 | 0 | 0.010 | 1,633,576 |
| regularized cold | 0.0517715 | 1.97e-5 / 2e-5 | 83 | 336 | 6,175 | 332 | 0.312 | 17,080,904 |
| regularized warm | 0.0517715 | 1.97e-5 / 2e-5 | 0 | 4 | 0 | 4 | 0.004 | 228,960 |

Every row converged, met its declared KKT tolerance, and satisfied the simplex
constraint. The warm solve reproduced the cold weights exactly on this
fixture, required no outer or component iterations, and reused all four
eligible component states. Timings and recorded allocations are descriptive
local observations, not cross-platform performance guarantees. The certified
claim is parity at matched answer quality with no greater iterative work.
