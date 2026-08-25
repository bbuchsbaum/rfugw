# Solver-protocol continuation benchmark

Run `Rscript inst/bench/benchmark_solver_protocol.R OUTPUT.csv` against an
installed rfugw package. The script holds the 40-by-32 rectangular cost,
probability measures, `epsilon = 0.65`, log-domain backend, and requested
answer quality fixed. It compares cold solves with a continuation sequence at
tolerances `1e-5`, `1e-7`, and `1e-9`; only the preceding certified protocol
state changes.

The retained run used rfugw 0.0.1.9000 from the working tree based on
`4b701b6`, R 4.5.1, Apple arm64, and macOS 14.3. The raw measurements are in
`solver-protocol-baseline.csv`.

| tolerance | warm iterations | cold iterations | maximum plan difference | regularized-objective difference |
|---:|---:|---:|---:|---:|
| 1e-5 | 40 | 40 | 0 | 0 |
| 1e-7 | 30 | 70 | 2.43e-17 | 4.44e-16 |
| 1e-9 | 20 | 90 | 2.78e-17 | 4.44e-16 |

The continuation path used 90 iterations in total versus 200 for cold
restarts. Both paths met each requested marginal tolerance, and the final
plans and regularized objectives agreed to near floating-point roundoff. The
elapsed-time columns are retained as observational measurements, not a speed
guarantee; the certified claim is no greater iteration work at matched answer
quality on this deterministic fixture.

`verify_solver_protocol_state_roundtrip.R` separately saves the public problem
and state objects, launches a fresh installed-package R process, reloads them,
and requires a newly certified warm solve with plan parity.
