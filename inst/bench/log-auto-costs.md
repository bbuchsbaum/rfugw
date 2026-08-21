# Log/auto correctness and cost record

`benchmark_log_auto.R` gates timing on agreement between scaling- and
log-domain UOT in a moderate regime (`max(abs(plan_scaling - plan_log)) <=
5e-7` and objective difference `<= 5e-6`). It then records median elapsed time
over three repetitions, result size, certificate status, effective method, and
dispatch criterion.

The checked-in baseline was collected on 2026-08-20 with R 4.5.1 on Darwin
23.3.0 arm64. The source tree was dirty because this certification work was not
committed; these values are engineering evidence, not a release or publication
claim. See `log-auto-baseline.csv` for the machine-readable rows.

At `n = 64`, the genuine R log-domain UOT path took 0.135 seconds versus 0.003
seconds for the C++ scaling path and used 87,664 versus 86,056 result bytes.
The safe auto case selected scaling. The adversarial auto case (dynamic-range
metric 99,958) selected log, converged in 1.274 seconds, and used 87,712 result
bytes. The moderate scaling run's fixed-point certificate remained `max_iter`;
that row is retained transparently as a cost comparison and is not eligible as
quality evidence. The log rows did converge.

Replay from the package root:

```sh
Rscript --vanilla inst/bench/benchmark_log_auto.R \
  inst/bench/log-auto-baseline.csv
```
