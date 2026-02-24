#!/usr/bin/env bash
set -euo pipefail

OUT_CSV="${1:-inst/bench/results/benchmark_sparse_sampled_latest.csv}"
SEED="${2:-123}"
ITERS="${3:-2}"
N_VALUES="${4:-400,800}"
SINKHORN_MAX_ITER="${5:-60}"
SINKHORN_TOL="${6:-1e-6}"
MIN_N="${7:-400}"
MIN_RATIO_SOLVER="${8:-1.0}"
MIN_RATIO_E2E="${9:-1.0}"
MAX_OBJ_GAP="${10:-0.005}"

RFUGW_PYTHON_BIN="${RFUGW_PYTHON:-python3}"

OUT_MD="${OUT_CSV%.csv}_report.md"

echo "== sparse sampled perf gate run =="
echo "python: ${RFUGW_PYTHON_BIN}"
echo "csv: ${OUT_CSV}"
echo "report: ${OUT_MD}"

echo "-- benchmark"
RFUGW_PYTHON="${RFUGW_PYTHON_BIN}" \
  Rscript inst/bench/benchmark_sparse_sampled.R \
    "${ITERS}" \
    "${OUT_CSV}" \
    "${SEED}" \
    "${N_VALUES}" \
    "${SINKHORN_MAX_ITER}" \
    "${SINKHORN_TOL}"

echo "-- report"
Rscript inst/bench/compare_sparse_sampled_benchmark.R \
  "${OUT_CSV}" \
  "${OUT_MD}"

echo "-- gate"
Rscript inst/bench/gate_sparse_sampled_perf.R \
  "${OUT_CSV}" \
  "${MIN_N}" \
  "${MIN_RATIO_SOLVER}" \
  "${MIN_RATIO_E2E}" \
  "${MAX_OBJ_GAP}"

cat <<MSG
Done.
- CSV: ${OUT_CSV}
- Report: ${OUT_MD}
MSG
