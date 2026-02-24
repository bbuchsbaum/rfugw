#!/usr/bin/env bash
set -euo pipefail

if [[ -f "inst/bench/benchmark_thread_scaling.R" ]]; then
  BENCH_ROOT="inst/bench"
  FIXTURE_DEFAULT="inst/extdata/fixtures"
  OUT_DEFAULT="inst/bench/results/thread_scaling_latest.csv"
else
  BENCH_ROOT="rfugw/inst/bench"
  FIXTURE_DEFAULT="rfugw/inst/extdata/fixtures"
  OUT_DEFAULT="rfugw/inst/bench/results/thread_scaling_latest.csv"
fi

OUT_CSV="${1:-$OUT_DEFAULT}"
ITERS="${2:-3}"
SEED="${3:-42}"
THREADS="${4:-1 2 4 8}"
THREADS="${THREADS//,/ }"
RFUGW_RLIB_PATH="${RFUGW_RLIB:-}"

mkdir -p "$(dirname "$OUT_CSV")"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

if [[ "${RFUGW_SKIP_ACCURACY_GATE:-0}" != "1" ]]; then
  echo "Running strict accuracy gate once before thread sweep"
  RFUGW_RLIB="$RFUGW_RLIB_PATH" \
  RFUGW_FIXTURE_DIR="${RFUGW_FIXTURE_DIR:-$FIXTURE_DEFAULT}" \
  Rscript "$BENCH_ROOT/accuracy_gate.R"
fi

first=1
for t in $THREADS; do
  part_csv="$TMP_DIR/part_${t}.csv"
  echo "Running thread benchmark with OMP_NUM_THREADS=$t"
  RFUGW_RLIB="$RFUGW_RLIB_PATH" \
  RFUGW_SKIP_ACCURACY_GATE=1 \
  OMP_NUM_THREADS="$t" \
  OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  VECLIB_MAXIMUM_THREADS=1 \
  Rscript "$BENCH_ROOT/benchmark_thread_scaling.R" "$part_csv" "$ITERS" "$SEED"

  if [[ $first -eq 1 ]]; then
    cat "$part_csv" > "$OUT_CSV"
    first=0
  else
    tail -n +2 "$part_csv" >> "$OUT_CSV"
  fi
done

echo "Wrote merged thread scaling CSV: $OUT_CSV"
