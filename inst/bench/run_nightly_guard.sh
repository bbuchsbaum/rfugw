#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUT_DIR="${1:-$SCRIPT_DIR/results/nightly}"
SEED="${2:-20260222}"
THREADS="${3:-2}"

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export OMP_NUM_THREADS="$THREADS"

cd "$ROOT_DIR"
Rscript "$SCRIPT_DIR/nightly_guard.R" "$OUT_DIR" "$SEED" "$THREADS"
