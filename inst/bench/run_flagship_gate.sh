#!/usr/bin/env bash
# Run the canonical protocol and apply the quality-first CI gate.
# Exit 0 ok, 1 solver/regression, 2 infrastructure.
set -u

if [[ ! -f DESCRIPTION ]]; then
  echo "INFRA: run inst/bench/run_flagship_gate.sh from the package root." >&2
  exit 2
fi

OUT_DIR="${1:-inst/bench/results/current}"
SUITES="${2:-fgw,fugw,semirelaxed}"
SCALE="${3:-pr}"
REPS="${4:-1}"
SEED="${5:-20260816}"
THREADS="${6:-1}"
CAPS="${7:-inst/bench/ci_time_caps.json}"

export RFUGW_PROTOCOL_SCALE="${SCALE}"
export RFUGW_PROTOCOL_SUITES="${SUITES}"

set +e
Rscript inst/bench/run_protocol.R "${REPS}" "${SEED}" "${OUT_DIR}" "${THREADS}" "${SUITES}"
proto_status=$?
set -e

if [[ ! -f "${OUT_DIR}/runs.csv" || ! -f "${OUT_DIR}/meta.json" ]]; then
  echo "INFRA: protocol did not write runs.csv/meta.json (exit ${proto_status})." >&2
  exit 2
fi

set +e
Rscript inst/bench/gate_protocol.R "${OUT_DIR}" "${CAPS}" "${SCALE}"
gate_status=$?
set -e

if [[ ${gate_status} -ne 0 ]]; then
  exit "${gate_status}"
fi
if [[ ${proto_status} -ne 0 ]]; then
  echo "SOLVER: protocol exited ${proto_status} after writing artifacts." >&2
  exit 1
fi
exit 0
