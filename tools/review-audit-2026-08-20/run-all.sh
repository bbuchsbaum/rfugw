#!/bin/sh
set -eu

audit_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

for script in \
  01_asymmetric_gradient.R \
  02_transport_simplex.R \
  03_mixed_precision_contract.R \
  04_status_contract.R \
  05_ot_kl.R \
  06_asymmetric_solver_paths.R \
  07_secondary_contracts.R
do
  Rscript --vanilla "$audit_dir/$script"
done
