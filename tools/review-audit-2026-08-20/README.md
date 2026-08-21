# rfugw review verification

This directory contains fresh-process reproduction scripts and the final
evidence ledger for the source audit reported against commit
`41b41ed713e68fe36a69ad088f34385e78f5699e`.

The current package source was verified byte-for-byte unchanged from that
commit outside `.mote/` tracker metadata before these probes were designed.

No package implementation changes are made by this audit.

## Contents

- `REPORT.md`: dispositions, evidence summary, and remediation plan.
- `01_asymmetric_gradient.R`: brute-force objective and gradient law.
- `02_transport_simplex.R`: exact-OT differential fuzzing and partial status
  propagation.
- `03_mixed_precision_contract.R`: requested/effective tolerance and precision
  behavior.
- `04_status_contract.R`: outer/inner convergence and result-helper laws.
- `05_ot_kl.R`: generalized KL and zero-support behavior.
- `06_asymmetric_solver_paths.R`: independent one-step solver differentials.
- `07_secondary_contracts.R`: validation, scale, and log-domain claims.
- `run-all.sh`: fresh-process runner for the seven probes.

## Run

Install the audited source into an isolated library, then run:

```sh
RFUGW_AUDIT_LIB=/path/to/library \
  tools/review-audit-2026-08-20/run-all.sh
```

`02_transport_simplex.R` and `06_asymmetric_solver_paths.R` require `lpSolve`.
