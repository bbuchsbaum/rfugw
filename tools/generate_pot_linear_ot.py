#!/usr/bin/env python3
"""Generate POT reference fixtures for rfugw linear-OT certification.

Provenance is recorded in the JSON `meta` block. Re-run from the package
root with the project POT venv:

    .venv_pot/bin/python tools/generate_pot_linear_ot.py
"""

from __future__ import annotations

import json
import pathlib
from datetime import datetime, timezone

import numpy as np
import ot

ROOT = pathlib.Path(__file__).resolve().parents[1]
out_path = ROOT / "inst" / "extdata" / "fixtures" / "linear_ot_fixture.json"
out_path.parent.mkdir(parents=True, exist_ok=True)


def tolist(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    return float(x) if np.isscalar(x) else x


rng = np.random.RandomState(20260816)

# Balanced 6 x 7 Sinkhorn / EMD problem
ns, nt = 6, 7
X = rng.normal(size=(ns, 2))
Y = rng.normal(size=(nt, 2))
M = ot.dist(X, Y)
M = M / M.max()
p = rng.dirichlet(np.ones(ns))
q = rng.dirichlet(np.ones(nt))
epsilon = 0.08

T_sink, log_sink = ot.sinkhorn(
    p, q, M, reg=epsilon, numItermax=2000, stopThr=1e-9, log=True
)
T_log, log_log = ot.sinkhorn(
    p, q, M, reg=epsilon, method="sinkhorn_log", numItermax=2000, stopThr=1e-9, log=True
)
T_emd = ot.emd(p, q, M, numItermax=100000)

# Unbalanced KL Sinkhorn
rho = 2.0
T_uot, log_uot = ot.unbalanced.sinkhorn_unbalanced(
    p, q, M, reg=epsilon, reg_m=rho, numItermax=800, stopThr=1e-9, log=True
)

# 2 x 2 analytic assignment (zero diagonal)
M2 = np.array([[0.0, 2.0], [2.0, 0.0]])
p2 = np.array([0.4, 0.6])
q2 = np.array([0.4, 0.6])
T2 = ot.emd(p2, q2, M2)

fixture = {
    "meta": {
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "pot_version": ot.__version__,
        "numpy_version": np.__version__,
        "seed": 20260816,
        "script": "tools/generate_pot_linear_ot.py",
        "notes": (
            "POT ot.sinkhorn / ot.emd / ot.unbalanced.sinkhorn_unbalanced "
            "references for rfugw ot_sinkhorn, ot_emd, ot_sinkhorn_unbalanced. "
            "Reported rfugw ot_dist is unregularized <M, plan>."
        ),
    },
    "balanced": {
        "inputs": {"M": tolist(M), "p": tolist(p), "q": tolist(q)},
        "params": {"epsilon": epsilon, "max_iter": 2000, "tol": 1e-9},
        "outputs": {
            "sinkhorn_plan": tolist(T_sink),
            "sinkhorn_cost": float(np.sum(M * T_sink)),
            "sinkhorn_log_plan": tolist(T_log),
            "sinkhorn_log_cost": float(np.sum(M * T_log)),
            "emd_plan": tolist(T_emd),
            "emd_cost": float(np.sum(M * T_emd)),
        },
    },
    "unbalanced": {
        "inputs": {"M": tolist(M), "p": tolist(p), "q": tolist(q)},
        "params": {"epsilon": epsilon, "rho": rho, "max_iter": 800, "tol": 1e-9},
        "outputs": {
            "plan": tolist(T_uot),
            "ot_dist": float(np.sum(M * T_uot)),
            "mass": float(np.sum(T_uot)),
        },
    },
    "analytic_2x2": {
        "inputs": {"M": tolist(M2), "p": tolist(p2), "q": tolist(q2)},
        "outputs": {
            "emd_plan": tolist(T2),
            "emd_cost": float(np.sum(M2 * T2)),
        },
    },
}

out_path.write_text(json.dumps(fixture, indent=2) + "\n")
print(f"wrote {out_path} (POT {ot.__version__})")
