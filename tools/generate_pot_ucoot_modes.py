#!/usr/bin/env python3
"""POT joint and independent UCOOT / across-spaces objectives."""
import json
import pathlib
import numpy as np
import ot

ROOT = pathlib.Path(__file__).resolve().parents[1]
out_path = ROOT / "inst" / "extdata" / "fixtures" / "ucoot_pot_modes.json"

rng = np.random.RandomState(20260816)
X = rng.normal(size=(7, 5))
Y = rng.normal(size=(6, 4))
params = dict(
    reg_marginals=(10.0, 8.0),
    epsilon=(0.05, 0.03),
    divergence="kl",
    unbalanced_solver="sinkhorn",
    alpha=(0.0, 0.0),
    max_iter=40,
    tol=1e-7,
    max_iter_ot=200,
    tol_ot=1e-7,
    log=True,
)

pi_ind, q_ind, log_ind = ot.gromov.unbalanced_co_optimal_transport(X, Y, **params)
pi_joint, q_joint, log_joint = ot.gromov.fused_unbalanced_across_spaces_divergence(
    X, Y, reg_type="joint", **params
)

out_path.write_text(
    json.dumps(
        {
            "pot_version": ot.__version__,
            "seed": 20260816,
            "inputs": {"X": X.tolist(), "Y": Y.tolist()},
            "params": {
                "reg_marginals": [10.0, 8.0],
                "epsilon": [0.05, 0.03],
                "max_iter": 40,
                "tol": 1e-7,
                "max_iter_ot": 200,
                "tol_ot": 1e-7,
            },
            "independent": {
                "pi_samp": pi_ind.tolist(),
                "pi_feat": q_ind.tolist(),
                "ucoot_cost": float(log_ind["ucoot_cost"]),
            },
            "joint": {
                "pi_samp": pi_joint.tolist(),
                "pi_feat": q_joint.tolist(),
                "ucoot_cost": float(log_joint["ucoot_cost"]),
            },
        },
        indent=2,
    )
)
print("wrote", out_path)
