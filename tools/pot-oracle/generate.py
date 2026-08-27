#!/usr/bin/env python3
"""Generate replayable live-POT evidence for rfugw differential checks.

This program never writes package fixtures.  It emits a self-contained JSON
receipt that the R comparator consumes from a temporary CI directory.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import pathlib
import platform
import sys
import warnings
from collections.abc import Callable
from datetime import datetime, timezone
from typing import Any

import numpy as np
import ot
import scipy

SCHEMA_VERSION = "rfugw-pot-oracle-v1"
BASELINE_POT_VERSION = "0.9.7.post1"


def to_jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    return value


def stable_digest(value: Any) -> str:
    encoded = json.dumps(
        to_jsonable(value), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def gkl(value: np.ndarray, reference: np.ndarray) -> float:
    value = np.asarray(value, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if np.any((value > 0) & (reference == 0)):
        return float("inf")
    positive = value > 0
    return float(
        np.sum(value[positive] * np.log(value[positive] / reference[positive]))
        - np.sum(value)
        + np.sum(reference)
    )


def product_kl(plan: np.ndarray, source: np.ndarray, target: np.ndarray) -> float:
    return gkl(plan, np.outer(source, target))


def capture_call(call: Callable[[], Any]) -> tuple[Any, dict[str, Any]]:
    stdout = io.StringIO()
    stderr = io.StringIO()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            result = call()
    return result, {
        "warnings": [
            {
                "category": item.category.__name__,
                "message": str(item.message),
            }
            for item in caught
        ],
        "stdout": stdout.getvalue().strip(),
        "stderr": stderr.getvalue().strip(),
    }


def normalize_cost(cost: np.ndarray) -> np.ndarray:
    maximum = float(np.max(cost))
    return cost / maximum if maximum > 0 else cost


def make_linear_case(seed: int) -> dict[str, Any]:
    rng = np.random.RandomState(seed)
    ns, nt = 5, 6
    source_points = rng.normal(size=(ns, 2))
    target_points = rng.normal(size=(nt, 2))
    cost = normalize_cost(ot.dist(source_points, target_points))
    source = rng.dirichlet(np.ones(ns))
    target = rng.dirichlet(np.ones(nt))
    epsilon = 0.12
    rho = (1.7, 2.4)

    exact, exact_diag = capture_call(
        lambda: ot.emd(source, target, cost, numItermax=100000)
    )
    scaling, scaling_diag = capture_call(
        lambda: ot.sinkhorn(
            source,
            target,
            cost,
            reg=epsilon,
            method="sinkhorn",
            numItermax=5000,
            stopThr=1e-10,
        )
    )
    logarithmic, logarithmic_diag = capture_call(
        lambda: ot.sinkhorn(
            source,
            target,
            cost,
            reg=epsilon,
            method="sinkhorn_log",
            numItermax=5000,
            stopThr=1e-10,
        )
    )
    unbalanced, unbalanced_diag = capture_call(
        lambda: ot.unbalanced.sinkhorn_unbalanced(
            source,
            target,
            cost,
            reg=epsilon,
            reg_m=rho,
            reg_type="kl",
            numItermax=5000,
            stopThr=1e-10,
        )
    )

    inputs = {"cost": cost, "source": source, "target": target}
    outputs = {
        "exact_plan": exact,
        "exact_linear_cost": np.sum(cost * exact),
        "scaling_plan": scaling,
        "scaling_linear_cost": np.sum(cost * scaling),
        "scaling_regularized_objective": np.sum(cost * scaling)
        + epsilon * product_kl(scaling, source, target),
        "log_plan": logarithmic,
        "log_linear_cost": np.sum(cost * logarithmic),
        "log_regularized_objective": np.sum(cost * logarithmic)
        + epsilon * product_kl(logarithmic, source, target),
        "unbalanced_plan": unbalanced,
        "unbalanced_linear_cost": np.sum(cost * unbalanced),
        "unbalanced_mass": np.sum(unbalanced),
        "unbalanced_regularized_objective": np.sum(cost * unbalanced)
        + epsilon * product_kl(unbalanced, source, target)
        + rho[0] * gkl(np.sum(unbalanced, axis=1), source)
        + rho[1] * gkl(np.sum(unbalanced, axis=0), target),
    }
    return make_case(
        "balanced_linear",
        seed,
        "strict",
        "ot.emd; ot.sinkhorn; ot.unbalanced.sinkhorn_unbalanced",
        inputs,
        {"epsilon": epsilon, "rho": rho, "max_iter": 5000, "tol": 1e-10},
        outputs,
        {
            "exact": exact_diag,
            "scaling": scaling_diag,
            "log": logarithmic_diag,
            "unbalanced": unbalanced_diag,
        },
        atol=2e-8,
        rtol=2e-7,
    )


def make_partial_linear_case(seed: int) -> dict[str, Any]:
    rng = np.random.RandomState(seed + 1000)
    cost = normalize_cost(rng.uniform(0.05, 1.0, size=(4, 5)))
    source = np.array([0.1, 0.2, 0.3, 0.4])
    target = np.array([0.1, 0.15, 0.2, 0.25, 0.3])
    mass = 0.67
    epsilon = 0.18
    discard_penalty = 0.34

    exact, exact_diag = capture_call(
        lambda: ot.partial.partial_wasserstein(source, target, cost, m=mass)
    )
    scaling, scaling_diag = capture_call(
        lambda: ot.partial.entropic_partial_wasserstein(
            source,
            target,
            cost,
            reg=epsilon,
            m=mass,
            method="sinkhorn",
            numItermax=10000,
            stopThr=1e-12,
        )
    )
    logarithmic, logarithmic_diag = capture_call(
        lambda: ot.partial.entropic_partial_wasserstein(
            source,
            target,
            cost,
            reg=epsilon,
            m=mass,
            method="sinkhorn_log",
            numItermax=10000,
            stopThr=1e-12,
        )
    )
    penalized, penalized_diag = capture_call(
        lambda: ot.partial.partial_wasserstein_lagrange(
            source,
            target,
            cost,
            reg_m=2 * discard_penalty,
        )
    )

    def partial_entropy_objective(plan: np.ndarray) -> float:
        positive = plan > 0
        entropy_minus_one = np.sum(plan[positive] * np.log(plan[positive])) - np.sum(
            plan
        )
        return float(np.sum(cost * plan) + epsilon * entropy_minus_one)

    penalized_mass = float(np.sum(penalized))
    penalized_objective = float(
        np.sum(cost * penalized)
        + discard_penalty * (np.sum(source) + np.sum(target) - 2 * penalized_mass)
    )
    inputs = {"cost": cost, "source": source, "target": target}
    outputs = {
        "exact_plan": exact,
        "exact_linear_cost": np.sum(cost * exact),
        "scaling_plan": scaling,
        "scaling_objective": partial_entropy_objective(scaling),
        "log_plan": logarithmic,
        "log_objective": partial_entropy_objective(logarithmic),
        "penalized_plan": penalized,
        "penalized_mass": penalized_mass,
        "penalized_objective": penalized_objective,
    }
    return make_case(
        "partial_linear",
        seed,
        "strict",
        "ot.partial.partial_wasserstein; ot.partial.entropic_partial_wasserstein; ot.partial.partial_wasserstein_lagrange",
        inputs,
        {
            "mass": mass,
            "epsilon": epsilon,
            "discard_penalty": discard_penalty,
            "max_iter": 10000,
            "tol": 1e-12,
        },
        outputs,
        {
            "exact": exact_diag,
            "scaling": scaling_diag,
            "log": logarithmic_diag,
            "penalized": penalized_diag,
        },
        atol=2e-8,
        rtol=2e-7,
    )


def relational_inputs(seed: int, asymmetric: bool) -> dict[str, np.ndarray]:
    rng = np.random.RandomState(seed + 2000)
    ns, nt = 6, 5
    source_points = rng.normal(size=(ns, 3))
    target_points = rng.normal(size=(nt, 3))
    source_features = rng.normal(size=(ns, 2))
    target_features = rng.normal(size=(nt, 2))
    source_cost = normalize_cost(ot.dist(source_points, source_points))
    target_cost = normalize_cost(ot.dist(target_points, target_points))
    if asymmetric:
        source_cost = normalize_cost(
            source_cost + 0.07 * rng.uniform(size=source_cost.shape)
        )
        target_cost = normalize_cost(
            target_cost + 0.07 * rng.uniform(size=target_cost.shape)
        )
        np.fill_diagonal(source_cost, 0)
        np.fill_diagonal(target_cost, 0)
    feature_cost = normalize_cost(ot.dist(source_features, target_features))
    return {
        "source_cost": source_cost,
        "target_cost": target_cost,
        "feature_cost": feature_cost,
        "source": ot.unif(ns),
        "target": ot.unif(nt),
    }


def make_balanced_gromov_case(seed: int, asymmetric: bool) -> dict[str, Any]:
    data = relational_inputs(seed, asymmetric)
    c1, c2, feature = (
        data["source_cost"],
        data["target_cost"],
        data["feature_cost"],
    )
    source, target = data["source"], data["target"]
    init = np.outer(source, target)
    epsilon = 0.16
    alpha = 0.58

    gw_exact, d1 = capture_call(
        lambda: ot.gromov.gromov_wasserstein(
            c1,
            c2,
            source,
            target,
            symmetric=not asymmetric,
            G0=init,
            max_iter=500,
            tol_rel=1e-9,
            tol_abs=1e-9,
            log=True,
        )
    )
    fgw_exact, d2 = capture_call(
        lambda: ot.gromov.fused_gromov_wasserstein(
            feature,
            c1,
            c2,
            source,
            target,
            symmetric=not asymmetric,
            alpha=alpha,
            G0=init,
            max_iter=500,
            tol_rel=1e-9,
            tol_abs=1e-9,
            log=True,
        )
    )
    gw_entropic, d3 = capture_call(
        lambda: ot.gromov.entropic_gromov_wasserstein(
            c1,
            c2,
            source,
            target,
            symmetric=not asymmetric,
            G0=init,
            epsilon=epsilon,
            solver="PGD",
            max_iter=1000,
            tol=1e-9,
            numItermax=5000,
            stopThr=1e-12,
            log=True,
        )
    )
    fgw_entropic, d4 = capture_call(
        lambda: ot.gromov.entropic_fused_gromov_wasserstein(
            feature,
            c1,
            c2,
            source,
            target,
            symmetric=not asymmetric,
            alpha=alpha,
            G0=init,
            epsilon=epsilon,
            solver="PGD",
            max_iter=1000,
            tol=1e-9,
            numItermax=5000,
            stopThr=1e-12,
            log=True,
        )
    )
    gw_exact_plan, gw_exact_log = gw_exact
    fgw_exact_plan, fgw_exact_log = fgw_exact
    gw_entropic_plan, gw_entropic_log = gw_entropic
    fgw_entropic_plan, fgw_entropic_log = fgw_entropic
    inputs = {**data, "init_plan": init}
    outputs = {
        "gw_exact_plan": gw_exact_plan,
        "gw_exact_objective": gw_exact_log["gw_dist"],
        "fgw_exact_plan": fgw_exact_plan,
        "fgw_exact_objective": fgw_exact_log["fgw_dist"],
        "gw_entropic_plan": gw_entropic_plan,
        "gw_entropic_objective": gw_entropic_log["gw_dist"],
        "fgw_entropic_plan": fgw_entropic_plan,
        "fgw_entropic_objective": fgw_entropic_log["fgw_dist"],
    }
    return make_case(
        "balanced_gromov",
        seed,
        "strict",
        "ot.gromov.gromov_wasserstein; ot.gromov.fused_gromov_wasserstein; entropic variants",
        inputs,
        {
            "epsilon": epsilon,
            "alpha": alpha,
            "symmetric": not asymmetric,
            "max_iter": 1000,
            "tol": 1e-9,
        },
        outputs,
        {"gw_exact": d1, "fgw_exact": d2, "gw_entropic": d3, "fgw_entropic": d4},
        atol=8e-7 if asymmetric else 2e-7,
        rtol=2e-6,
    )


def make_partial_gromov_case(seed: int, asymmetric: bool) -> dict[str, Any]:
    data = relational_inputs(seed + 10, asymmetric)
    c1, c2, feature = (
        data["source_cost"],
        data["target_cost"],
        data["feature_cost"],
    )
    source, target = data["source"], data["target"]
    mass = 0.68
    alpha = 0.61
    epsilon = 0.22
    init = np.outer(source, target) * mass

    pgw, d1 = capture_call(
        lambda: ot.gromov.partial_gromov_wasserstein(
            c1,
            c2,
            source,
            target,
            m=mass,
            symmetric=not asymmetric,
            G0=init,
            numItermax=500,
            tol=1e-9,
            log=True,
        )
    )
    pfgw, d2 = capture_call(
        lambda: ot.gromov.partial_fused_gromov_wasserstein(
            feature,
            c1,
            c2,
            source,
            target,
            m=mass,
            alpha=alpha,
            symmetric=not asymmetric,
            G0=init,
            numItermax=500,
            tol=1e-9,
            log=True,
        )
    )
    epgw, d3 = capture_call(
        lambda: ot.gromov.entropic_partial_gromov_wasserstein(
            c1,
            c2,
            source,
            target,
            reg=epsilon,
            m=mass,
            symmetric=not asymmetric,
            G0=init,
            numItermax=300,
            tol=1e-8,
            log=True,
        )
    )
    epfgw, d4 = capture_call(
        lambda: ot.gromov.entropic_partial_fused_gromov_wasserstein(
            feature,
            c1,
            c2,
            source,
            target,
            reg=epsilon,
            m=mass,
            alpha=alpha,
            symmetric=not asymmetric,
            G0=init,
            numItermax=300,
            tol=1e-8,
            log=True,
        )
    )
    pgw_plan, pgw_log = pgw
    pfgw_plan, pfgw_log = pfgw
    epgw_plan, epgw_log = epgw
    epfgw_plan, epfgw_log = epfgw
    return make_case(
        "partial_gromov",
        seed,
        "strict_with_exception",
        "ot.gromov partial GW and partial FGW APIs",
        {**data, "init_plan": init},
        {
            "mass": mass,
            "alpha": alpha,
            "epsilon": epsilon,
            "symmetric": not asymmetric,
            "max_iter": 500,
            "tol": 1e-9,
        },
        {
            "partial_gw_plan": pgw_plan,
            "partial_gw_objective": pgw_log["partial_gw_dist"],
            "partial_fgw_plan": pfgw_plan,
            "partial_fgw_objective": pfgw_log["partial_fgw_dist"],
            "entropic_partial_gw_plan": epgw_plan,
            "entropic_partial_gw_objective": epgw_log["partial_gw_dist"],
            "entropic_partial_fgw_plan_monitor": epfgw_plan,
            "entropic_partial_fgw_objective_monitor": epfgw_log["partial_fgw_dist"],
        },
        {
            "partial_gw": d1,
            "partial_fgw": d2,
            "entropic_partial_gw": d3,
            "entropic_partial_fgw": d4,
        },
        atol=2e-6,
        rtol=5e-6,
        exception={
            "path": "entropic_partial_fused_gromov_wasserstein",
            "policy": "independent_rfugw_objective_authoritative",
            "reason": "POT 0.9.7.post1 adds a scalar feature term to the entropic partial-FGW gradient instead of the feature-cost matrix.",
        },
    )


def make_semirelaxed_case(seed: int, asymmetric: bool) -> dict[str, Any]:
    data = relational_inputs(seed + 20, asymmetric)
    c1, c2, feature = (
        data["source_cost"],
        data["target_cost"],
        data["feature_cost"],
    )
    source = data["source"]
    init = np.outer(source, ot.unif(c2.shape[0]))
    epsilon = 0.17
    alpha = 0.57
    calls = {}
    diagnostics = {}

    calls["gw_exact"], diagnostics["gw_exact"] = capture_call(
        lambda: ot.gromov.semirelaxed_gromov_wasserstein(
            c1,
            c2,
            p=source,
            symmetric=not asymmetric,
            G0=init,
            max_iter=500,
            tol_rel=1e-9,
            tol_abs=1e-9,
            random_state=0,
            log=True,
        )
    )
    calls["fgw_exact"], diagnostics["fgw_exact"] = capture_call(
        lambda: ot.gromov.semirelaxed_fused_gromov_wasserstein(
            feature,
            c1,
            c2,
            p=source,
            symmetric=not asymmetric,
            alpha=alpha,
            G0=init,
            max_iter=500,
            tol_rel=1e-9,
            tol_abs=1e-9,
            random_state=0,
            log=True,
        )
    )
    calls["gw_entropic"], diagnostics["gw_entropic"] = capture_call(
        lambda: ot.gromov.entropic_semirelaxed_gromov_wasserstein(
            c1,
            c2,
            p=source,
            symmetric=not asymmetric,
            epsilon=epsilon,
            G0=init,
            max_iter=2000,
            tol=1e-10,
            random_state=0,
            log=True,
        )
    )
    calls["fgw_entropic"], diagnostics["fgw_entropic"] = capture_call(
        lambda: ot.gromov.entropic_semirelaxed_fused_gromov_wasserstein(
            feature,
            c1,
            c2,
            p=source,
            symmetric=not asymmetric,
            epsilon=epsilon,
            alpha=alpha,
            G0=init,
            max_iter=2000,
            tol=1e-10,
            random_state=0,
            log=True,
        )
    )
    outputs: dict[str, Any] = {}
    for key, distance_key in (
        ("gw_exact", "srgw_dist"),
        ("fgw_exact", "srfgw_dist"),
        ("gw_entropic", "srgw_dist"),
        ("fgw_entropic", "srfgw_dist"),
    ):
        plan, log = calls[key]
        outputs[f"{key}_plan"] = plan
        outputs[f"{key}_objective"] = log[distance_key]
    return make_case(
        "semirelaxed_gromov",
        seed,
        "strict",
        "ot.gromov semirelaxed exact and entropic APIs",
        {**data, "init_plan": init},
        {
            "epsilon": epsilon,
            "alpha": alpha,
            "symmetric": not asymmetric,
            "max_iter": 2000,
            "tol": 1e-10,
        },
        outputs,
        diagnostics,
        atol=2e-6,
        rtol=5e-6,
    )


def make_ti_uot_case(seed: int) -> dict[str, Any]:
    rng = np.random.RandomState(seed + 3000)
    cost = normalize_cost(rng.uniform(size=(3, 4)))
    source = np.array([0.4, 1.1, 0.7])
    target = np.array([0.2, 0.7, 0.6, 1.4])
    epsilon = 0.31
    rho = (1.4, 2.3)
    classic, d1 = capture_call(
        lambda: ot.unbalanced.sinkhorn_unbalanced(
            source,
            target,
            cost,
            reg=epsilon,
            reg_m=rho,
            reg_type="kl",
            method="sinkhorn",
            numItermax=20000,
            stopThr=1e-12,
        )
    )
    specialized, d2 = capture_call(
        lambda: ot.unbalanced.sinkhorn_unbalanced_translation_invariant(
            source,
            target,
            cost,
            reg=epsilon,
            reg_m=rho,
            reg_type="kl",
            numItermax=20000,
            stopThr=1e-12,
        )
    )

    def objective(plan: np.ndarray) -> float:
        return float(
            np.sum(cost * plan)
            + epsilon * product_kl(plan, source, target)
            + rho[0] * gkl(np.sum(plan, axis=1), source)
            + rho[1] * gkl(np.sum(plan, axis=0), target)
        )

    return make_case(
        "translation_invariant_uot",
        seed,
        "strict_with_exception",
        "ot.unbalanced.sinkhorn_unbalanced (reference estimand); translation-invariant POT solver monitored",
        {"cost": cost, "source": source, "target": target},
        {"epsilon": epsilon, "rho": rho, "max_iter": 20000, "tol": 1e-12},
        {
            "reference_plan": classic,
            "reference_objective": objective(classic),
            "specialized_plan_monitor": specialized,
            "specialized_objective_monitor": objective(specialized),
        },
        {"classic": d1, "translation_invariant": d2},
        atol=2e-7,
        rtol=2e-6,
        exception={
            "path": "ot.unbalanced.sinkhorn_unbalanced_translation_invariant with unequal reg_m",
            "policy": "ordinary_pot_uot_and_independent_objective_authoritative",
            "reason": "The specialized POT path does not reach the ordinary solver objective for unequal source and target marginal penalties.",
        },
    )


def make_sinkhorn_divergence_case(seed: int) -> dict[str, Any]:
    rng = np.random.RandomState(seed + 4000)
    source_points = rng.normal(size=(4, 2))
    target_points = rng.normal(size=(5, 2))
    source = np.array([0.1, 0.2, 0.3, 0.4])
    target = np.array([0.1, 0.15, 0.2, 0.25, 0.3])
    epsilon = 1.5
    cross_cost = ot.dist(source_points, target_points)
    source_cost = ot.dist(source_points, source_points)
    target_cost = ot.dist(target_points, target_points)
    diagnostics = {}
    plans = {}
    for key, a, b, cost in (
        ("cross", source, target, cross_cost),
        ("source_self", source, source, source_cost),
        ("target_self", target, target, target_cost),
    ):
        plans[key], diagnostics[key] = capture_call(
            lambda a=a, b=b, cost=cost: ot.sinkhorn(
                a,
                b,
                cost,
                reg=epsilon,
                method="sinkhorn_log",
                numItermax=20000,
                stopThr=1e-9,
            )
        )
    values = {
        "cross": np.sum(cross_cost * plans["cross"])
        + epsilon * product_kl(plans["cross"], source, target),
        "source_self": np.sum(source_cost * plans["source_self"])
        + epsilon * product_kl(plans["source_self"], source, source),
        "target_self": np.sum(target_cost * plans["target_self"])
        + epsilon * product_kl(plans["target_self"], target, target),
    }
    divergence = values["cross"] - 0.5 * (values["source_self"] + values["target_self"])
    return make_case(
        "sinkhorn_divergence",
        seed,
        "strict",
        "three ot.sinkhorn solves with product-reference KL recomputation",
        {
            "source_points": source_points,
            "target_points": target_points,
            "source": source,
            "target": target,
            "cross_cost": cross_cost,
            "source_cost": source_cost,
            "target_cost": target_cost,
        },
        {"epsilon": epsilon, "power": 2, "max_iter": 20000, "tol": 1e-9},
        {"plans": plans, "component_values": values, "divergence": divergence},
        diagnostics,
        atol=3e-8,
        rtol=3e-7,
    )


def make_fugw_case(seed: int) -> dict[str, Any]:
    rng = np.random.RandomState(seed + 5000)
    nx, ny = 6, 7
    source_points = rng.normal(size=(nx, 3))
    target_points = rng.normal(size=(ny, 3))
    source_cost = normalize_cost(ot.dist(source_points, source_points))
    target_cost = normalize_cost(ot.dist(target_points, target_points))
    source = ot.unif(nx)
    target = ot.unif(ny)
    feature_cost = normalize_cost(rng.uniform(0.05, 1.0, size=(nx, ny)))
    rho = (35.0, 19.0)
    epsilon = 0.025
    alpha = 0.53
    result, diagnostic = capture_call(
        lambda: ot.gromov.fused_unbalanced_gromov_wasserstein(
            source_cost,
            target_cost,
            wx=source,
            wy=target,
            reg_marginals=rho,
            epsilon=epsilon,
            divergence="kl",
            unbalanced_solver="sinkhorn",
            alpha=alpha,
            M=feature_cost,
            init_pi=None,
            max_iter=70,
            tol=1e-9,
            max_iter_ot=700,
            tol_ot=1e-9,
            log=True,
        )
    )
    sample_plan, feature_plan, log = result
    return make_case(
        "fugw",
        seed,
        "strict",
        "ot.gromov.fused_unbalanced_gromov_wasserstein",
        {
            "source_cost": source_cost,
            "target_cost": target_cost,
            "feature_cost": feature_cost,
            "source": source,
            "target": target,
        },
        {
            "rho": rho,
            "epsilon": epsilon,
            "alpha": alpha,
            "max_iter": 70,
            "tol": 1e-9,
            "max_iter_ot": 700,
            "tol_ot": 1e-9,
        },
        {
            "sample_plan": sample_plan,
            "feature_plan": feature_plan,
            "objective": log["fugw_cost"],
        },
        {"fugw": diagnostic},
        atol=8e-4,
        rtol=3e-3,
    )


def make_ucoot_case(seed: int) -> dict[str, Any]:
    rng = np.random.RandomState(seed + 6000)
    source = rng.normal(size=(6, 5))
    target = rng.normal(size=(7, 4))
    rho = (12.0, 8.0)
    epsilon = (0.06, 0.04)
    result, diagnostic = capture_call(
        lambda: ot.gromov.unbalanced_co_optimal_transport(
            source,
            target,
            reg_marginals=rho,
            epsilon=epsilon,
            divergence="kl",
            unbalanced_solver="sinkhorn",
            alpha=(0.0, 0.0),
            max_iter=60,
            tol=1e-8,
            max_iter_ot=700,
            tol_ot=1e-9,
            log=True,
        )
    )
    sample_plan, feature_plan, log = result
    return make_case(
        "ucoot",
        seed,
        "strict",
        "ot.gromov.unbalanced_co_optimal_transport",
        {"source": source, "target": target},
        {
            "rho": rho,
            "epsilon": epsilon,
            "max_iter": 60,
            "tol": 1e-8,
            "max_iter_ot": 700,
            "tol_ot": 1e-9,
        },
        {
            "sample_plan": sample_plan,
            "feature_plan": feature_plan,
            "objective": log["ucoot_cost"],
        },
        {"ucoot": diagnostic},
        atol=2e-2,
        rtol=3e-2,
    )


def square_gw_objective(
    source_cost: np.ndarray, target_cost: np.ndarray, plan: np.ndarray
) -> float:
    loss = (source_cost[:, :, None, None] - target_cost[None, None, :, :]) ** 2
    product = plan[:, None, :, None] * plan[None, :, None, :]
    return float(np.sum(loss * product))


def make_sampled_quality_case(seed: int) -> dict[str, Any]:
    rng = np.random.RandomState(seed + 7000)
    n = 12
    source_points = rng.normal(size=(n, 3))
    target_points = rng.normal(size=(n, 3))
    source_cost = normalize_cost(ot.dist(source_points, source_points))
    target_cost = normalize_cost(ot.dist(target_points, target_points))
    source = ot.unif(n)
    target = ot.unif(n)
    epsilon = 0.5
    max_iter = 80
    random_states = [int(seed % (2**31 - 1)) + offset for offset in range(6)]
    init = np.outer(source, target)

    dense_result, dense_diagnostic = capture_call(
        lambda: ot.gromov.entropic_gromov_wasserstein(
            source_cost,
            target_cost,
            source,
            target,
            epsilon=epsilon,
            solver="PGD",
            G0=init,
            max_iter=max_iter,
            tol=1e-9,
            numItermax=5000,
            stopThr=1e-12,
            log=True,
        )
    )

    def sampled(
        budget: tuple[int, int], random_state: int
    ) -> tuple[Any, dict[str, Any]]:
        return capture_call(
            lambda: ot.gromov.sampled_gromov_wasserstein(
                source_cost,
                target_cost,
                source,
                target,
                loss_fun=lambda a, b: (a - b) ** 2,
                nb_samples_grad=budget,
                epsilon=epsilon,
                max_iter=max_iter,
                random_state=random_state,
                log=True,
            )
        )

    tiny_results = []
    tiny_diagnostics = []
    for random_state in random_states:
        result, diagnostic = sampled((2, 1), random_state)
        tiny_results.append(result)
        tiny_diagnostics.append(diagnostic)
    high_result, high_diagnostic = sampled((8, 2), random_states[0])
    dense_plan, dense_log = dense_result
    tiny_plans = [result[0] for result in tiny_results]
    tiny_logs = [result[1] for result in tiny_results]
    high_plan, high_log = high_result
    tiny_objectives = [
        square_gw_objective(source_cost, target_cost, plan) for plan in tiny_plans
    ]
    tiny_distances = [float(np.linalg.norm(plan - dense_plan)) for plan in tiny_plans]
    objectives = {
        "dense": square_gw_objective(source_cost, target_cost, dense_plan),
        "tiny": tiny_objectives,
        "high": square_gw_objective(source_cost, target_cost, high_plan),
    }
    plan_distances = {
        "tiny_to_dense": tiny_distances,
        "high_to_dense": float(np.linalg.norm(high_plan - dense_plan)),
    }
    return make_case(
        "sampled_gromov_quality",
        seed,
        "strict_quality_not_plan_parity",
        "ot.gromov.sampled_gromov_wasserstein; entropic_gromov_wasserstein",
        {
            "source_cost": source_cost,
            "target_cost": target_cost,
            "source": source,
            "target": target,
            "init_plan": init,
        },
        {
            "epsilon": epsilon,
            "max_iter": max_iter,
            "random_states": random_states,
            "tiny_budget": (2, 1),
            "pot_high_budget": (8, 2),
            "rfugw_high_budget": (n, n),
        },
        {
            "dense_plan": dense_plan,
            "dense_log_objective": dense_log["gw_dist"],
            "tiny_plans_monitor": tiny_plans,
            "tiny_reported_estimates": [log["gw_dist_estimated"] for log in tiny_logs],
            "high_plan_monitor": high_plan,
            "high_reported_estimate": high_log["gw_dist_estimated"],
            "independent_objectives": objectives,
            "plan_distances": plan_distances,
            "pot_high_beats_median_tiny_objective_gap": abs(
                objectives["high"] - objectives["dense"]
            )
            < float(
                np.median(
                    [abs(value - objectives["dense"]) for value in objectives["tiny"]]
                )
            ),
            "pot_high_beats_median_tiny_plan_distance": plan_distances["high_to_dense"]
            < float(np.median(plan_distances["tiny_to_dense"])),
        },
        {
            "dense": dense_diagnostic,
            **{
                f"tiny_{index + 1}": diagnostic
                for index, diagnostic in enumerate(tiny_diagnostics)
            },
            "high": high_diagnostic,
        },
        atol=2e-6,
        rtol=5e-6,
        exception={
            "path": "ot.gromov.sampled_gromov_wasserstein inner Sinkhorn diagnostics",
            "policy": "quality_distribution_and_feasibility_authoritative",
            "reason": "POT's sampled GW API does not expose inner Sinkhorn iteration or tolerance controls and can emit convergence warnings; exact stochastic plan parity is not certified.",
        },
    )


def make_barycenter_semantics_case(seed: int) -> dict[str, Any]:
    rng = np.random.RandomState(seed + 8000)
    costs = []
    features = []
    weights = []
    for n in (4, 5):
        points = rng.normal(size=(n, 2))
        costs.append(normalize_cost(ot.dist(points, points)))
        features.append(rng.normal(size=(n, 2)))
        weights.append(ot.unif(n))
    n_barycenter = 3
    barycenter_weights = ot.unif(n_barycenter)
    init_points = rng.normal(size=(n_barycenter, 2))
    init_cost = normalize_cost(ot.dist(init_points, init_points))
    init_features = rng.normal(size=(n_barycenter, 2))
    lambdas = (0.4, 0.6)
    epsilon = 0.4

    gw_result, gw_diagnostic = capture_call(
        lambda: ot.gromov.entropic_gromov_barycenters(
            n_barycenter,
            costs,
            ps=weights,
            p=barycenter_weights,
            lambdas=lambdas,
            epsilon=epsilon,
            max_iter=1,
            tol=1e-12,
            init_C=init_cost,
            log=True,
            numItermax=5000,
            stopThr=1e-12,
        )
    )
    fused_result, fused_diagnostic = capture_call(
        lambda: ot.gromov.entropic_fused_gromov_barycenters(
            n_barycenter,
            features,
            costs,
            ps=weights,
            p=barycenter_weights,
            lambdas=lambdas,
            epsilon=epsilon,
            alpha=0.55,
            max_iter=1,
            tol=1e-12,
            init_C=init_cost,
            init_Y=init_features,
            log=True,
            numItermax=5000,
            stopThr=1e-12,
        )
    )
    gw_cost, gw_log = gw_result
    fused_features, fused_cost, fused_log = fused_result
    return make_case(
        "barycenter_semantics",
        seed,
        "strict_off_diagonal_with_zero_diagonal_exception",
        "ot.gromov.entropic_gromov_barycenters; entropic_fused_gromov_barycenters",
        {
            "costs": costs,
            "features": features,
            "weights": weights,
            "barycenter_weights": barycenter_weights,
            "init_cost": init_cost,
            "init_features": init_features,
        },
        {
            "lambdas": lambdas,
            "epsilon": epsilon,
            "alpha": 0.55,
            "max_iter": 1,
            "tol": 1e-12,
            "sinkhorn_max_iter": 5000,
            "sinkhorn_tol": 1e-12,
        },
        {
            "gw_cost": gw_cost,
            "gw_couplings": gw_log["T"],
            "fused_features": fused_features,
            "fused_cost": fused_cost,
            "fused_couplings": fused_log["T"],
        },
        {"gw": gw_diagnostic, "fused": fused_diagnostic},
        atol=3e-5,
        rtol=8e-5,
        exception={
            "path": "POT entropic barycenter output semantics",
            "policy": "closed_form_updates_and_rfugw_zero_diagonal_invariant_are_authoritative",
            "reason": "POT preserves positive entropic self-costs on the barycenter diagonal and its fused routine returns the stale initial Y variable after updating X; rfugw enforces diag(C) = 0 and returns the closed-form feature update.",
        },
    )


def make_case(
    family: str,
    seed: int,
    policy: str,
    oracle_function: str,
    inputs: dict[str, Any],
    params: dict[str, Any],
    outputs: dict[str, Any],
    diagnostics: dict[str, Any],
    *,
    atol: float,
    rtol: float,
    exception: dict[str, Any] | None = None,
) -> dict[str, Any]:
    payload = {
        "family": family,
        "seed": seed,
        "policy": policy,
        "oracle_function": oracle_function,
        "inputs": to_jsonable(inputs),
        "params": to_jsonable(params),
        "outputs": to_jsonable(outputs),
        "diagnostics": to_jsonable(diagnostics),
        "tolerance": {"atol": atol, "rtol": rtol},
    }
    if exception is not None:
        payload["exception"] = exception
    payload["input_sha256"] = stable_digest(payload["inputs"])
    payload["output_sha256"] = stable_digest(payload["outputs"])
    return payload


def strict_diagnostic_violations(case: dict[str, Any]) -> list[str]:
    violations = []
    ignored = {
        "partial_gromov": {"entropic_partial_fgw"},
        "translation_invariant_uot": {"translation_invariant"},
        "sampled_gromov_quality": set(case["diagnostics"]),
    }.get(case["family"], set())
    for name, diagnostic in case["diagnostics"].items():
        if name in ignored:
            continue
        if diagnostic.get("warnings"):
            violations.append(f"{name}: warning")
        if diagnostic.get("stdout"):
            violations.append(f"{name}: stdout")
        if diagnostic.get("stderr"):
            violations.append(f"{name}: stderr")
    return violations


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--profile", choices=("pr", "nightly", "release"), default="pr")
    parser.add_argument("--channel", choices=("baseline", "latest"), default="baseline")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.channel == "baseline" and ot.__version__ != BASELINE_POT_VERSION:
        raise RuntimeError(
            f"baseline channel requires POT {BASELINE_POT_VERSION}, found {ot.__version__}"
        )

    seeds = [20260826]
    if args.profile in ("nightly", "release"):
        seeds.extend([20260827, 20260828])

    cases: list[dict[str, Any]] = []
    for index, seed in enumerate(seeds):
        asymmetric = index > 0
        cases.extend(
            [
                make_linear_case(seed),
                make_partial_linear_case(seed),
                make_balanced_gromov_case(seed, asymmetric),
                make_partial_gromov_case(seed, asymmetric),
                make_semirelaxed_case(seed, asymmetric),
                make_ti_uot_case(seed),
                make_sinkhorn_divergence_case(seed),
                make_fugw_case(seed),
                make_ucoot_case(seed),
                make_sampled_quality_case(seed),
                make_barycenter_semantics_case(seed),
            ]
        )

    for case in cases:
        violations = strict_diagnostic_violations(case)
        if violations:
            raise RuntimeError(
                f"strict POT path in {case['family']} emitted diagnostics: {', '.join(violations)}"
            )

    receipt = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "profile": args.profile,
        "channel": args.channel,
        "baseline_pot_version": BASELINE_POT_VERSION,
        "runtime": {
            "pot_version": ot.__version__,
            "numpy_version": np.__version__,
            "scipy_version": scipy.__version__,
            "python_version": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "platform": platform.platform(),
        },
        "case_count": len(cases),
        "cases": cases,
    }
    receipt["cases_sha256"] = stable_digest(cases)
    output = pathlib.Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(
        f"wrote {output} with {len(cases)} cases "
        f"(POT {ot.__version__}, profile={args.profile}, channel={args.channel})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
