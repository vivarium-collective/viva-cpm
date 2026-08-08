"""Workspace-registered behavior-contract evaluators for cpm-studies.

The workbench's ``study_evaluator`` calls ``register_evaluators(registry)`` and
dispatches any ``measure.kind`` it doesn't natively know to the callable we
register here. This is the concrete instance of the Biological Claim Layer's
"backend-neutral behavior contract routed through code": a study declares

    behavior_tests:
      - measure:  {kind: recruitment_index, condition: recruitment-baseline, stat: final}
        pass_if:  {op: gt, value: 0.5}

and the CPM-specific observable extractor below turns the recorded run into a
deterministic PASS/FAIL — the same contract could be re-pointed at another
backend by registering a different extractor for the same ``kind``.
"""
from __future__ import annotations

import json
import os
from typing import Any, Callable

_OPS: dict[str, Callable[[float, float], bool]] = {
    "gt": lambda a, b: a > b,
    "gte": lambda a, b: a >= b,
    "lt": lambda a, b: a < b,
    "lte": lambda a, b: a <= b,
}


def _results_dir(ws_root: Any) -> str | None:
    for cand in (
        os.path.join(str(ws_root), "workspace", "chemotaxis_data", "results"),
        os.path.join(str(ws_root), "chemotaxis_data", "results"),
    ):
        if os.path.isdir(cand):
            return cand
    return None


def _code(result, measured_value, operator, detail):
    return {"result": result, "measured_value": measured_value,
            "evaluated_by": "code", "operator": operator, "detail": detail}


def _recruitment_index(test: dict, reader: Any, ws_root: Any) -> dict:
    """Evaluate a recruitment_index contract against the recorded series-JSON.

    Ignores ``reader`` (this workspace persists canonical runs as series-JSON,
    not the standard run store) and reads the condition's recorded trajectory.
    """
    measure = test.get("measure") or {}
    pass_if = test.get("pass_if") or {}
    condition = measure.get("condition")
    stat = measure.get("stat", "final")
    op = pass_if.get("op")
    threshold = pass_if.get("value")

    if condition is None or op not in _OPS or threshold is None:
        return {"evaluated_by": "agent",
                "reason": "recruitment_index needs measure.condition, pass_if.op in "
                          "{gt,gte,lt,lte}, pass_if.value"}
    rdir = _results_dir(ws_root)
    if rdir is None:
        return {"evaluated_by": "needs_rerun",
                "reason": "no chemotaxis_data/results — run pbg_cpm_studies.chemotaxis.run"}
    path = os.path.join(rdir, f"{condition}_series.json")
    if not os.path.isfile(path):
        return {"evaluated_by": "needs_rerun", "reason": f"missing series for {condition}"}

    with open(path) as fh:
        series = json.load(fh).get("recruitment_index") or []
    if not series:
        return {"evaluated_by": "needs_rerun", "reason": f"empty series for {condition}"}

    value = series[-1] if stat == "final" else (max(series) if stat == "max" else min(series))
    ok = _OPS[op](float(value), float(threshold))
    return _code(
        "PASS" if ok else "FAIL",
        round(float(value), 4),
        f"recruitment_index[{stat}]/{op}",
        f"{stat} recruitment_index = {value:.3f} {op} {threshold} for {condition}",
    )


def register_evaluators(registry: dict) -> None:
    registry["recruitment_index"] = _recruitment_index
