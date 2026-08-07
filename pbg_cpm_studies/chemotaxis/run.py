"""Canonical runner for the chemotactic-recruitment slice.

Runs the three conditions (baseline / inhibited / adversarial), equilibrates the
secreted gradient, then steps the CPM world recording the recruitment index and
per-frame cell positions. Writes one series-JSON per condition under
``workspace/chemotaxis_data/results/`` (the ``runs[].emitter`` convention) plus a
combined ``recruitment_summary.json`` the visualizations read.

Run (from the worktree root, using the viva-cpm venv)::

    PYTHONPATH=. python -m pbg_cpm_studies.chemotaxis.run
"""
from __future__ import annotations

import json
import os

from cpm import load_world

from ..composites.chemotaxis import (
    build_spec, meta, CUE_RATE, CHEMO_LAMBDA, SOURCE_TYPE, RESPONDER_TYPE,
)
from . import metrics as M

WARMUP = 600           # advance_fields sweeps to equilibrate the gradient
SWEEPS = 500           # MC sweeps of the coupled cell+field dynamics
SAMPLE = 25            # record every SAMPLE sweeps
RADIUS = meta()["recruitment_radius"]

CONDITIONS = {
    "recruitment-baseline":    dict(cue_rate=CUE_RATE, chemo_lambda=CHEMO_LAMBDA),
    "recruitment-inhibited":   dict(cue_rate=CUE_RATE, chemo_lambda=0.0),
    "recruitment-adversarial": dict(cue_rate=0.0,      chemo_lambda=CHEMO_LAMBDA),
}


def _frame(world):
    coms = list(world.cell_coms())
    types = list(world.cell_types())
    cells = [[cid, int(types[cid]), round(coms[cid][0], 1), round(coms[cid][1], 1)]
             for cid in range(1, len(types)) if int(types[cid]) != 0]
    return {
        "recruitment_index": round(M.recruitment_index_from_coms(
            coms, types, radius=RADIUS), 4),
        "mean_distance": round(M.mean_distance_to_source_from_coms(coms, types), 2),
        "cells": cells,
    }


def run_condition(name, *, cue_rate, chemo_lambda, seed=17):
    world = load_world(build_spec(cue_rate=cue_rate, chemo_lambda=chemo_lambda, seed=seed))
    for _ in range(WARMUP):
        world.advance_fields(1)
    series_mcs, recruitment, distance, frames = [], [], [], []
    f0 = _frame(world)
    series_mcs.append(0); recruitment.append(f0["recruitment_index"])
    distance.append(f0["mean_distance"]); frames.append({"mcs": 0, **f0})
    for s in range(SAMPLE, SWEEPS + 1, SAMPLE):
        world.step(SAMPLE)
        f = _frame(world)
        series_mcs.append(s); recruitment.append(f["recruitment_index"])
        distance.append(f["mean_distance"]); frames.append({"mcs": s, **f})
    approach = M.mean_approach(distance[0], distance[-1])
    return {
        "condition": name,
        "params": {"cue_rate": cue_rate, "chemo_lambda": chemo_lambda,
                   "radius": RADIUS, "warmup": WARMUP, "sweeps": SWEEPS, "seed": seed},
        "series_mcs": series_mcs,
        "recruitment_index": recruitment,
        "mean_distance": distance,
        "final_recruitment_index": recruitment[-1],
        "mean_approach": round(approach, 3),
        "frames": frames,
    }


def main(out_dir=None):
    here = os.path.dirname(os.path.abspath(__file__))
    ws = os.path.abspath(os.path.join(here, "..", "..", "workspace"))
    out_dir = out_dir or os.path.join(ws, "chemotaxis_data", "results")
    os.makedirs(out_dir, exist_ok=True)
    summary = {"radius": RADIUS, "source_type": SOURCE_TYPE,
               "responder_type": RESPONDER_TYPE, "conditions": {}}
    print(f"{'condition':26s} {'final recruit':>13s} {'approach':>9s}")
    for name, kw in CONDITIONS.items():
        res = run_condition(name, **kw)
        with open(os.path.join(out_dir, f"{name}_series.json"), "w") as fh:
            json.dump(res, fh)
        summary["conditions"][name] = {
            "params": res["params"],
            "final_recruitment_index": res["final_recruitment_index"],
            "mean_approach": res["mean_approach"],
            "series_mcs": res["series_mcs"],
            "recruitment_index": res["recruitment_index"],
        }
        print(f"{name:26s} {res['final_recruitment_index']:13.3f} {res['mean_approach']:9.3f}")
    with open(os.path.join(out_dir, "recruitment_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\nwrote {len(CONDITIONS)} series + summary to {out_dir}")
    return summary


if __name__ == "__main__":
    main()
