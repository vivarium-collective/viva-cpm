"""Multi-seed run harness for the receptor-level recruitment composite.

Unlike ``pbg_cpm_studies.chemotaxis.run`` (which uses raw ``cpm.load_world``
and steps the CPM world directly), this module builds and steps the
``recruitment_receptor`` composite via ``process_bigraph.Composite`` -- the
only way the per-responder ``ReceptorSubcell`` (Task 1) actually gets wired
and invoked. ``load_world`` alone bypasses the subcells entirely: the fates
port would never be written and no cell would ever activate.

The recruitment readout counts ALL responders (naive(2) AND activated(3)),
not activated-only -- gating recruitment on "reached activated" would be
circular, since activation and recruitment would both derive from the same
underlying occupancy signal. See ``pbg_cpm_studies.chemotaxis.metrics.
recruitment_index_all_responders``.

Run (from the worktree root, using the viva-cpm venv)::

    PYTHONPATH=. python -c "from pbg_cpm_studies.chemotaxis.run_receptor import run_receptor; \\
        run_receptor(blocked=False, seeds=[17,29,43,61,89], steps=500, kd=2.9); \\
        run_receptor(blocked=True,  seeds=[17,29,43,61,89], steps=500, kd=2.9)"
"""
from __future__ import annotations

import json
import math
import os

import process_bigraph as pb

from ..composites import chemotaxis_receptor as CR
from . import metrics as M

SAMPLE = 25  # record every SAMPLE composite ticks (matches run.py's SAMPLE)
RESPONDER_TYPES = (CR.NAIVE_TYPE, CR.ACTIVATED_TYPE)


def _default_out_dir():
    here = os.path.dirname(os.path.abspath(__file__))
    ws = os.path.abspath(os.path.join(here, "..", "..", "workspace"))
    return os.path.join(ws, "chemotaxis_data", "results")


def _frame(world, responder_ids, radius):
    coms = list(world.cell_coms())
    types = list(world.cell_types())
    idx = M.recruitment_index_all_responders(
        coms, types, source_type=CR.CT.SOURCE_TYPE,
        responder_types=RESPONDER_TYPES, radius=radius)
    activation = {str(cid): int(types[cid] == CR.ACTIVATED_TYPE) for cid in responder_ids}
    return idx, activation


def _run_seed(*, blocked, seed, steps, kd, radius):
    """Build the ``recruitment_receptor`` composite for one seed, step it via
    ``process_bigraph.Composite``, and sample the all-responders recruitment
    index + per-responder activation every ``SAMPLE`` ticks."""
    core = pb.allocate_core()
    doc = CR.recruitment_receptor(core=core, blocked=blocked, seed=seed, kd=kd)
    comp = pb.Composite({"state": doc}, core=core)
    responder_ids = CR._responder_ids(CR.build_receptor_spec(seed=seed))

    world = comp.state["cpm"]["instance"].world
    mcs_series, idx_series, activation_series = [], [], []
    idx0, act0 = _frame(world, responder_ids, radius)
    mcs_series.append(0); idx_series.append(idx0); activation_series.append(act0)
    for s in range(SAMPLE, steps + 1, SAMPLE):
        comp.run(SAMPLE)
        world = comp.state["cpm"]["instance"].world
        idx, act = _frame(world, responder_ids, radius)
        mcs_series.append(s); idx_series.append(idx); activation_series.append(act)
    return mcs_series, idx_series, activation_series, responder_ids


def _mean_ci(values):
    """Mean + 95% CI across seeds (normal approximation on the seed spread)."""
    n = len(values)
    mean = sum(values) / n
    if n < 2:
        return mean, [mean, mean]
    var = sum((v - mean) ** 2 for v in values) / (n - 1)
    margin = 1.96 * math.sqrt(var / n)
    return mean, [mean - margin, mean + margin]


def run_receptor(*, blocked: bool, seeds: list[int], steps: int, kd: float,
                 out_dir: str | None = None) -> dict:
    """Run ``recruitment_receptor`` for each seed via the process-bigraph
    Composite, aggregate the all-responders recruitment index across seeds,
    and write a summary JSON to
    ``workspace/chemotaxis_data/results/receptor_<condition>.json``.

    Returns the summary dict (also written to disk).
    """
    radius = CR.meta()["recruitment_radius"]
    per_seed = [_run_seed(blocked=blocked, seed=seed, steps=steps, kd=kd, radius=radius)
                for seed in seeds]
    mcs_series = per_seed[0][0]
    responder_ids = per_seed[0][3]
    n_frames = len(mcs_series)

    recruitment_index_mean, recruitment_index_ci = [], []
    for f in range(n_frames):
        vals = [idx_series[f] for _, idx_series, _, _ in per_seed]
        mean, ci = _mean_ci(vals)
        recruitment_index_mean.append(round(mean, 4))
        recruitment_index_ci.append([round(ci[0], 4), round(ci[1], 4)])

    # per-responder activation fraction over time, averaged across seeds
    # (for the spatial figure): {cid: [frac_at_t0, frac_at_t1, ...]}
    activation_by_cell: dict[str, list[float]] = {str(cid): [] for cid in responder_ids}
    for f in range(n_frames):
        for cid in responder_ids:
            key = str(cid)
            frac = sum(act_series[f][key] for _, _, act_series, _ in per_seed) / len(per_seed)
            activation_by_cell[key].append(round(frac, 4))

    final_mean, final_ci = recruitment_index_mean[-1], recruitment_index_ci[-1]

    condition = "blocked" if blocked else "baseline"
    summary = {
        "condition": condition,
        "mcs": mcs_series,
        "recruitment_index_mean": recruitment_index_mean,
        "recruitment_index_ci": recruitment_index_ci,
        "final_mean": final_mean,
        "final_ci": final_ci,
        "activation_by_cell": activation_by_cell,
        "kd": float(kd),
        "seeds": list(seeds),
        "params": {
            "blocked": bool(blocked),
            "steps": int(steps),
            "radius": radius,
            "conc_scale": CR.CONC_SCALE_DEFAULT,
            "hill": CR.HILL_DEFAULT,
            "activate_occupancy": CR.ACTIVATE_OCCUPANCY_DEFAULT,
            "responder_types": list(RESPONDER_TYPES),
            "sample": SAMPLE,
            "cite": "nasser2009",
        },
    }

    out_dir = out_dir or _default_out_dir()
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"receptor_{condition}.json")
    with open(out_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    return summary
