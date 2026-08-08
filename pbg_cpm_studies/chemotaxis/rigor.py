"""Multi-seed replication for the chemotactic-recruitment slice.

Runs each condition across several independent seeds and aggregates the final
recruitment index (mean +/- sd), plus the effect size separating the baseline
from each negative condition. This turns the single-seed findings into replicated
ones — the evidence the rigor scorecard asks for.

Run (from the worktree root, viva-cpm venv)::

    PYTHONPATH=. python -m pbg_cpm_studies.chemotaxis.rigor
"""
from __future__ import annotations

import json
import math
import os

from .run import run_condition, CONDITIONS

SEEDS = [17, 29, 43, 61, 89]


def _mean_sd(xs):
    n = len(xs)
    m = sum(xs) / n
    var = sum((x - m) ** 2 for x in xs) / n if n > 1 else 0.0
    return m, math.sqrt(var)


def _cohens_d(a, b):
    """Effect size between two samples (pooled-sd). Large sd-free separations
    (a moves, b is pinned at 0) return a big finite number, capped for display."""
    ma, sa = _mean_sd(a)
    mb, sb = _mean_sd(b)
    sp = math.sqrt((sa ** 2 + sb ** 2) / 2) or 1e-6
    d = (ma - mb) / sp
    return round(min(d, 999.0), 2)


def main(out_dir=None):
    here = os.path.dirname(os.path.abspath(__file__))
    ws = os.path.abspath(os.path.join(here, "..", "..", "workspace"))
    out_dir = out_dir or os.path.join(ws, "chemotaxis_data", "results")
    os.makedirs(out_dir, exist_ok=True)

    per_cond = {}
    print(f"{'condition':26s} {'n':>2s} {'mean final':>11s} {'sd':>6s}  seeds")
    for name, kw in CONDITIONS.items():
        finals = []
        for seed in SEEDS:
            res = run_condition(name, seed=seed, **kw)
            finals.append(res["final_recruitment_index"])
        m, sd = _mean_sd(finals)
        per_cond[name] = {"seeds": SEEDS, "final_recruitment_index": finals,
                          "mean": round(m, 4), "sd": round(sd, 4)}
        print(f"{name:26s} {len(SEEDS):2d} {m:11.3f} {sd:6.3f}  {[round(x,2) for x in finals]}")

    base = per_cond["recruitment-baseline"]["final_recruitment_index"]
    effects = {
        "baseline_vs_inhibited":
            _cohens_d(base, per_cond["recruitment-inhibited"]["final_recruitment_index"]),
        "baseline_vs_adversarial":
            _cohens_d(base, per_cond["recruitment-adversarial"]["final_recruitment_index"]),
    }
    summary = {"n_seeds": len(SEEDS), "seeds": SEEDS, "per_condition": per_cond,
               "effect_size_cohens_d": effects}
    with open(os.path.join(out_dir, "replication_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print("\neffect size (Cohen's d):", effects)
    print(f"wrote replication_summary.json -> {out_dir}")
    return summary


if __name__ == "__main__":
    main()
