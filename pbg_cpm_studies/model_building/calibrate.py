"""Calibration as a program, not a conversation.

The closed-loop methods review flagged the hand grid sweep as the single biggest
token cost of a model-build run (~275k tokens, twice) and diagnosed two fixable
pathologies: (1) not knowing which parameters mattered before searching, and
(2) unpaired random seeds letting simulator noise swamp the signal. This module
retires both, numpy-only:

  * ``elementary_effects`` — a Morris-style sensitivity screen that ranks
    parameters by influence BEFORE any search, so a screen of a handful of runs
    tells you which knobs to sweep and which are empirically irrelevant.
  * ``common_random_numbers`` — wrap a stochastic simulator so every objective
    evaluation reuses the SAME seed set; comparisons become paired and the
    seed-to-seed noise cancels, instead of forcing a redundant re-run.
  * ``refine`` / ``calibrate`` — a coarse deterministic grid over ONLY the
    influential parameters, emitting a ``report_card_verdict``-style margin table
    (which configs clear the acceptance band, and by how much).

The objective is any ``f(params: dict[str, float]) -> float`` — mechanism-agnostic;
the recruitment mechanism ladder is one caller (see the module tests). Escalate to
SALib Sobol indices, a Gaussian-process surrogate, or Bayesian optimization only
past the triggers in the review: a global variance decomposition when a screen is
ambiguous; a surrogate when a single run exceeds ~1 hour; BO when a calibration
needs more than ~50 objective calls after screening. Until then this is enough.
"""
from __future__ import annotations

import itertools
from typing import Callable, Dict, Iterable, Sequence

import numpy as np

Params = Dict[str, float]
Bounds = Dict[str, "tuple[float, float]"]
Objective = Callable[[Params], float]


def _scale(unit: Params, bounds: Bounds) -> Params:
    """Map a point in the unit cube (each coord in [0,1]) to real parameter values."""
    return {n: lo + u * (hi - lo) for n, u in unit.items() for (lo, hi) in [bounds[n]]}


def common_random_numbers(sim: Callable[..., float], seeds: Sequence[int]) -> Objective:
    """Wrap a per-seed simulator into a paired-seed objective.

    ``sim(params, seed)`` returns one noisy scalar; the returned ``objective(params)``
    averages ``sim`` over the FIXED ``seeds`` for every call. Because the same seeds
    are reused at every parameter setting, the simulator's seed-to-seed variance is
    common to all comparisons and cancels — the reason the source run had to sweep
    twice was unpaired seeds, and this is the fix. Deterministic given ``seeds``."""
    seeds = tuple(int(s) for s in seeds)

    def objective(params: Params) -> float:
        return float(np.mean([sim(params, s) for s in seeds]))

    return objective


def elementary_effects(objective: Objective, bounds: Bounds, *,
                       levels: int = 4, trajectories: int = 8, seed: int = 0) -> Dict[str, float]:
    """Morris-style elementary-effects screen. Returns ``{param: mu_star}`` — the
    mean ABSOLUTE elementary effect per parameter, i.e. an influence ranking.

    Cost is ``trajectories * (len(bounds) + 1)`` objective calls — a handful, far
    cheaper than a full grid, and it tells you which knobs a subsequent ``refine``
    should sweep. ``objective`` should already use common random numbers so its own
    noise does not confound the effects. The design RNG is seeded for reproducibility."""
    names = list(bounds)
    if not names:
        return {}
    rng = np.random.default_rng(seed)
    step = 1.0 / (levels - 1)
    grid = np.linspace(0.0, 1.0, levels)
    effects: Dict[str, list] = {n: [] for n in names}
    for _ in range(trajectories):
        point = {n: float(rng.choice(grid)) for n in names}
        f_prev = objective(_scale(point, bounds))
        for n in rng.permutation(np.array(names, dtype=object)):
            n = str(n)
            nxt = dict(point)
            d = step if nxt[n] + step <= 1.0 + 1e-9 else -step
            nxt[n] = min(1.0, max(0.0, nxt[n] + d))
            f = objective(_scale(nxt, bounds))
            effects[n].append(abs((f - f_prev) / d))
            point, f_prev = nxt, f
    return {n: (float(np.mean(v)) if v else 0.0) for n, v in effects.items()}


def rank_influential(mu_star: Dict[str, float], *, top_k: int | None = None,
                     rel_threshold: float = 0.1) -> list[str]:
    """Parameters worth sweeping: sorted by influence, keeping those whose mu_star
    is at least ``rel_threshold`` of the maximum (or the top_k, if given)."""
    if not mu_star:
        return []
    ordered = sorted(mu_star, key=lambda n: mu_star[n], reverse=True)
    if top_k is not None:
        return ordered[:top_k]
    hi = mu_star[ordered[0]] or 1.0
    return [n for n in ordered if mu_star[n] >= rel_threshold * hi]


def _margin(value: float, band: "tuple[float, float]") -> float:
    """Signed distance INTO the acceptance band (>=0 = inside), matching the
    signed-margin convention of test_contract.band."""
    lo, hi = band
    return min(value - lo, hi - value)


def refine(objective: Objective, bounds: Bounds, params: Iterable[str], *,
           band: "tuple[float, float]", levels: int = 5,
           fixed: Params | None = None) -> list[dict]:
    """Coarse grid over ONLY ``params`` (the influential ones), holding the rest at
    ``fixed`` (or each bound's midpoint). Returns a margin table sorted best-first:
    ``[{"params", "value", "margin", "within"}]`` — configs that clear ``band`` have
    ``margin >= 0``. Deterministic. This is the small, legible replacement for the
    hand grid: it only sweeps what the screen said matters."""
    params = list(params)
    mids = {n: 0.5 * (lo + hi) for n, (lo, hi) in bounds.items()}
    base = {**mids, **(fixed or {})}
    axes = [np.linspace(bounds[n][0], bounds[n][1], levels) for n in params]
    table = []
    for combo in itertools.product(*axes) if params else [()]:
        cfg = {**base, **dict(zip(params, (float(x) for x in combo)))}
        v = objective(cfg)
        table.append({"params": cfg, "value": v, "margin": _margin(v, band),
                      "within": _margin(v, band) >= 0.0})
    table.sort(key=lambda r: r["margin"], reverse=True)
    return table


def calibrate(objective: Objective, bounds: Bounds, *, band: "tuple[float, float]",
              top_k: int | None = None, screen_trajectories: int = 8,
              refine_levels: int = 5, seed: int = 0) -> dict:
    """Screen → refine over the influential parameters → report.

    Returns ``{"mu_star", "influential", "best", "within", "table", "n_calls"}``.
    ``best`` is the top-margin config; ``within`` is whether it clears ``band``.
    ``n_calls`` counts objective evaluations, so the token/compute win over a full
    grid is auditable. Escalate (SALib/BO/surrogate) only past the review's triggers."""
    calls = {"n": 0}

    def counted(p: Params) -> float:
        calls["n"] += 1
        return objective(p)

    mu = elementary_effects(counted, bounds, trajectories=screen_trajectories, seed=seed)
    influential = rank_influential(mu, top_k=top_k) or list(bounds)
    table = refine(counted, bounds, influential, band=band, levels=refine_levels)
    best = table[0] if table else None
    return {
        "mu_star": mu,
        "influential": influential,
        "best": best,
        "within": bool(best and best["within"]),
        "table": table,
        "n_calls": calls["n"],
    }
