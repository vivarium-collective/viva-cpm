"""Scale-agnostic acceptance-band comparison harness (Increment 9, Task 9.0).

Compares a model time-series against the Increment-0 digitized acceptance
bands (`targets/{fig3b,fig5,fig7}.json`, loaded via
:func:`viva_cpm_studies.influenza.targets.load_target`): value = ODE
reference, lo/hi = 50-replica spatial-model spread. This module does not
know or care about the scale (small crux-driver runs vs. paper-scale
Increment-9 runs) -- it only interpolates a `[(t_days, value)]` series to
each target observation's `t_days` and checks containment.
"""
from __future__ import annotations

from . import targets

# fig3b (dossier notes, Sec.3.1) flags Type-I IFN + downstream signals/immune
# cells as the figure's most notable spatial-vs-ODE spread; those observables
# carry a per-entry "soft": true in the target JSON. `soft=True` widens each
# target's [lo, hi] band about its `value` by this factor before checking
# containment -- documented, not tuned per-call.
SOFT_WIDEN_FACTOR = 1.5


def load(fig: str) -> dict:
    """Load the digitized acceptance-band target for figure ``fig``
    (e.g. ``"fig3b"``, ``"fig5"``, ``"fig7"``); thin wrapper over
    :func:`targets.load_target`."""
    return targets.load_target(fig)


def _interp(series: list[tuple[float, float]], t: float) -> float:
    """Linearly interpolate ``series`` (a list of ``(t_days, value)``,
    assumed sorted ascending by ``t_days``) at time ``t``. Clamps to the
    series' first/last value outside its time range (flat extrapolation)."""
    ts = [p[0] for p in series]
    vs = [p[1] for p in series]
    if t <= ts[0]:
        return vs[0]
    if t >= ts[-1]:
        return vs[-1]
    for i in range(1, len(ts)):
        if ts[i] >= t:
            t0, t1 = ts[i - 1], ts[i]
            v0, v1 = vs[i - 1], vs[i]
            if t1 == t0:
                return v1
            frac = (t - t0) / (t1 - t0)
            return v0 + frac * (v1 - v0)
    return vs[-1]


def series_in_band(series: list[tuple[float, float]], target_obs: list[dict],
                    *, soft: bool = False) -> dict:
    """Interpolate the model ``series`` (``[(t_days, value)]``) to each
    ``target_obs`` entry's ``t_days`` and check whether it falls in that
    entry's ``[lo, hi]`` band.

    ``soft=True`` widens each entry's band about its ``value`` by
    :data:`SOFT_WIDEN_FACTOR` before checking (``new_lo = value - (value -
    lo) * factor``, ``new_hi = value + (hi - value) * factor``) -- for
    observables the source itself flags as high-spread (fig3b's Type-I IFN
    + downstream immune cells; fig5/fig7's spatial-vs-ODE divergence
    scenarios).

    Returns ``{"in_band": bool, "n_checked": int, "n_in": int,
    "worst_miss": float}``; ``worst_miss`` is the largest relative distance
    of any out-of-band point outside its band (relative to the band width,
    or to the target value when the band has zero width), 0.0 when every
    checked point is in-band.
    """
    n_checked = len(target_obs)
    if not series or n_checked == 0:
        return {"in_band": n_checked == 0, "n_checked": n_checked, "n_in": 0,
                "worst_miss": 0.0 if n_checked == 0 else 1.0}

    series_sorted = sorted(series, key=lambda p: p[0])
    n_in = 0
    worst_miss = 0.0
    for obs in target_obs:
        t = obs["t_days"]
        value = obs["value"]
        lo = obs["lo"]
        hi = obs["hi"]
        if soft:
            lo = value - (value - lo) * SOFT_WIDEN_FACTOR
            hi = value + (hi - value) * SOFT_WIDEN_FACTOR

        v = _interp(series_sorted, t)
        if lo <= v <= hi:
            n_in += 1
        else:
            miss_abs = max(lo - v, 0.0) + max(v - hi, 0.0)
            width = hi - lo
            if width > 0:
                rel = miss_abs / width
            elif value:
                rel = miss_abs / abs(value)
            else:
                rel = miss_abs
            worst_miss = max(worst_miss, rel)

    return {"in_band": n_in == n_checked, "n_checked": n_checked, "n_in": n_in,
            "worst_miss": worst_miss}


def _get_path(d: dict, key_path):
    """Navigate ``d`` by ``key_path`` (a dotted string ``"a.b.c"`` or an
    iterable of keys) and return the value found there."""
    keys = key_path.split(".") if isinstance(key_path, str) else list(key_path)
    cur = d
    for k in keys:
        cur = cur[k]
    return cur


def aggregate_replicas(runs: list[dict], key_path) -> list[tuple]:
    """Given ``runs`` (one dict per replica) and a ``key_path`` locating a
    ``[(t_days, value)]`` series within each replica's dict, return the
    ensemble MEAN series across replicas (elementwise mean of ``value`` at
    each index; ``t_days`` is taken from the first replica, all replicas are
    assumed to share the same time grid)."""
    if not runs:
        return []
    series_list = [_get_path(r, key_path) for r in runs]
    n = len(series_list[0])
    out = []
    for i in range(n):
        t = series_list[0][i][0]
        mean_v = sum(s[i][1] for s in series_list) / len(series_list)
        out.append((t, mean_v))
    return out


def evaluate_study(ensemble: dict, fig: str) -> dict:
    """Evaluate an ``ensemble`` (``{observable_name: [(t_days, value)]}``,
    typically the output of :func:`aggregate_replicas` per observable)
    against figure ``fig``'s acceptance bands.

    Honors each observable's ``soft`` flag as recorded in the target JSON
    (per Increment-0 digitization notes; e.g. fig3b's Type-I IFN + downstream
    immune observables). Observables absent from ``ensemble`` are skipped
    (not counted against ``passed``).

    Returns ``{"figure": fig, "observables": {name: series_in_band(...)},
    "passed": bool}`` -- ``passed`` is True iff every evaluated observable's
    ``in_band`` is True (vacuously True if no observable overlaps).
    """
    target = load(fig)
    observables = target["observables"]
    results = {}
    passed = True
    for name, obs_list in observables.items():
        if name not in ensemble:
            continue
        soft = any(o.get("soft", False) for o in obs_list)
        verdict = series_in_band(ensemble[name], obs_list, soft=soft)
        results[name] = verdict
        if not verdict["in_band"]:
            passed = False
    return {"figure": fig, "observables": results, "passed": passed}
