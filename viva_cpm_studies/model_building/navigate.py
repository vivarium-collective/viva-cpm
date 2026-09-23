"""Deterministic NAVIGATE policy for the recruitment mechanism ladder.

Climbs ``mechanisms.MECHANISMS`` (static_lambda -> hill_occupancy ->
adaptive_receptor) from a graded contract: given the currently-installed
mechanism and the graded axes produced by the test harness, pick the next
mechanism to install, or ``None`` once all hard tests pass.
"""

# Maps the failing hard test -> the mechanism that fixes it (the ladder order).
_FIXES = {
    "receptor_gating": "hill_occupancy",
    "recruits_high": "adaptive_receptor",
    "recruits_low": "static_lambda",
}
_ORDER = ("static_lambda", "hill_occupancy", "adaptive_receptor")


def next_mechanism(active, graded_axes):
    """Given the currently-installed mechanism and the graded contract axes
    (each a dict with id/severity/verdict/margin), return the next mechanism to
    install, or None when all hard tests pass. Deterministic policy: take the
    most-negative-margin failing HARD axis; install the mechanism it pulls on;
    if that mechanism is already active (still failing), advance to the next rung
    in _ORDER so the loop never stalls installing what it already has.
    """
    hard_fails = [a for a in graded_axes
                  if a.get("severity") == "hard" and a.get("verdict") == "mismatch"]
    if not hard_fails:
        return None
    worst = min(hard_fails, key=lambda a: a.get("margin", 0.0))
    target = _FIXES.get(worst["id"])
    if target is None:
        return None
    if target == active:
        i = _ORDER.index(active)
        return _ORDER[i + 1] if i + 1 < len(_ORDER) else None
    return target
