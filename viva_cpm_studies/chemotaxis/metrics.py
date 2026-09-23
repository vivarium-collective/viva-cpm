"""Recruitment observables for the chemotactic-recruitment slice.

The core readout is **recruitment index** — the fraction of responder cells whose
centre of mass has reached the source region (within ``radius`` of the source
centroid). It rises 0 -> ~1 as responders climb the secreted gradient, and stays
near 0 when the chemotactic response is inhibited. A companion **mean approach**
readout (normalised reduction in responder-to-source distance) is scale-free and
convenient for thresholds.

All functions come in two forms:

* ``*_from_coms(coms, types, ...)`` — pure functions over emitted run data
  (a list of centre-of-mass triples + a parallel list of integer cell types).
  These are what the behavior-contract evaluator calls on recorded frames.
* ``*(world, ...)`` — thin wrappers that pull ``world.cell_coms()`` /
  ``world.cell_types()`` from a live CPM world (for run scripts / notebooks).

Cell type 0 is Medium and is always ignored.
"""
from __future__ import annotations

import math
from typing import Iterable, Sequence

SOURCE_TYPE = 1
RESPONDER_TYPE = 2
DEFAULT_RADIUS = 15.0


def _centroid(points: Sequence[Sequence[float]]) -> tuple[float, float, float]:
    n = len(points)
    if n == 0:
        return (math.nan, math.nan, math.nan)
    sx = sum(p[0] for p in points)
    sy = sum(p[1] for p in points)
    sz = sum(p[2] if len(p) > 2 else 0.0 for p in points)
    return (sx / n, sy / n, sz / n)


def _dist(a: Sequence[float], b: Sequence[float]) -> float:
    return math.sqrt(
        (a[0] - b[0]) ** 2
        + (a[1] - b[1]) ** 2
        + ((a[2] if len(a) > 2 else 0.0) - (b[2] if len(b) > 2 else 0.0)) ** 2
    )


def _responder_type_set(responder_type, responder_types) -> set:
    """Resolve the effective set of responder type ids. ``responder_types``
    (an iterable of type ids), when given, takes precedence over the single
    ``responder_type`` -- letting a receptor-level composite count naive(2)
    AND activated(3) cells as responders without disturbing the single-type
    default every existing chemotaxis study call relies on."""
    if responder_types is not None:
        return {int(t) for t in responder_types}
    return {int(responder_type)}


def _split(coms, types, source_type, responder_type_set):
    """Return (source_coms, responder_coms) skipping Medium (type 0).

    ``coms``/``types`` may be indexed by cell id (index 0 == Medium, as
    ``world.cell_coms()`` returns) or be a compact parallel pair — both work as
    long as they are the same length.
    """
    src, resp = [], []
    for com, t in zip(coms, types):
        ti = int(t)
        if ti == source_type:
            src.append(com)
        elif ti in responder_type_set:
            resp.append(com)
    return src, resp


def recruitment_index_from_coms(
    coms, types, *, source_type=SOURCE_TYPE, responder_type=RESPONDER_TYPE,
    responder_types=None, radius=DEFAULT_RADIUS,
) -> float:
    """Fraction of responder cells within ``radius`` of the source centroid.

    Counts only ``responder_type`` by default (backward compatible with every
    existing chemotaxis study). Pass ``responder_types`` — an iterable of type
    ids, e.g. ``{NAIVE_TYPE, ACTIVATED_TYPE}`` — to count several responder
    types at once, which a receptor-level composite needs: gating recruitment
    on "reached activated (type 3)" alone is circular (activation and
    recruitment would share the same underlying occupancy signal), so
    receptor-level recruitment counts ALL responders regardless of fate.
    """
    types_set = _responder_type_set(responder_type, responder_types)
    src, resp = _split(coms, types, source_type, types_set)
    if not resp or not src:
        return 0.0
    c = _centroid(src)
    within = sum(1 for r in resp if _dist(r, c) <= radius)
    return within / len(resp)


def recruitment_index_all_responders(
    coms, types, *, source_type=SOURCE_TYPE,
    responder_types=(RESPONDER_TYPE, RESPONDER_TYPE + 1), radius=DEFAULT_RADIUS,
) -> float:
    """Thin wrapper: recruitment index counting ALL responder cell types
    (naive + activated), not activated-only. Defaults to ``{2, 3}`` — the
    naive/activated pair used by ``recruitment_receptor``."""
    return recruitment_index_from_coms(
        coms, types, source_type=source_type, responder_types=responder_types,
        radius=radius)


def mean_distance_to_source_from_coms(
    coms, types, *, source_type=SOURCE_TYPE, responder_type=RESPONDER_TYPE,
    responder_types=None,
) -> float:
    """Mean responder centre-of-mass distance to the source centroid."""
    types_set = _responder_type_set(responder_type, responder_types)
    src, resp = _split(coms, types, source_type, types_set)
    if not resp or not src:
        return math.nan
    c = _centroid(src)
    return sum(_dist(r, c) for r in resp) / len(resp)


def mean_approach(first_distance: float, last_distance: float) -> float:
    """Scale-free recruitment: fraction of the initial gap that was closed.

    1.0 == responders reached the source; 0.0 == no net approach; negative ==
    responders moved away.
    """
    if not first_distance or math.isnan(first_distance):
        return 0.0
    return (first_distance - last_distance) / first_distance


# --- live-world convenience wrappers -------------------------------------------

def recruitment_index(world, *, source_type=SOURCE_TYPE,
                      responder_type=RESPONDER_TYPE, responder_types=None,
                      radius=DEFAULT_RADIUS) -> float:
    return recruitment_index_from_coms(
        list(world.cell_coms()), list(world.cell_types()),
        source_type=source_type, responder_type=responder_type,
        responder_types=responder_types, radius=radius)


def mean_distance_to_source(world, *, source_type=SOURCE_TYPE,
                            responder_type=RESPONDER_TYPE) -> float:
    return mean_distance_to_source_from_coms(
        list(world.cell_coms()), list(world.cell_types()),
        source_type=source_type, responder_type=responder_type)
