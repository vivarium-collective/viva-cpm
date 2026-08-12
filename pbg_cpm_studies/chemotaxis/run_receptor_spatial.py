"""Spatial-frames recorder for the receptor-recruitment composite.

The scalar summaries written by :mod:`pbg_cpm_studies.chemotaxis.run_receptor`
(recruitment index + per-cell activation fraction over time) throw away *where*
the responders are and *what* the chemokine field looks like. Those two things
are exactly what makes the chemotaxis mechanism legible, so this module records,
for ONE representative seed per condition, a light per-frame spatial artifact:

* every responder cell's centre of mass ``(x, y)``, its CPM ``type`` (naive 2 /
  activated 3), its living ``volume``, the local chemokine concentration it
  senses (``field_mean_at_cell`` * ``conc_scale``) and the resulting receptor
  **occupancy** ``theta`` (recomputed here from the SAME Hill relation the
  per-cell ``ReceptorSubcell`` uses, so the recorded ``theta`` is the quantity
  that actually gated activation);
* the source cell(s) centre of mass + bounding-box slab extent (the source is
  NOT anchored -- it random-walks across the lattice, and the chemokine blob
  tracks it, so both are recorded per frame);
* a block-mean-**downsampled** chemokine field grid from ``world.field_conc(0)``
  (native 80x40 -> 40x20), the diffusing gradient the responders climb.

Built and stepped through a ``process_bigraph.Composite`` (NOT raw
``load_world``) so the per-responder subcells fire and cells actually activate
-- mirrors :func:`run_receptor._run_seed`.

Run (from the worktree root, viva-cpm venv)::

    PYTHONPATH=. python -c "from pbg_cpm_studies.chemotaxis.run_receptor_spatial \\
        import record_spatial; record_spatial()"

Writes ``workspace/chemotaxis_data/results/receptor_spatial_baseline.json`` and
``receptor_spatial_blocked.json`` (seed 17, every 25 ticks, 21 frames). Those
JSONs are baked figure data (committed, like the existing summaries) and are
inlined into ``visualizations/receptor_spatial_studies.py``.
"""
from __future__ import annotations

import json
import os

import process_bigraph as pb

from ..composites import chemotaxis as CT
from ..composites import chemotaxis_receptor as CR

SAMPLE = 25          # record every SAMPLE composite ticks (matches run_receptor)
FIELD_DS = 2         # field downsample factor: native 80x40 -> 40x20
REPR_SEED = 17       # the one representative seed we record spatially


def _default_out_dir():
    here = os.path.dirname(os.path.abspath(__file__))
    ws = os.path.abspath(os.path.join(here, "..", "..", "workspace"))
    return os.path.join(ws, "chemotaxis_data", "results")


def _occupancy(field_mean, *, kd, hill, conc_scale):
    """Fractional receptor occupancy from the local field mean -- identical to
    ``cpm.subcellular.receptor.ReceptorSubcell.occupancy`` so the recorded
    theta is the exact quantity that gated the naive->activated switch."""
    c = max(0.0, float(field_mean)) * conc_scale
    if c <= 0.0:
        return 0.0, 0.0
    ch = c ** hill
    return ch / (kd ** hill + ch), c


def _downsample_field(fc, nx, ny, factor):
    """Block-mean downsample the flat field (owner index x + y*nx) into a
    rows(y) x cols(x) grid, plus the lattice-unit centre coord of each block."""
    dnx, dny = nx // factor, ny // factor
    grid = []
    for j in range(dny):
        row = []
        for i in range(dnx):
            acc = 0.0
            for dy in range(factor):
                base = (j * factor + dy) * nx + i * factor
                for dx in range(factor):
                    acc += fc[base + dx]
            row.append(round(acc / (factor * factor), 1))
        grid.append(row)
    xs = [i * factor + (factor - 1) / 2.0 for i in range(dnx)]
    ys = [j * factor + (factor - 1) / 2.0 for j in range(dny)]
    return grid, xs, ys


def _source_extent(world, nx, ny, source_ids):
    """Bounding box [x0, y0, x1, y1] of the sites owned by any source cell."""
    lab = world.snapshot()
    sids = set(source_ids)
    x0 = y0 = 10 ** 9
    x1 = y1 = -1
    for idx, owner in enumerate(lab):
        if owner in sids:
            x = idx % nx
            y = (idx // nx) % ny
            if x < x0: x0 = x
            if y < y0: y0 = y
            if x > x1: x1 = x
            if y > y1: y1 = y
    if x1 < 0:
        return None
    return [x0, y0, x1, y1]


def _frame(world, responder_ids, nx, ny, *, kd, hill, conc_scale, activate_occupancy):
    types = list(world.cell_types())
    vols = list(world.cell_volumes())
    coms = list(world.cell_coms())
    source_ids = [cid for cid in range(len(types)) if types[cid] == CT.SOURCE_TYPE]

    cells = []
    for cid in responder_ids:
        vol = int(vols[cid])
        if vol <= 0:
            # cell has been squeezed out of the lattice: it has no meaningful
            # position this frame, so record it as dead rather than plot it at
            # the (0,0,0) sentinel the engine returns for an empty cell.
            cells.append({"id": cid, "alive": False})
            continue
        fmac = world.field_mean_at_cell(0, cid)
        theta, conc = _occupancy(fmac, kd=kd, hill=hill, conc_scale=conc_scale)
        x, y = coms[cid][0], coms[cid][1]
        cells.append({
            "id": cid, "alive": True, "x": round(x, 2), "y": round(y, 2),
            "type": int(types[cid]), "vol": vol,
            "conc": round(conc, 4), "theta": round(theta, 4),
            "activated": bool(theta >= activate_occupancy),
        })

    # source centre of mass (average over source cells; there is one slab here)
    src_coms = [coms[cid] for cid in source_ids]
    if src_coms:
        sx = round(sum(c[0] for c in src_coms) / len(src_coms), 2)
        sy = round(sum(c[1] for c in src_coms) / len(src_coms), 2)
    else:
        sx = sy = None
    extent = _source_extent(world, nx, ny, source_ids)

    grid, _, _ = _downsample_field(world.field_conc(0), nx, ny, FIELD_DS)
    return {
        "cells": cells,
        "source_com": [sx, sy],
        "source_extent": extent,
        "field": grid,
    }


def record_spatial_condition(*, blocked, seed=REPR_SEED, steps=500, kd=2.9,
                             out_dir=None):
    """Build + step ``recruitment_receptor`` for one seed and record a spatial
    frame every ``SAMPLE`` ticks. Writes ``receptor_spatial_<condition>.json``
    and returns the artifact dict."""
    kd = float(kd)
    hill = CR.HILL_DEFAULT
    conc_scale = CR.CONC_SCALE_DEFAULT
    activate = CR.ACTIVATE_OCCUPANCY_DEFAULT
    radius = CR.meta()["recruitment_radius"]

    core = pb.allocate_core()
    doc = CR.recruitment_receptor(core=core, blocked=blocked, seed=seed, kd=kd)
    comp = pb.Composite({"state": doc}, core=core)
    responder_ids = CR._responder_ids(CR.build_receptor_spec(seed=seed))

    world = comp.state["cpm"]["instance"].world
    nx, ny, _ = world.dims()

    def snap(mcs):
        f = _frame(world, responder_ids, nx, ny, kd=kd, hill=hill,
                   conc_scale=conc_scale, activate_occupancy=activate)
        f["mcs"] = mcs
        return f

    frames = [snap(0)]
    for s in range(SAMPLE, steps + 1, SAMPLE):
        comp.run(SAMPLE)
        world = comp.state["cpm"]["instance"].world
        frames.append(snap(s))

    # downsampled-grid axis coords (constant across frames)
    _, field_x, field_y = _downsample_field(world.field_conc(0), nx, ny, FIELD_DS)

    condition = "blocked" if blocked else "baseline"
    artifact = {
        "condition": condition,
        "seed": int(seed),
        "kd": kd, "hill": hill, "conc_scale": conc_scale,
        "activate_occupancy": activate, "recruitment_radius": radius,
        "lattice_dims": [nx, ny],
        "field_downsample": FIELD_DS,
        "field_x": [round(v, 2) for v in field_x],
        "field_y": [round(v, 2) for v in field_y],
        "responder_ids": list(responder_ids),
        "naive_type": CR.NAIVE_TYPE,
        "activated_type": CR.ACTIVATED_TYPE,
        "source_type": CT.SOURCE_TYPE,
        "sample": SAMPLE,
        "cite": "nasser2009",
        "frames": frames,
    }
    out_dir = out_dir or _default_out_dir()
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"receptor_spatial_{condition}.json")
    with open(out_path, "w") as fh:
        json.dump(artifact, fh, indent=1)
    return artifact


def record_spatial(*, seed=REPR_SEED, steps=500, kd=2.9, out_dir=None):
    """Record BOTH the baseline and blocked spatial artifacts (seed 17)."""
    base = record_spatial_condition(blocked=False, seed=seed, steps=steps, kd=kd, out_dir=out_dir)
    blk = record_spatial_condition(blocked=True, seed=seed, steps=steps, kd=kd, out_dir=out_dir)
    return {"baseline": base, "blocked": blk}


if __name__ == "__main__":
    record_spatial()
