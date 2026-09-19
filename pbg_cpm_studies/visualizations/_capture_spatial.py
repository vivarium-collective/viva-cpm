"""Capture CPM 2D spatial-state frames for every influenza-sego2022 study.

For each study we run its REAL ``run.run_*`` driver (unmodified) and record the
per-site cell-owner lattice (``world.snapshot()``) plus the per-cell type table
(``world.cell_types()``) BEFORE every ``world.step`` (which reflects the prior
update's fully-resolved transitions), then one final post-run snapshot. Each
captured owner grid is mapped to a per-site cell-TYPE grid (MEDIUM=0, H=1, I=2,
D=3, M=4, K=5, E=6) and mode-coarsened so no side exceeds ~60. The result is
baked to ``_influenza_spatial.py`` as the ``INFLUENZA_SPATIAL`` JSON blob that
the ``InfluenzaSpatial*`` dashboard cards render from (no live simulation at
render time, identical in the live and published read-only dashboards).

We do NOT touch any driver's logic: the recorder wraps the ``cpm_core.World``
the driver builds (via a temporary monkeypatch of ``build.world_from_spec``) and
only observes it. Deterministic seeds match each study's reported config.

Run (from the worktree root, viva-cpm venv)::

    .venv/bin/python -m pbg_cpm_studies.visualizations._capture_spatial
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from pbg_cpm_studies.influenza import build, run

TARGET_MAX_SIDE = 60   # coarsen so max(nx, ny) <= this
N_FRAMES = 9           # evenly-spaced frames per study (incl. first + last)


# ── recorder: wrap the world the driver builds, observe it, never mutate ─────
class _RecWorld:
    """Transparent proxy over a ``cpm_core.World``: forwards every attribute to
    the real world, and on each ``.step()`` records the PRE-step lattice
    (owner grid + per-cell type table) — which reflects the previous update's
    resolved transitions."""

    def __init__(self, world):
        object.__setattr__(self, "_w", world)
        object.__setattr__(self, "_frames", [])  # list[(owner_flat, types_list)]

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, "_w"), name)

    def _snap(self):
        w = object.__getattribute__(self, "_w")
        self._frames.append((list(w.snapshot()), list(w.cell_types())))

    def step(self, *a, **k):
        self._snap()
        return object.__getattribute__(self, "_w").step(*a, **k)


_CREATED: list[_RecWorld] = []
_ORIG_WFS = build.world_from_spec


def _patched_world_from_spec(*a, **k):
    w = _ORIG_WFS(*a, **k)
    rec = _RecWorld(w)
    _CREATED.append(rec)
    return rec


def _run_capture(driver, /, **kwargs):
    """Run a driver under the recorder; return the PRIMARY world's frames
    (owner_flat, types_list) as a list, plus its dims (nx, ny). The primary
    world is the one that was stepped the most (drivers may build auxiliary
    close-contact/control worlds; the main scene has the most steps)."""
    _CREATED.clear()
    build.world_from_spec = _patched_world_from_spec
    try:
        driver(**kwargs)
    finally:
        build.world_from_spec = _ORIG_WFS
    primary = max(_CREATED, key=lambda r: len(r._frames))
    primary._snap()  # final post-run state (after the last step's transitions)
    frames = list(primary._frames)
    nx, ny, _nz = object.__getattribute__(primary, "_w").dims()
    return frames, int(nx), int(ny)


# ── coarsening (categorical mode-pool; keeps lesions/immune clusters) ────────
def _coarsen_mode(grid2d: np.ndarray, k: int, n_labels: int) -> np.ndarray:
    if k <= 1:
        return grid2d
    ny, nx = grid2d.shape
    ny2, nx2 = ny // k, nx // k
    g = grid2d[: ny2 * k, : nx2 * k].reshape(ny2, k, nx2, k)
    out = np.zeros((ny2, nx2), dtype=np.int64)
    for i in range(ny2):
        for j in range(nx2):
            block = g[i, :, j, :].ravel()
            out[i, j] = np.bincount(block, minlength=n_labels).argmax()
    return out


def _subsample(frames: list, n: int) -> list[int]:
    """Indices of ~n evenly-spaced frames, always including first and last."""
    total = len(frames)
    if total <= n:
        return list(range(total))
    return sorted({round(i * (total - 1) / (n - 1)) for i in range(n)})


def _type_frames(frames, nx, ny, mcs_of):
    """Map each captured owner grid -> per-site cell-TYPE grid, coarsened."""
    k = max(1, math.ceil(max(nx, ny) / TARGET_MAX_SIDE))
    idxs = _subsample(frames, N_FRAMES)
    out = []
    for fi in idxs:
        owners_flat, types_list = frames[fi]
        owners = np.asarray(owners_flat, dtype=np.int64).reshape(ny, nx)
        tg = np.asarray(types_list, dtype=np.int64)[owners]      # per-site type
        tg = _coarsen_mode(tg, k, 7)
        out.append({"mcs": int(mcs_of(fi)), "grid": tg.astype(int).ravel().tolist()})
    return out, nx // k, ny // k, k


def _owner_frames(frames, nx, ny, mcs_of):
    """Coarsened OWNER grid (cell-label mosaic) — used for the all-H substrate
    relaxation, where a type grid would be a single flat colour."""
    k = max(1, math.ceil(max(nx, ny) / TARGET_MAX_SIDE))
    n_labels = max((max(f[0]) for f in frames), default=0) + 1
    idxs = _subsample(frames, N_FRAMES)
    out = []
    for fi in idxs:
        owners_flat, _types = frames[fi]
        owners = np.asarray(owners_flat, dtype=np.int64).reshape(ny, nx)
        og = _coarsen_mode(owners, k, n_labels)
        out.append({"mcs": int(mcs_of(fi)), "grid": og.astype(int).ravel().tolist()})
    return out, nx // k, ny // k, k


# ── per-study capture specs ──────────────────────────────────────────────────
def capture_all() -> dict:
    spatial: dict = {}

    # ---- epithelial-sheet-baseline: substrate geometry relaxation (owner) ----
    # No transition driver — build the confluent sheet and relax it under the
    # Potts temperature, recording the tile mosaic (owner grid) shifting.
    from pbg_cpm_studies.influenza import sheet as _sheet
    spec = _sheet.build_sheet_spec(0.3, seed=17)
    world = _RecWorld(_ORIG_WFS(spec, finalize=True))
    relax_per = 25
    n_relax = 8
    for _ in range(n_relax):
        world.step(relax_per)
    world._snap()
    frames = list(world._frames)
    nx, ny, _nz = object.__getattribute__(world, "_w").dims()
    of, gnx, gny, k = _owner_frames(frames, int(nx), int(ny),
                                    mcs_of=lambda fi: fi * relax_per)
    spatial["epithelial-sheet-baseline"] = {
        "kind": "owner", "nx": gnx, "ny": gny, "coarsen": k,
        "seed": 17, "frames": of,
    }
    print(f"sheet: {len(of)} frames {gnx}x{gny} (k={k})")

    # ---- the transition/immune studies: per-site cell-TYPE mosaic ----
    STUDIES = [
        ("virus-field-infection", run.run_virus_infection,
         dict(patch_mm=0.3, steps=60, seed=17), 10),
        ("ifn-resistance", run.run_virus_infection_with_ifn,
         dict(patch_mm=0.3, steps=60, seed=17), 10),
        ("epithelial-fate", run.run_epithelial_fate,
         dict(patch_mm=0.3, steps=200, seed=17), 10),
        ("macrophage-response", run.run_macrophage_response,
         dict(steps=40, seed=17), 10),
        ("signaling-fields", run.run_macrophage_signaling,
         dict(steps=20, seed=17), 10),
        ("cytotoxic-killing", run.run_cytotoxic_response,
         dict(steps=60, seed=17), 10),
        ("global-coupling", run.run_global_coupling,
         dict(side=30, steps=20, seed=3), 1),
        # repro studies share one full-model scene (reduced-scale, seed 0)
        ("repro-full-model", run.run_full_model,
         dict(cells_per_side=15, steps=20, seed=0, init_infection_frac=0.05), 7),
    ]
    for slug, driver, kwargs, mcs_per in STUDIES:
        frames, nx, ny = _run_capture(driver, **kwargs)
        tf, gnx, gny, k = _type_frames(frames, nx, ny,
                                       mcs_of=lambda fi: fi * mcs_per)
        spatial[slug] = {
            "kind": "type", "nx": gnx, "ny": gny, "coarsen": k,
            "seed": int(kwargs.get("seed", 17)), "mcs_per": mcs_per, "frames": tf,
        }
        print(f"{slug}: {len(tf)} frames {gnx}x{gny} (k={k}) "
              f"from {len(frames)} captured")
    return spatial


def main() -> int:
    spatial = capture_all()
    blob = json.dumps(spatial, separators=(",", ":"))
    out = Path(__file__).resolve().parent / "_influenza_spatial.py"
    header = (
        '"""Baked CPM 2D spatial-state frames for the influenza-sego2022 '
        'studies.\n\n'
        "REAL engine output captured by ``_capture_spatial.py`` (run its module "
        "to\nregenerate). Each study maps to a per-site cell-TYPE grid "
        "(MEDIUM=0, H=1,\nI=2, D=3, M=4, K=5, E=6), mode-coarsened; the sheet "
        "substrate stores the\nowner (cell-label) mosaic instead. Rendered by "
        "``influenza_studies.py``'s\n``InfluenzaSpatial*`` cards.\n"
        '"""\n'
    )
    out.write_text(header + "INFLUENZA_SPATIAL = " + repr(blob) + "\n",
                   encoding="utf-8")
    print(f"\nwrote {out} — {len(blob):,} bytes JSON")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
