"""Capture CPM 2D spatial-state frames for every influenza-sego2022 study.

For each study we run its REAL ``run.run_*`` driver (unmodified) and record the
per-site cell-owner lattice (``world.snapshot()``) plus the per-cell type table
(``world.cell_types()``) BEFORE every ``world.step`` (which reflects the prior
update's fully-resolved transitions), then one final post-run snapshot. Each
captured owner grid is mapped to a per-site cell-TYPE grid (MEDIUM=0, H=1, I=2,
D=3, M=4, K=5, E=6, dormant-reserve=7).

Unlike the first cut, these captures are **paper-realistic**:

* **Full resolution** — the type-grid is stored at the REAL lattice size, with
  NO coarsening, so each five-by-five epithelial cell reads as a dense tile
  (like Sego 2022 Fig 2/3), not a blocky down-sampled mosaic.
* **Full scale + full duration** — the capstone ``run_full_model`` scene is
  captured at ``cells_per_side=35`` (a 35×35 = 1225-cell / 0.3 mm epithelial
  patch, the paper's shipped scale) for the WHOLE ~3.5-day run, so the infection
  sweep and the ROS-driven mass death (green H → orange I → grey D) are visible.
* **Compressed frames** — each frame's ``uint8`` type-grid is zlib-compressed
  and base64-encoded (a confluent sheet compresses ~40×), so a full-res
  multi-frame, multi-study bake stays small. Decode in the render with::

      import zlib, base64, numpy as np
      np.frombuffer(zlib.decompress(base64.b64decode(s)), np.uint8).reshape(ny, nx)

The result is baked to ``_influenza_spatial.py`` as the ``INFLUENZA_SPATIAL``
JSON blob that the ``InfluenzaSpatial*`` dashboard cards render from (no live
simulation at render time, identical in the live and published read-only
dashboards).

We do NOT touch any driver's logic: the recorder wraps the ``cpm_core.World``
the driver builds (via a temporary monkeypatch of ``build.world_from_spec``) and
only observes it. Deterministic seeds match each study's reported config.

Run (from the worktree root, viva-cpm venv)::

    .venv/bin/python -m viva_cpm_studies.visualizations._capture_spatial

This is HEAVY: the full-scale full-duration capstone capture takes ~10 min on a
laptop. That is expected — the realism is the point.
"""
from __future__ import annotations

import base64
import json
import zlib
from pathlib import Path

import numpy as np

from viva_cpm_studies.influenza import build, run

N_FRAMES = 18          # evenly-spaced frames per study (incl. first + last)


def _encode(grid2d: np.ndarray) -> str:
    """zlib+base64 a small-int lattice (row-major); decoded in the render."""
    return base64.b64encode(zlib.compress(grid2d.tobytes(), 9)).decode("ascii")


# ── recorder: wrap the world the driver builds, observe it, never mutate ─────
class _RecWorld:
    """Transparent proxy over a ``cpm_core.World``: forwards every attribute to
    the real world, and — at a fixed MCS ``stride`` — records the PRE-step
    lattice (owner grid + per-cell type table), which reflects the previous
    update's resolved transitions. Frames are compressed on capture so a
    full-resolution, full-duration run stays memory-light."""

    def __init__(self, world, stride, *, owner=False):
        object.__setattr__(self, "_w", world)
        object.__setattr__(self, "_frames", [])     # list[{mcs, grid}]
        object.__setattr__(self, "_mcs", 0)          # cumulative MCS stepped
        object.__setattr__(self, "_next", 0)         # next MCS to snapshot at
        object.__setattr__(self, "_stride", int(stride))
        object.__setattr__(self, "_owner", bool(owner))
        object.__setattr__(self, "_nx", None)
        object.__setattr__(self, "_ny", None)
        object.__setattr__(self, "_n_epi", 0)    # epithelial cell count (first snap)
        object.__setattr__(self, "_n_imm", 0)    # immune cell count (peak over run)

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, "_w"), name)

    def _dims(self):
        if self._nx is None:
            nx, ny, _nz = object.__getattribute__(self, "_w").dims()
            object.__setattr__(self, "_nx", int(nx))
            object.__setattr__(self, "_ny", int(ny))
        return self._nx, self._ny

    def _snap(self):
        w = object.__getattribute__(self, "_w")
        nx, ny = self._dims()
        owners = np.asarray(w.snapshot(), dtype=np.int64).reshape(ny, nx)
        # Cell counts: epithelial (types 1-3 = H/I/D) fixed at the intact sheet
        # (first snap); immune (types 4-6 = M/K/E) tracked at its peak (they
        # recruit/move over the run). cell_types() is per-cell (index 0 = medium).
        ct = np.asarray(w.cell_types(), dtype=np.int64)[1:]
        if not self._frames:
            object.__setattr__(self, "_n_epi", int(np.isin(ct, (1, 2, 3)).sum()))
        object.__setattr__(self, "_n_imm",
                           max(self._n_imm, int(np.isin(ct, (4, 5, 6)).sum())))
        if self._owner:
            # cell-label mosaic (substrate relaxation) — labels can exceed 255
            grid = owners.astype(np.uint16)
            self._frames.append({"mcs": int(self._mcs), "grid": _encode(grid)})
        else:
            types_arr = np.asarray(w.cell_types(), dtype=np.int64)
            grid = types_arr[owners].astype(np.uint8)   # per-site cell TYPE
            # Outline every individual cell (adjacent same-TYPE cells share a
            # colour but differ in owner): paint the cell sites whose right or
            # down neighbour has a different owner with the EDGE sentinel (8),
            # giving thin 1px tessellation lines. Baked INTO the type grid (not a
            # separate overlay) so the card embeds one heatmap, not two -> no
            # size increase even at full-sheet density. Only cell sites (owner
            # != 0) are marked, so the outline stays inside cells, not the medium.
            edges = np.zeros_like(owners, dtype=bool)
            edges[:, :-1] |= owners[:, :-1] != owners[:, 1:]
            edges[:-1, :] |= owners[:-1, :] != owners[1:, :]
            edges &= owners != 0
            grid[edges] = 8
            self._frames.append({"mcs": int(self._mcs), "grid": _encode(grid)})

    def step(self, n=1, *a, **k):
        if self._mcs >= self._next:
            self._snap()
            object.__setattr__(self, "_next", self._mcs + self._stride)
        r = object.__getattribute__(self, "_w").step(n, *a, **k)
        object.__setattr__(self, "_mcs", self._mcs + int(n))
        return r

    def _finish(self):
        """Record the final post-run state (after the last step's transitions),
        de-duplicated against the last snapshot's MCS."""
        if not self._frames or self._frames[-1]["mcs"] != self._mcs:
            self._snap()


_CREATED: list[_RecWorld] = []
_ORIG_WFS = build.world_from_spec
_STRIDE = 1


def _patched_world_from_spec(*a, **k):
    w = _ORIG_WFS(*a, **k)
    rec = _RecWorld(w, _STRIDE)
    _CREATED.append(rec)
    return rec


def _run_capture(driver, total_mcs, /, **kwargs):
    """Run a driver under the recorder; return the PRIMARY world's compressed
    frames + dims. The primary world is the one stepped the most (drivers may
    build auxiliary close-contact/control worlds; the main scene has the most
    frames). ``total_mcs`` sets the snapshot stride for ~``N_FRAMES`` frames."""
    global _STRIDE
    _CREATED.clear()
    _STRIDE = max(1, int(round(total_mcs / (N_FRAMES - 1))))
    build.world_from_spec = _patched_world_from_spec
    try:
        driver(**kwargs)
    finally:
        build.world_from_spec = _ORIG_WFS
    primary = max(_CREATED, key=lambda r: len(r._frames))
    primary._finish()
    nx, ny = primary._dims()
    return list(primary._frames), nx, ny, primary._n_epi, primary._n_imm


def _subsample(frames: list, n: int) -> list:
    """~n evenly-spaced frames, always including the first and last."""
    total = len(frames)
    if total <= n:
        return list(frames)
    idxs = sorted({round(i * (total - 1) / (n - 1)) for i in range(n)})
    return [frames[i] for i in idxs]


# ── per-study capture specs ──────────────────────────────────────────────────
def capture_all() -> dict:
    global _STRIDE
    spatial: dict = {}

    # ---- epithelial-sheet-baseline: substrate geometry relaxation (owner) ----
    # No transition driver — build the confluent sheet and relax it under the
    # Potts temperature, recording the full-resolution tile mosaic (owner grid).
    from viva_cpm_studies.influenza import sheet as _sheet
    spec = _sheet.build_sheet_spec(0.3, seed=17)
    relax_per, n_relax = 25, 12
    _STRIDE = max(1, int(round(relax_per * n_relax / (N_FRAMES - 1))))
    world = _RecWorld(_ORIG_WFS(spec, finalize=True), _STRIDE, owner=True)
    for _ in range(n_relax):
        world.step(relax_per)
    world._finish()
    of = _subsample(list(world._frames), N_FRAMES)
    gnx, gny = world._dims()
    spatial["epithelial-sheet-baseline"] = {
        "kind": "owner", "enc": "zlib+b64", "dtype": "uint16",
        "nx": gnx, "ny": gny, "seed": 17,
        "n_epi": int(world._n_epi), "n_imm": 0, "frames": of,
    }
    print(f"sheet: {len(of)} frames {gnx}x{gny} (owner mosaic)")

    # ---- the transition/immune studies: per-site cell-TYPE mosaic ----
    # (slug, driver, kwargs, total_mcs) — total_mcs sets the snapshot stride.
    STUDIES = [
        ("virus-field-infection", run.run_virus_infection,
         dict(patch_mm=0.3, steps=60, seed=17), 600),
        ("ifn-resistance", run.run_virus_infection_with_ifn,
         dict(patch_mm=0.3, steps=60, seed=17), 600),
        ("epithelial-fate", run.run_epithelial_fate,
         dict(patch_mm=0.3, steps=200, seed=17), 2000),
        # Consolidated immune-response scene at TISSUE scale (~900 epithelial
        # cells): macrophage/NK/CD8 chemotaxis toward an infected patch + the
        # chemokine signalling field + contact-killing, replacing the three tiny
        # (~10-cell) single-mechanism demos (Incr 5/6/7).
        ("immune-response", run.run_cytotoxic_response,
         dict(epithelial_cells_per_side=30, n_infected=12, n_macrophages=20,
              n_nk=20, n_cd8=20, margin_sites=6, separation_sites=2,
              steps=60, seed=17), 600),
        ("global-coupling", run.run_global_coupling,
         dict(side=30, steps=20, seed=3), 20),
        # CAPSTONE money-shot: full paper scale (35×35 = 1225-cell epithelial
        # patch), full duration (720 records × 7 MCS = 5040 MCS ≈ 3.5 days), on
        # THIS branch (Increment-10 ROS death active) so the sheet visibly dies.
        # Shared by all three repro-fig* cards (same full_model scene).
        ("repro-full-model", run.run_full_model,
         dict(cells_per_side=35, steps=720, seed=0, init_infection_frac=0.05,
              mcs_per_step=7), 5040),
    ]
    for slug, driver, kwargs, total_mcs in STUDIES:
        frames, nx, ny, n_epi, n_imm = _run_capture(driver, total_mcs, **kwargs)
        frames = _subsample(frames, N_FRAMES)
        # Cell boundaries are baked into the type grid as the EDGE sentinel (8),
        # so every study is outlined with no extra per-frame layer.
        spatial[slug] = {
            "kind": "type", "enc": "zlib+b64", "dtype": "uint8",
            "nx": nx, "ny": ny, "seed": int(kwargs.get("seed", 17)),
            "n_epi": int(n_epi), "n_imm": int(n_imm), "frames": frames,
        }
        span = f'MCS {frames[0]["mcs"]}->{frames[-1]["mcs"]}'
        print(f"{slug}: {len(frames)} frames {nx}x{ny}  ({span}) "
              f"[{n_epi} epi + {n_imm} immune cells]")
    return spatial


def main() -> int:
    spatial = capture_all()
    blob = json.dumps(spatial, separators=(",", ":"))
    out = Path(__file__).resolve().parent / "_influenza_spatial.py"
    header = (
        '"""Baked CPM 2D spatial-state frames for the influenza-sego2022 '
        'studies.\n\n'
        "REAL engine output captured by ``_capture_spatial.py`` (run its module "
        "to\nregenerate). Each study maps to a FULL-RESOLUTION per-site cell-TYPE "
        "grid\n(MEDIUM=0, H=1, I=2, D=3, M=4, K=5, E=6, dormant-reserve=7); the "
        "sheet\nsubstrate stores the owner (cell-label) mosaic instead. Every "
        "frame's grid\nis a zlib-compressed, base64-encoded row-major "
        "``uint8``/``uint16`` lattice\n(decode: "
        "``np.frombuffer(zlib.decompress(base64.b64decode(s)), dtype).reshape"
        "(ny, nx)``).\nRendered by ``influenza_studies.py``'s ``InfluenzaSpatial*"
        "`` cards.\n"
        '"""\n'
    )
    out.write_text(header + "INFLUENZA_SPATIAL = " + repr(blob) + "\n",
                   encoding="utf-8")
    print(f"\nwrote {out} — {len(blob):,} bytes JSON ({len(blob)/1e6:.2f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
