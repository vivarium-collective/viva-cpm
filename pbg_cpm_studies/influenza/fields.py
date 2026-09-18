"""Extracellular virus field (Increment 2, Task 2.1).

Adds a single diffusing/decaying "virus" field to the CPM world and wires
secretion to InfectedReleasing (`types.I`) cells, matching the Sego 2022
source convention documented in `params.yaml`'s `virus:` block.

`world.step(mcs)` advances the fields (secretion + diffusion + decay)
together with the Potts MCS loop -- callers must NOT call advance_fields
separately.
"""
from __future__ import annotations

from .params import load_params
from .sheet import CELL_SIDE_SITES
from . import types

# Expected virus spread length in lattice sites: the source's diffusion
# length is `diffusion_length_cell_diam` (5) cell-diameters, and one cell
# diameter is CELL_SIDE_SITES (5) sites -> 5 * 5 = 25 sites.
_params_at_import = load_params()
VIRUS_DIFFUSION_LENGTH_SITES = int(
    _params_at_import["virus"]["diffusion_length_cell_diam"] * CELL_SIDE_SITES
)

# Per-substep diffusion dt and the number of sub-steps per MCS. Stability for
# explicit forward-Euler diffusion on a 2D lattice requires dt*D*2*ndim < 1
# (ndim=2 here). With D ~ 0.179 lat^2/MCS, dt=1.0 alone would already satisfy
# this (1*0.179*4 ~ 0.72 < 1), but we split each MCS into several smaller
# sub-steps for a more accurate integration, while keeping substeps*dt == 1.0
# so one world.step() MCS still corresponds to exactly one physical unit of
# diffusion time.
_FIELD_DT = 0.2
_FIELD_SUBSTEPS = 5


def add_virus_field(world) -> int:
    """Add the virus field to `world`, wire I-cell secretion, return its index.

    Derivation of the per-pixel secretion rate (see `params.yaml` `virus:`
    block for the full chain back to the CC3D source):
      - `secretion_g_vi` (0.38620664...) is the source's per-cell release
        amount, POST-multiplied by the source's z=2 lattice slab thickness
        (`dim.z`). Our sheet is a single-layer z=1 (pure 2D) field, so we use
        the PRE-dim.z per-cell release: `secretion_g_vi / 2`.
      - That per-cell release is spread (CC3D `secreteInsideCell` semantics)
        evenly over the cell's `cell_sites` (25) pixels, one increment added
        per pixel per MCS: `per_pixel_per_mcs = (secretion_g_vi / 2) / cell_sites`.
    """
    params = load_params()
    virus = params["virus"]
    cell_sites = int(params["cpm"]["cell_sites"])  # 25

    diffusion = float(virus["diffusion_lat2_per_mcs"])
    decay = float(virus["decay_per_mcs"])

    field_idx = world.add_field("virus", diffusion, decay)

    dt, substeps = _FIELD_DT, _FIELD_SUBSTEPS
    assert dt * diffusion * 2 * 2 < 1.0, "virus field diffusion dt is unstable"
    world.set_field_dynamics(field_idx, dt, substeps)

    pre_dimz_g_vi = float(virus["secretion_g_vi"]) / 2.0  # z=1 sheet, see docstring
    per_pixel_per_mcs = pre_dimz_g_vi / cell_sites

    # Engine quirk (crates/cpm-core/src/field.rs `advance_fields`): secretion
    # is added as `rate * dt` exactly ONCE per world.step() MCS (not once per
    # sub-step), using the same `dt` passed to set_field_dynamics above. To
    # land the intended `per_pixel_per_mcs` amount per pixel per MCS
    # regardless of the dt/substeps split chosen for diffusion stability,
    # compensate here: rate = per_pixel_per_mcs / dt, so rate * dt cancels
    # back to per_pixel_per_mcs.
    secretion_rate = per_pixel_per_mcs / dt

    world.set_secretion(field_idx, types.I, secretion_rate)

    return field_idx
