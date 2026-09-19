"""Extracellular virus + type-I IFN fields (Increment 2 Task 2.1, Increment 3
Task 3.2).

Adds diffusing/decaying fields to the CPM world and wires secretion to
InfectedReleasing (`types.I`) cells, matching the Sego 2022 source
convention documented in `params.yaml`'s `virus:` / `ifn:` blocks.

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

# Expected IFN spread length in lattice sites: 2 cell-diameters * 5 sites
# per cell-diameter -> 10 sites (shorter than the virus field's 25, since
# IFN's much larger diffusion constant is offset by its much larger decay).
IFN_DIFFUSION_LENGTH_SITES = int(
    _params_at_import["ifn"]["diffusion_length_cell_diam"] * CELL_SIDE_SITES
)

# Per-substep diffusion dt and the number of sub-steps per MCS, for the virus
# field. Stability for explicit forward-Euler diffusion on a 2D lattice
# requires dt*D*2*ndim < 1 (ndim=2 here). With D ~ 0.179 lat^2/MCS, dt=1.0
# alone would already satisfy this (1*0.179*4 ~ 0.72 < 1), but we split each
# MCS into several smaller sub-steps for a more accurate integration, while
# keeping substeps*dt == 1.0 so one world.step() MCS still corresponds to
# exactly one physical unit of diffusion time.
_FIELD_DT = 0.2
_FIELD_SUBSTEPS = 5

# Same idea for the IFN field, but IFN's diffusion constant (D ~ 7.79
# lat^2/MCS) is ~44x the virus field's, so dt=0.2 would be badly unstable
# (0.2*7.79*4 ~ 6.2 >> 1). Use a much smaller dt / more substeps, still with
# substeps*dt == 1.0 MCS. dt=0.025 -> 0.025*7.79*4 ~ 0.779 < 1 (stable, with
# margin), substeps=40 -> 40*0.025 == 1.0.
_IFN_FIELD_DT = 0.025
_IFN_FIELD_SUBSTEPS = 40


def _per_pixel_secretion_rate(secretion_dimz_baked: float, cell_sites: int, dt: float) -> float:
    """Shared secretion-rate derivation for viva-cpm's single-layer (z=1)
    fields, from a CC3D-source per-cell release value that bakes in the
    source's dim.z=2 lattice slab thickness (see `virus.secretion_g_vi` /
    `ifn.secretion_g_fp` in params.yaml for the full derivation chain back to
    the CC3D source). Used by both the virus field (Increment 2) and the IFN
    field (Increment 3) so neither duplicates this convention as a separate
    magic-number computation.

    Two steps:
      1. Halve the z=2-baked per-cell release to get the PRE-dim.z per-cell
         release appropriate for a z=1 (pure 2D) field, then spread it evenly
         over the cell's `cell_sites` pixels (CC3D `secreteInsideCell`
         semantics): one increment added per pixel per MCS.
      2. Engine quirk (crates/cpm-core/src/field.rs `advance_fields`):
         secretion is added as `rate * dt` exactly ONCE per world.step() MCS
         (not once per sub-step), using the same `dt` passed to
         `set_field_dynamics`. To land the intended per-pixel-per-MCS amount
         regardless of the dt/substeps split chosen for diffusion stability,
         divide by dt here so `rate * dt` cancels back to the intended
         per-pixel-per-MCS amount.
    """
    per_pixel_per_mcs = (secretion_dimz_baked / 2.0) / cell_sites
    return per_pixel_per_mcs / dt


def _virus_field_params():
    """Compute (diffusion, decay, dt, substeps, secretion_rate) for the virus
    field. Shared by `add_virus_field` (post-hoc `world.add_field`, used by
    callers that need `finalize=False` -> add field -> `finalize`) and
    `virus_field_spec_entry` (declarative `load_world` spec `fields:` entry,
    used by composites built from a plain spec dict) so both routes apply
    the exact same numbers.

    See `_per_pixel_secretion_rate` for the secretion-rate derivation.
    """
    params = load_params()
    virus = params["virus"]
    cell_sites = int(params["cpm"]["cell_sites"])  # 25

    diffusion = float(virus["diffusion_lat2_per_mcs"])
    decay = float(virus["decay_per_mcs"])

    dt, substeps = _FIELD_DT, _FIELD_SUBSTEPS
    assert dt * diffusion * 2 * 2 < 1.0, "virus field diffusion dt is unstable"
    assert substeps * dt == 1.0, "virus field substeps*dt must equal one MCS"

    secretion_rate = _per_pixel_secretion_rate(float(virus["secretion_g_vi"]), cell_sites, dt)

    return diffusion, decay, dt, substeps, secretion_rate


def _ifn_field_params():
    """Compute (diffusion, decay, dt, substeps, secretion_rate) for the
    type-I IFN field. Mirrors `_virus_field_params`, but IFN's much larger
    diffusion constant (D ~ 7.79 vs virus' ~0.179 lat^2/MCS) requires a
    smaller dt / more substeps for stability -- see `_IFN_FIELD_DT` /
    `_IFN_FIELD_SUBSTEPS` above. Secretion here is only the LOCAL
    infected-cell basal release (`ifn.secretion_g_fp`, source symbol
    `b_fi`) -- NOT the global APC-amplification `b_fp*P` term, which is out
    of scope (deferred to Increment 8 per params.yaml's `ifn:` block notes).

    See `_per_pixel_secretion_rate` for the secretion-rate derivation (same
    z=1/per-pixel convention as the virus field).
    """
    params = load_params()
    ifn = params["ifn"]
    cell_sites = int(params["cpm"]["cell_sites"])  # 25

    diffusion = float(ifn["diffusion_lat2_per_mcs"])
    decay = float(ifn["decay_per_mcs"])

    dt, substeps = _IFN_FIELD_DT, _IFN_FIELD_SUBSTEPS
    assert dt * diffusion * 2 * 2 < 1.0, "IFN field diffusion dt is unstable"
    assert substeps * dt == 1.0, "IFN field substeps*dt must equal one MCS"

    secretion_rate = _per_pixel_secretion_rate(float(ifn["secretion_g_fp"]), cell_sites, dt)

    return diffusion, decay, dt, substeps, secretion_rate


def add_virus_field(world) -> int:
    """Add the virus field to `world`, wire I-cell secretion, return its index.

    See `_virus_field_params` for the secretion-rate derivation.
    """
    diffusion, decay, dt, substeps, secretion_rate = _virus_field_params()

    field_idx = world.add_field("virus", diffusion, decay)
    world.set_field_dynamics(field_idx, dt, substeps)
    world.set_secretion(field_idx, types.I, secretion_rate)

    return field_idx


def add_ifn_field(world) -> int:
    """Add the type-I IFN field to `world`, wire I-cell secretion, return its
    index. Mirrors `add_virus_field`; see `_ifn_field_params` for the
    diffusion/decay/dt/substeps/secretion-rate derivation.
    """
    diffusion, decay, dt, substeps, secretion_rate = _ifn_field_params()

    field_idx = world.add_field("ifn", diffusion, decay)
    world.set_field_dynamics(field_idx, dt, substeps)
    world.set_secretion(field_idx, types.I, secretion_rate)

    return field_idx


def virus_field_spec_entry() -> dict:
    """Declarative `load_world` spec `fields:` entry for the virus field --
    equivalent to `add_virus_field`, for composites/drivers that build the
    world from a plain spec dict (`cpm.schema.load_world`) rather than
    driving `cpm_core.World` directly (see `pbg_cpm_studies.composites.
    influenza.virus_infection`)."""
    diffusion, decay, dt, substeps, secretion_rate = _virus_field_params()
    return {
        "name": "virus", "d": diffusion, "decay": decay,
        "secretion": [{"type": types.I, "rate": secretion_rate}],
        "dynamics": {"dt": dt, "substeps": substeps},
    }
