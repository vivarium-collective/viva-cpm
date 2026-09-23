"""Extracellular virus + type-I IFN fields (Increment 2 Task 2.1, Increment 3
Task 3.2), plus the macrophage-released chemokine + IL-10 fields (Increment 6
Task 6.1).

Adds diffusing/decaying fields to the CPM world and wires secretion to
InfectedReleasing (`types.I`) cells (virus/IFN) or macrophage (`types.M`) /
uninfected (`types.H`) cells (chemokine/IL-10), matching the Sego 2022
source convention documented in `params.yaml`'s `virus:` / `ifn:` /
`chemokine:` / `il10:` blocks.

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

# Chemokine field (Increment 6 Task 6.1): D ~ 15.608 lat^2/MCS is the LARGEST
# of any field so far (~2x IFN's ~7.79) -- dt*D*2*ndim<1 needs dt < 1/(15.608*4)
# ~ 0.016. dt=0.0125 -> 0.0125*15.608*4 ~ 0.780 < 1 (same stability margin as
# the IFN field's dt=0.025 choice above), substeps=80 -> 80*0.0125 == 1.0 MCS.
_CHEMOKINE_FIELD_DT = 0.0125
_CHEMOKINE_FIELD_SUBSTEPS = 80

# IL-10 field (Increment 6 Task 6.1): D ~ 4.901 lat^2/MCS, similar order to
# the IFN field's ~7.79 (brief: "similar to IFN"). dt=0.04 -> 0.04*4.901*4 ~
# 0.784 < 1 (same ~0.78 margin convention as the IFN/chemokine fields above),
# substeps=25 -> 25*0.04 == 1.0 MCS.
_IL10_FIELD_DT = 0.04
_IL10_FIELD_SUBSTEPS = 25


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


def _chemokine_field_params():
    """(diffusion, decay, dt, substeps, secretion_rate) for the chemokine
    field. Secretion here is only the macrophage BASE rate (`chemokine.b_c`,
    per-pixel via `_per_pixel_secretion_rate`) -- the IL-10-Hill per-cell
    REGULATION (`signaling.macrophage_secretion_scale`) is applied on top
    via `world.set_cell_secretion_scale` by the run driver each update, not
    baked into this base rate. No global boundary term (deferred to
    Increment 8, see `params.yaml`'s `il10.deferred.global_boundary_secretion`
    -- the SAME deferred item covers both chemokine and IL-10's global block).
    """
    params = load_params()
    chemo = params["chemokine"]
    cell_sites = int(params["cpm"]["cell_sites"])

    diffusion = float(chemo["diffusion_lat2_per_mcs"])
    decay = float(chemo["decay_per_mcs"])

    dt, substeps = _CHEMOKINE_FIELD_DT, _CHEMOKINE_FIELD_SUBSTEPS
    assert dt * diffusion * 2 * 2 < 1.0, "chemokine field diffusion dt is unstable"
    assert substeps * dt == 1.0, "chemokine field substeps*dt must equal one MCS"

    secretion_rate = _per_pixel_secretion_rate(float(chemo["b_c"]), cell_sites, dt)

    return diffusion, decay, dt, substeps, secretion_rate


def _il10_field_params():
    """(diffusion, decay, dt, substeps, macrophage_secretion_rate,
    uninfected_secretion_rate) for the IL-10 field. IL-10 has TWO base
    secreting types (source: `ChemokineSecretionSteppable`'s sibling
    `IL10SecretionSteppable`): macrophages (`il10.b_l`, Hill-regulated same
    as chemokine) AND uninfected (H) cells (`il10.mu_l * il10.b_lh`, gated
    by `signaling.uninfected_il10_scale` = `1 - resist`). Both base rates go
    through the same `_per_pixel_secretion_rate` helper -- `mu_l*b_lh`
    carries exactly one dim.z factor (baked into `b_lh` alone, `mu_l` has
    none per params.yaml's comment), same convention as `b_c`/`b_l` alone.
    """
    params = load_params()
    il10 = params["il10"]
    cell_sites = int(params["cpm"]["cell_sites"])

    diffusion = float(il10["diffusion_lat2_per_mcs"])
    decay = float(il10["decay_per_mcs"])

    dt, substeps = _IL10_FIELD_DT, _IL10_FIELD_SUBSTEPS
    assert dt * diffusion * 2 * 2 < 1.0, "IL-10 field diffusion dt is unstable"
    assert substeps * dt == 1.0, "IL-10 field substeps*dt must equal one MCS"

    macrophage_rate = _per_pixel_secretion_rate(float(il10["b_l"]), cell_sites, dt)
    uninfected_rate = _per_pixel_secretion_rate(
        float(il10["mu_l"]) * float(il10["b_lh"]), cell_sites, dt)

    return diffusion, decay, dt, substeps, macrophage_rate, uninfected_rate


def il10_hill_constants() -> tuple[float, float, float, float]:
    """`(sig_1, g_1, g_2, d_2)`, all CELLULARIZED, for `signaling.
    macrophage_secretion_scale` -- the shared IL-10-Hill regulation factor
    gating both macrophage chemokine and macrophage IL-10 secretion.

    - `sig_1` = `il10.sig_1_stub` (see that key's params.yaml comment for the
      Task-6.1 sig_1-STUB ruling: the source's true `sig_1 := a_11*T+a_12*D`
      is a dynamic ODE quantity, deferred to Increment 8; the L-dependent
      Hill self-regulation shape and the resulting chemokine gradient shape
      are robust to this stub's exact value).
    - `g_1`, `g_2` are recorded RAW (pre-cellularization) in params.yaml
      because they carry a scenario-dependent `s_v` (= eta =
      num_epithelial/tot_ec_ODE) factor; `d_2` is recorded RAW too (needs
      `s_l` only, no `s_v`, but multiplying by 1.0 first keeps one formula).
      This driver defaults to the `scaling.eta_0p3mm` (0.0049) scaling
      factor -- the SAME "default to the 0.3mm scenario" convention every
      other driver in this package already uses (e.g. `run.
      run_virus_infection`'s `patch_mm: float = 0.3` default) -- rather than
      computing a scenario-exact eta for the non-confluent macrophage
      scenario (`immune.build_macrophage_scenario_spec`, which has no
      `patch_mm` of its own): a documented approximation, consistent with
      that scenario's own "2D approximation" status (see `immune.py`'s
      module docstring).
    """
    params = load_params()
    il10 = params["il10"]
    scaling = params["scaling"]

    s_l = 1.0 / (float(scaling["ode_epithelial_population"]) * float(params["cpm"]["cell_sites"]))
    s_v = float(scaling["eta_0p3mm"])

    sig_1 = float(il10["sig_1_stub"])
    g_1 = float(il10["g_1"]) * s_v
    g_2 = float(il10["g_2"]) * s_v * s_l
    d_2 = float(il10["d_2"]) * s_l

    return sig_1, g_1, g_2, d_2


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


def add_chemokine_field(world) -> int:
    """Add the chemokine field to `world`, wire the macrophage BASE
    secretion rate, return its index. Mirrors `add_virus_field`/
    `add_ifn_field`; see `_chemokine_field_params` for the diffusion/decay/
    dt/substeps/secretion-rate derivation.

    Callers must additionally apply the per-cell IL-10-Hill regulation each
    update via `world.set_cell_secretion_scale(field_idx, macrophage_cid,
    signaling.macrophage_secretion_scale(...))` (see `run.py`'s signaling
    driver) -- the base rate set here is the UNREGULATED (scale=1.0)
    per-pixel rate; `world`'s default per-cell scale is 1.0 until a caller
    sets it otherwise, so a caller that never calls `set_cell_secretion_scale`
    gets the unregulated (Hill=1) rate, not zero.
    """
    diffusion, decay, dt, substeps, secretion_rate = _chemokine_field_params()

    field_idx = world.add_field("chemokine", diffusion, decay)
    world.set_field_dynamics(field_idx, dt, substeps)
    world.set_secretion(field_idx, types.M, secretion_rate)

    return field_idx


def add_il10_field(world) -> int:
    """Add the IL-10 field to `world`, wire BOTH base secretion sources
    (macrophage `types.M` -- Hill-regulated, same convention as chemokine;
    uninfected `types.H` -- gated by `1 - resist`), return its index. See
    `_il10_field_params` for the diffusion/decay/dt/substeps/secretion-rate
    derivation. As with `add_chemokine_field`, the per-cell regulation scales
    (Hill for M, `1-resist` for H) are applied by the run driver each update
    via `world.set_cell_secretion_scale` -- this only wires the unregulated
    (scale=1.0) base per-pixel rates.
    """
    diffusion, decay, dt, substeps, macrophage_rate, uninfected_rate = _il10_field_params()

    field_idx = world.add_field("il10", diffusion, decay)
    world.set_field_dynamics(field_idx, dt, substeps)
    world.set_secretion(field_idx, types.M, macrophage_rate)
    world.set_secretion(field_idx, types.H, uninfected_rate)

    return field_idx


def virus_field_spec_entry() -> dict:
    """Declarative `load_world` spec `fields:` entry for the virus field --
    equivalent to `add_virus_field`, for composites/drivers that build the
    world from a plain spec dict (`cpm.schema.load_world`) rather than
    driving `cpm_core.World` directly (see `viva_cpm_studies.composites.
    influenza.viral_infection`)."""
    diffusion, decay, dt, substeps, secretion_rate = _virus_field_params()
    return {
        "name": "virus", "d": diffusion, "decay": decay,
        "secretion": [{"type": types.I, "rate": secretion_rate}],
        "dynamics": {"dt": dt, "substeps": substeps},
    }


def chemokine_field_spec_entry() -> dict:
    """Declarative `load_world` spec `fields:` entry for the chemokine field
    -- equivalent to `add_chemokine_field`, for composites/drivers that
    build the world from a plain spec dict (`cpm.schema.load_world`) rather
    than driving `cpm_core.World` directly (see `viva_cpm_studies.composites.
    influenza.cytotoxic_immunity`, Task 7.3's live-demo wrapper). This is
    the UNREGULATED (IL-10-Hill scale=1.0) base macrophage secretion rate,
    same caveat as `add_chemokine_field` -- a caller that wants the
    Increment-6 per-cell regulation must still apply it via
    `set_cell_secretion_scale` after world construction (not available
    declaratively in the `load_world` spec format)."""
    diffusion, decay, dt, substeps, secretion_rate = _chemokine_field_params()
    return {
        "name": "chemokine", "d": diffusion, "decay": decay,
        "secretion": [{"type": types.M, "rate": secretion_rate}],
        "dynamics": {"dt": dt, "substeps": substeps},
    }
