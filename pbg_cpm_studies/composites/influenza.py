"""Modular process-bigraph composites for the influenza-sego2022 model.

The composites are organized as **biological subsystems that compose
bottom-up**, not one composite per figure/study. Each is a reusable,
parameterized ``@composite_generator``; studies reference the subsystem they
need with config rather than getting a bespoke composite:

    epithelium          -- CPM epithelial sheet (the tissue substrate)
      |
    viral_infection     -- + virus field + stochastic H->I infection (on the sheet)
      |
    innate_immunity     -- + macrophages: chemotaxis up virus toward the lesion
      |
    cytotoxic_immunity  -- + NK/CD8+ T: chemotaxis up the macrophage chemokine
      |
    systemic_ode        -- the Price-2015 global (non-spatial) compartment,
                           coupled to the epithelial-infection substrate
      |
    full_model          -- composition of all five subsystems

Each composite returns a process-bigraph *document* wrapping the real Rust CPM
engine (``cpm.processes.cpm_process.CPMProcess``, + ``InfectionProcess`` for the
H->I transition) over a ``load_world`` scene built by the ``influenza.sheet`` /
``influenza.immune`` spec builders. These are dashboard/live-demo scenes: the
quantitative, calibration-pending mechanism (per-cell IL-10-Hill secretion, the
hybrid Price-2015 ODE + dynamic sig_1 feedback, ODE-driven recruitment, and
NK/CD8 contact/nearby killing) is a custom per-MCS loop that is NOT expressible
in the declarative composite spec format -- it runs in the ``influenza.run``
drivers (``run_full_model`` and the per-subsystem ``run_*``), which are what the
studies' measured numbers and the ``influenza.viz`` figures come from. The
composites carry the SCENE (geometry + fields + chemotaxis wiring); the run
drivers carry the mechanism. Not a Fig-3B/5/7 reproduction (that is the capstone
reproduction studies, which drive ``run_full_model``).
"""
from __future__ import annotations

import numpy as np
from process_bigraph.composite_generator import composite_generator

from ..influenza import fields, immune, price_ode, sheet
from ..influenza import types as inf_types
from ..influenza.epithelium_process import _DEFAULT_ENABLE as _EPITHELIUM_DEFAULT_ENABLE
from ..influenza.params import load_params

CPM_ADDR = "local:!cpm.processes.cpm_process.CPMProcess"
INFECTION_ADDR = "local:!pbg_cpm_studies.influenza.infection_process.InfectionProcess"
EPITHELIUM_ADDR = "local:!pbg_cpm_studies.influenza.epithelium_process.EpitheliumProcess"
IMMUNE_ADDR = "local:!pbg_cpm_studies.influenza.immune_process.ImmuneProcess"
ODE_ADDR = "local:!pbg_cpm_studies.influenza.ode_process.SystemicODEProcess"

# Source-faithfulness fix (final-review gap): `run_full_model`'s default
# `enable` (`run.py`'s `_FULL_MODEL_SUBSYSTEMS`) includes "chemokine", which
# gates the uninfected-H IL-10 secretion loop inside
# `EpitheliumProcess._epithelial_fates` (see that module's docstring). The
# composite's epithelium node must enable the SAME epithelial-relevant token
# set or that H-cell IL-10 source is silently off, breaking the IL-10
# negative-feedback loop macrophage secretion depends on
# (`ImmuneProcess.macrophage_secretion_scale`). `EpitheliumProcess` itself
# only inspects a subset of `_FULL_MODEL_SUBSYSTEMS`'s tokens (it has no
# "macrophage"/"nk_cd8"/"recruitment"/"ode" logic -- those are handled by
# `ImmuneProcess`/`SystemicODEProcess` elsewhere in this composite), so the
# full epithelium-relevant set is exactly its own `_DEFAULT_ENABLE` plus the
# one token it additionally understands but excludes by default: "chemokine".
FULL_MODEL_EPITHELIUM_ENABLE = tuple(_EPITHELIUM_DEFAULT_ENABLE) + ("chemokine",)

# modest live-demo aggregate (full-scale 1mm^2 throughput is characterized
# separately by tests/test_influenza_perf.py, not run live from the dashboard)
DEMO_PATCH_MM = 0.1
DEMO_SEED = 17
DEMO_INIT_INFECTED_FRAC = 0.05


def _cpm_store(spec: dict, *, n_fields: int, secretory_types: list[int]) -> dict:
    """The shared CPMProcess store wrapping a ``load_world`` scene ``spec`` --
    identical wiring across every subsystem composite (only the scene, field
    count, and secretory types differ)."""
    return {
        "_type": "process",
        "address": CPM_ADDR,
        "config": {"spec": spec, "mcs_per_update": 10, "n_fields": n_fields,
                   "secretory_types": list(secretory_types)},
        "inputs": {"fates": ["fates"]},
        "outputs": {
            "volumes": ["volumes"],
            "types": ["types"],
            "positions": ["positions"],
            "field_at_cell": ["field_at_cell"],
            "neighbor_secretory": ["neighbor_secretory"],
        },
    }


# ===========================================================================
# 1. epithelium -- the CPM epithelial-sheet substrate every subsystem builds on
# ===========================================================================

def build_spec(patch_mm: float = DEMO_PATCH_MM) -> dict:
    """A ``load_world`` spec for a confluent epithelial sheet (CPM only)."""
    return sheet.build_sheet_spec(patch_mm)


def epithelium_composite_document(patch_mm: float = DEMO_PATCH_MM) -> dict:
    """Composite document: the CPM engine over the confluent epithelial sheet
    (no fields, no immune cells) -- the tissue substrate."""
    return {"cpm": _cpm_store(build_spec(patch_mm), n_fields=0, secretory_types=[]),
            "fates": {}}


@composite_generator(
    name="epithelium", default_n_steps=10,
    description=(
        "Influenza-sego2022 subsystem 1/5 -- the epithelial tissue substrate: a "
        "confluent CPM epithelial sheet (single CPMProcess, no fields, no immune "
        "cells). Every other subsystem composes on top of this. Validates the "
        "substrate geometry, not a biology reproduction."
    ),
    parameters={
        "patch_mm": {"type": "float", "default": DEMO_PATCH_MM,
                     "description": "square patch side length in mm (live demo; small by default)"},
    },
)
def epithelium(core=None, patch_mm: float = DEMO_PATCH_MM) -> dict:
    return epithelium_composite_document(patch_mm)


# ===========================================================================
# 2. viral_infection -- epithelium + virus field + stochastic H->I infection
# ===========================================================================

def build_viral_infection_spec(patch_mm: float = DEMO_PATCH_MM, seed: int = DEMO_SEED,
                               init_infected_frac: float = DEMO_INIT_INFECTED_FRAC) -> dict:
    """Extends the ``epithelium`` sheet scene with the virus field and
    ``init_infected_frac`` of the (seeded-RNG-chosen) H cells pre-set to
    ``types.I`` so the live demo starts from a lesion."""
    spec = sheet.build_sheet_spec(patch_mm, seed=seed)
    spec["fields"] = [fields.virus_field_spec_entry()]
    n_cells = len(spec["cells"])
    n_init = int(round(n_cells * init_infected_frac))
    if n_init > 0:
        rng = np.random.default_rng(seed)
        infected_ids = rng.choice(np.arange(1, n_cells + 1), size=n_init, replace=False)
        for cid in infected_ids:
            spec["cells"][int(cid) - 1]["type"] = inf_types.I
    return spec


def viral_infection_composite_document(patch_mm: float = DEMO_PATCH_MM, seed: int = DEMO_SEED,
                                       init_infected_frac: float = DEMO_INIT_INFECTED_FRAC) -> dict:
    """CPMProcess (epithelium + virus field) + InfectionProcess (the stochastic
    H -> I transition), wired through the shared ``fates`` map (the
    ``cpm.coupling`` down-scale convention). One population-level infection
    process (a single shared RNG for the whole sheet), since the transition is
    a population-wide stochastic draw."""
    spec = build_viral_infection_spec(patch_mm, seed=seed, init_infected_frac=init_infected_frac)
    g_hv = float(load_params()["virus"]["infection_g_hv"])
    return {
        "cpm": _cpm_store(spec, n_fields=1, secretory_types=[inf_types.I]),
        "infection": {
            "_type": "process",
            "address": INFECTION_ADDR,
            "config": {"g_hv": g_hv, "seed": seed},
            "inputs": {"types": ["types"], "field_at_cell": ["field_at_cell"]},
            "outputs": {"fates": ["fates"]},
        },
        "fates": {},
    }


@composite_generator(
    name="viral_infection", default_n_steps=20,
    description=(
        "Influenza-sego2022 subsystem 2/5 -- viral infection on the epithelium: "
        "composes the virus field (secretion/diffusion/decay from "
        "InfectedReleasing cells) + the stochastic H->I infection transition on "
        "top of the epithelial sheet. Live-demo scene; the lesion-spread + "
        "IFN/resistance measurements come from run.run_virus_infection / "
        "run_virus_infection_with_ifn."
    ),
    parameters={
        "patch_mm": {"type": "float", "default": DEMO_PATCH_MM,
                     "description": "square patch side length in mm (live demo; small by default)"},
        "seed": {"type": "int", "default": DEMO_SEED,
                 "description": "RNG seed for both the sheet build and the infection transition"},
        "init_infected_frac": {"type": "float", "default": DEMO_INIT_INFECTED_FRAC,
                               "description": "fraction of H cells seeded as I (InfectedReleasing) at t=0"},
    },
)
def viral_infection(core=None, patch_mm: float = DEMO_PATCH_MM, seed: int = DEMO_SEED,
                    init_infected_frac: float = DEMO_INIT_INFECTED_FRAC) -> dict:
    return viral_infection_composite_document(patch_mm, seed=seed,
                                              init_infected_frac=init_infected_frac)


# ===========================================================================
# 3. innate_immunity -- + macrophages chemotaxing up the virus toward the lesion
# ===========================================================================
# Interior-placed clusters (both the epithelial/infection patch and the
# macrophage cluster >= margin from every noflux wall) -- an unbiased lambda=0
# control geometry (see immune.build_macrophage_scenario_spec). Live-demo scene:
# the measured localization (chemotaxis-on closes the mean distance to the
# infection across 5/5 seeds vs a near-isotropic lambda=0 control) comes from
# run.run_macrophage_response, not this composite.

DEMO_MACROPHAGE_CELLS_PER_SIDE = 4
DEMO_N_INFECTED = 1
DEMO_N_MACROPHAGES = 6
DEMO_MARGIN_SITES = 30
DEMO_SEPARATION_SITES = 25
DEMO_CHEMOTAXIS_V_MACRO = 5000.0  # params.yaml macrophage.chemotaxis_v_macro; pass 0.0 for the control


def build_innate_immunity_spec(*, epithelial_cells_per_side: int = DEMO_MACROPHAGE_CELLS_PER_SIDE,
                               n_infected: int = DEMO_N_INFECTED,
                               n_macrophages: int = DEMO_N_MACROPHAGES,
                               margin_sites: int = DEMO_MARGIN_SITES,
                               separation_sites: int = DEMO_SEPARATION_SITES,
                               seed: int = DEMO_SEED,
                               chemotaxis_v_macro: float = DEMO_CHEMOTAXIS_V_MACRO) -> dict:
    """The macrophage scenario (interior epithelial/infection patch + interior
    macrophage cluster) + virus field with macrophage (``types.M``) chemotaxis
    on the virus field. ``chemotaxis_v_macro=0.0`` reproduces the lambda=0
    control."""
    spec = immune.build_macrophage_scenario_spec(
        epithelial_cells_per_side=epithelial_cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, margin_sites=margin_sites,
        separation_sites=separation_sites, seed=seed)
    field_entry = fields.virus_field_spec_entry()
    field_entry["chemotaxis"] = [{"type": inf_types.M, "lambda": float(chemotaxis_v_macro)}]
    spec["fields"] = [field_entry]
    return spec


def innate_immunity_composite_document(*, epithelial_cells_per_side=DEMO_MACROPHAGE_CELLS_PER_SIDE,
                                       n_infected=DEMO_N_INFECTED, n_macrophages=DEMO_N_MACROPHAGES,
                                       margin_sites=DEMO_MARGIN_SITES,
                                       separation_sites=DEMO_SEPARATION_SITES, seed=DEMO_SEED,
                                       chemotaxis_v_macro=DEMO_CHEMOTAXIS_V_MACRO) -> dict:
    spec = build_innate_immunity_spec(
        epithelial_cells_per_side=epithelial_cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, margin_sites=margin_sites,
        separation_sites=separation_sites, seed=seed,
        chemotaxis_v_macro=chemotaxis_v_macro)
    return {"cpm": _cpm_store(spec, n_fields=1, secretory_types=[inf_types.I]),
            "fates": {}}


@composite_generator(
    name="innate_immunity", default_n_steps=40,
    description=(
        "Influenza-sego2022 subsystem 3/5 -- innate (macrophage) response: adds "
        "macrophages (type M) chemotaxing up the virus field toward the infection "
        "on top of the viral-infection scene (interior-placed clusters, an "
        "unbiased lambda=0-control geometry). Live-demo scene; the measured "
        "localization + chemokine/IL-10 signaling numbers come from "
        "run.run_macrophage_response / run_macrophage_signaling. "
        "chemotaxis_v_macro=0.0 reproduces the lambda=0 control."
    ),
    parameters={
        "chemotaxis_v_macro": {"type": "float", "default": DEMO_CHEMOTAXIS_V_MACRO,
                               "description": "macrophage chemotaxis strength on the virus field (0.0 = lambda=0 control)"},
        "seed": {"type": "int", "default": DEMO_SEED,
                 "description": "RNG seed for the scenario build"},
    },
)
def innate_immunity(core=None, chemotaxis_v_macro: float = DEMO_CHEMOTAXIS_V_MACRO,
                    seed: int = DEMO_SEED) -> dict:
    return innate_immunity_composite_document(seed=seed, chemotaxis_v_macro=chemotaxis_v_macro)


# ===========================================================================
# 4. cytotoxic_immunity -- + NK/CD8+ T chemotaxing up the macrophage chemokine
# ===========================================================================
# Three interior-placed clusters (epithelial/infection patch | macrophage
# cluster | NK+CD8 cluster) + virus field (macrophage chemotaxis) + chemokine
# field (NK/CD8 chemotaxis). Live-demo scene: the contact/nearby-killing
# mechanism (killing.contact_kill_rate/nearby_kill_rate, a custom per-update
# Python calc over world.cell_contact_area_by_type) is NOT in this declarative
# scene -- the localization half is; the killing + measured numbers come from
# run.run_cytotoxic_response.

DEMO_N_NK = 6
DEMO_N_CD8 = 6
DEMO_CYTOTOXIC_MARGIN_SITES = 20
DEMO_CYTOTOXIC_SEPARATION_SITES = 20
DEMO_CHEMOTAXIS_V_NK = 500000.0    # params.yaml nk.chemotaxis_v_nk (5000) x 100 (run.NK_CD8_CHEMOTAXIS_ENGINE_SCALE)
DEMO_CHEMOTAXIS_V_CD8 = 1000000.0  # params.yaml cd8.chemotaxis_v_cd8 (10000) x 100 (same scale, preserves 2:1)


def build_cytotoxic_immunity_spec(*, epithelial_cells_per_side: int = DEMO_MACROPHAGE_CELLS_PER_SIDE,
                                  n_infected: int = DEMO_N_INFECTED,
                                  n_macrophages: int = DEMO_N_MACROPHAGES,
                                  n_nk: int = DEMO_N_NK, n_cd8: int = DEMO_N_CD8,
                                  margin_sites: int = DEMO_CYTOTOXIC_MARGIN_SITES,
                                  separation_sites: int = DEMO_CYTOTOXIC_SEPARATION_SITES,
                                  seed: int = DEMO_SEED,
                                  chemotaxis_v_macro: float = DEMO_CHEMOTAXIS_V_MACRO,
                                  chemotaxis_v_nk: float = DEMO_CHEMOTAXIS_V_NK,
                                  chemotaxis_v_cd8: float = DEMO_CHEMOTAXIS_V_CD8) -> dict:
    """The three-cluster scenario + virus field (macrophage chemotaxis) +
    chemokine field (NK/CD8 chemotaxis), every cluster interior-placed.
    ``chemotaxis_v_nk``/``chemotaxis_v_cd8=0.0`` reproduces the localization
    lambda=0 control; no killing step exists in this declarative scene."""
    spec = immune.build_cytotoxic_scenario_spec(
        epithelial_cells_per_side=epithelial_cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, n_nk=n_nk, n_cd8=n_cd8,
        margin_sites=margin_sites, separation_sites=separation_sites, seed=seed)
    virus_entry = fields.virus_field_spec_entry()
    virus_entry["chemotaxis"] = [{"type": inf_types.M, "lambda": float(chemotaxis_v_macro)}]
    chemo_entry = fields.chemokine_field_spec_entry()
    chemo_entry["chemotaxis"] = [
        {"type": inf_types.K, "lambda": float(chemotaxis_v_nk)},
        {"type": inf_types.E, "lambda": float(chemotaxis_v_cd8)},
    ]
    spec["fields"] = [virus_entry, chemo_entry]
    return spec


def cytotoxic_immunity_composite_document(*, epithelial_cells_per_side=DEMO_MACROPHAGE_CELLS_PER_SIDE,
                                          n_infected=DEMO_N_INFECTED, n_macrophages=DEMO_N_MACROPHAGES,
                                          n_nk=DEMO_N_NK, n_cd8=DEMO_N_CD8,
                                          margin_sites=DEMO_CYTOTOXIC_MARGIN_SITES,
                                          separation_sites=DEMO_CYTOTOXIC_SEPARATION_SITES, seed=DEMO_SEED,
                                          chemotaxis_v_macro=DEMO_CHEMOTAXIS_V_MACRO,
                                          chemotaxis_v_nk=DEMO_CHEMOTAXIS_V_NK,
                                          chemotaxis_v_cd8=DEMO_CHEMOTAXIS_V_CD8) -> dict:
    spec = build_cytotoxic_immunity_spec(
        epithelial_cells_per_side=epithelial_cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, n_nk=n_nk, n_cd8=n_cd8,
        margin_sites=margin_sites, separation_sites=separation_sites, seed=seed,
        chemotaxis_v_macro=chemotaxis_v_macro, chemotaxis_v_nk=chemotaxis_v_nk,
        chemotaxis_v_cd8=chemotaxis_v_cd8)
    return {"cpm": _cpm_store(spec, n_fields=2, secretory_types=[inf_types.I, inf_types.M]),
            "fates": {}}


@composite_generator(
    name="cytotoxic_immunity", default_n_steps=60,
    description=(
        "Influenza-sego2022 subsystem 4/5 -- cytotoxic (NK + CD8+ T) response: "
        "adds NK (type K) + CD8+ (type E) cells chemotaxing up the "
        "macrophage-released chemokine field on top of the innate-immunity "
        "scene (three interior-placed clusters). Live-demo scene -- LOCALIZATION "
        "geometry/chemotaxis only, NO contact-killing step (not expressible in "
        "the declarative scene; see killing.contact_kill_rate / "
        "run.run_cytotoxic_response for the killing mechanism + measured numbers). "
        "chemotaxis_v_nk=chemotaxis_v_cd8=0.0 reproduces the localization control."
    ),
    parameters={
        "chemotaxis_v_macro": {"type": "float", "default": DEMO_CHEMOTAXIS_V_MACRO,
                               "description": "macrophage chemotaxis strength on the virus field"},
        "chemotaxis_v_nk": {"type": "float", "default": DEMO_CHEMOTAXIS_V_NK,
                            "description": "NK chemotaxis strength on the chemokine field (0.0 = lambda=0 control)"},
        "chemotaxis_v_cd8": {"type": "float", "default": DEMO_CHEMOTAXIS_V_CD8,
                             "description": "CD8+ chemotaxis strength on the chemokine field (0.0 = lambda=0 control)"},
        "seed": {"type": "int", "default": DEMO_SEED,
                 "description": "RNG seed for the scenario build"},
    },
)
def cytotoxic_immunity(core=None, chemotaxis_v_macro: float = DEMO_CHEMOTAXIS_V_MACRO,
                       chemotaxis_v_nk: float = DEMO_CHEMOTAXIS_V_NK,
                       chemotaxis_v_cd8: float = DEMO_CHEMOTAXIS_V_CD8,
                       seed: int = DEMO_SEED) -> dict:
    return cytotoxic_immunity_composite_document(
        seed=seed, chemotaxis_v_macro=chemotaxis_v_macro,
        chemotaxis_v_nk=chemotaxis_v_nk, chemotaxis_v_cd8=chemotaxis_v_cd8)


# ===========================================================================
# 5. systemic_ode -- the Price-2015 global (non-spatial) compartment
# ===========================================================================
# The systemic compartment (antibodies, TNF, IL-12, type-II IFN, CD4+, B,
# neutrophils, ROS, APCs) coupled to the epithelial-infection substrate. The
# hybrid Price-2015 ODE integration + spatial<->ODE bidirectional coupling +
# dynamic sig_1 feedback are a custom per-MCS loop (price_ode.GlobalODE), NOT
# expressible in the declarative scene -- this composite carries the coupled
# spatial substrate (the epithelial-infection scene) as the systemic
# compartment's spatial anchor; the ODE mechanism + measured (calibration-
# pending) numbers come from run.run_global_coupling / run_full_model.


def systemic_ode_composite_document(patch_mm: float = DEMO_PATCH_MM, seed: int = DEMO_SEED,
                                    init_infected_frac: float = DEMO_INIT_INFECTED_FRAC) -> dict:
    """The epithelial-infection substrate the systemic Price-2015 ODE couples
    to (reuses the ``viral_infection`` scene). The ODE compartment + coupling
    run in run.run_global_coupling / run_full_model (see module note)."""
    return viral_infection_composite_document(patch_mm, seed=seed,
                                              init_infected_frac=init_infected_frac)


@composite_generator(
    name="systemic_ode", default_n_steps=40,
    description=(
        "Influenza-sego2022 subsystem 5/5 -- the systemic Price-2015 global "
        "(non-spatial) ODE compartment (TNF, ROS, antibody, APC, IL-12, IFN-gamma, "
        "CD4+, B, neutrophils) coupled to the epithelial-infection substrate. The "
        "hybrid ODE integration + spatial<->ODE coupling + dynamic sig_1 feedback "
        "are a custom per-MCS loop (price_ode.GlobalODE), NOT expressible in this "
        "declarative scene -- the composite carries the coupled spatial substrate; "
        "the ODE mechanism + measured (calibration-pending) numbers come from "
        "run.run_global_coupling / run_full_model and viz.global_coupling_figure."
    ),
    parameters={
        "patch_mm": {"type": "float", "default": DEMO_PATCH_MM,
                     "description": "square patch side length in mm (live demo; small by default)"},
        "seed": {"type": "int", "default": DEMO_SEED,
                 "description": "RNG seed for the scenario build"},
        "init_infected_frac": {"type": "float", "default": DEMO_INIT_INFECTED_FRAC,
                               "description": "fraction of H cells seeded as I at t=0"},
    },
)
def systemic_ode(core=None, patch_mm: float = DEMO_PATCH_MM, seed: int = DEMO_SEED,
                 init_infected_frac: float = DEMO_INIT_INFECTED_FRAC) -> dict:
    return systemic_ode_composite_document(patch_mm, seed=seed,
                                           init_infected_frac=init_infected_frac)


# ===========================================================================
# full_model -- composition of all five biological subsystems
# ===========================================================================
# Task 3.2 (influenza-immune-process-composite epic): unlike subsystems 1-5
# above -- which wrap the on-lattice CPMProcess directly and leave the full
# per-MCS mechanism to run.run_full_model's custom Python loop -- full_model
# now wires three REAL process-bigraph Processes (EpitheliumProcess,
# ImmuneProcess, SystemicODEProcess; Phases 1-3 of this epic) into one
# runnable Composite DOCUMENT. This is the "swappable" architecture the epic
# targets: the mechanism itself (fate transitions, off-lattice immune-agent
# chemotaxis/killing/secretion/recruitment, the Price-2015 ODE) now lives
# inside the processes, not a bespoke driver loop -- `run.run_full_model`
# remains the source-faithful reference the capstone reproduction studies
# (repro-fig3b/5/7) use; this composite is the new declarative alternative,
# not (yet) a drop-in replacement for it (see docstring below for what is
# deferred to Task 3.3/3.4).

DEMO_FULL_MODEL_CELLS_PER_SIDE = DEMO_MACROPHAGE_CELLS_PER_SIDE
# margin_sites/separation_sites/mcs_per_step/s_per_mcs match run_full_model's
# own defaults (run.py:1614-1623), not immune.py's build_cytotoxic_scenario_
# spec's larger standalone defaults (30/25) -- this composite is meant to be
# the same scenario run_full_model builds.
DEMO_FULL_MODEL_MARGIN_SITES = 10
DEMO_FULL_MODEL_SEPARATION_SITES = 8
DEMO_FULL_MODEL_MCS_PER_STEP = 7
DEMO_S_PER_MCS = 60.0


def full_model_composite_document(*, cells_per_side: int = DEMO_FULL_MODEL_CELLS_PER_SIDE,
                                  seed: int = DEMO_SEED,
                                  init_infection_frac: float = DEMO_INIT_INFECTED_FRAC,
                                  init_viral_load: float | None = None,
                                  n_macrophages: int = DEMO_N_MACROPHAGES,
                                  n_nk: int = DEMO_N_NK, n_cd8: int = DEMO_N_CD8,
                                  margin_sites: int = DEMO_FULL_MODEL_MARGIN_SITES,
                                  separation_sites: int = DEMO_FULL_MODEL_SEPARATION_SITES,
                                  mcs_per_step: int = DEMO_FULL_MODEL_MCS_PER_STEP,
                                  s_per_mcs: float = DEMO_S_PER_MCS) -> dict:
    """Task 3.2: the process-bigraph Composite document wiring EpitheliumProcess
    + ImmuneProcess + SystemicODEProcess -- see task-3.2-report.md for the full
    store/wiring map. Scene built the SAME way `run.run_full_model` does
    (`immune.build_cytotoxic_scenario_spec` at `epithelial_cells_per_side=
    cells_per_side`, then the `init_infection_frac` central-lesion ->
    sheet-wide-scatter re-seed at `seed+7`, verbatim from run.py:1727-1787).

    Two deliberate departures from `build_cytotoxic_scenario_spec`'s raw
    output, both documented in the report:
      - the epithelium node's CPM `spec` carries ONLY the H/I epithelial
        cells (the scenario's M/K/E CPM-cell blocks are dropped) -- immune
        cells are OFF-LATTICE `ImmuneProcess` agents here, replacing the
        on-lattice M/K/E cells `run_full_model` uses, so keeping both would
        double-count the immune population and leave inert, unchemotaxing
        CPM cells sitting on the lattice (EpitheliumProcess never wires
        `world.set_chemotaxis` for them).
      - those dropped M/K/E cluster cells' seed_block centers instead seed
        the initial `immune_agents` list (one off-lattice agent per CPM-cell
        block the scenario would have placed), so the initial immune
        population starts at the same scenario-designed cluster positions.
    """
    params = load_params()
    tot_cell = int(cells_per_side) * int(cells_per_side)
    # init_viral_load scenario (fig5/fig7): NO pre-infected cells -- infection
    # emerges from a ~uniform virus-field IC that EpitheliumProcess seeds at
    # construction. Mutually exclusive with init_infection_frac (as in
    # run_full_model); when a load is given, force n_infected=0 so the scatter
    # block below is skipped.
    if init_viral_load is not None and float(init_viral_load) > 0.0:
        n_infected = 0
    else:
        n_infected = max(0, min(int(round(float(init_infection_frac) * tot_cell)), tot_cell))

    spec = immune.build_cytotoxic_scenario_spec(
        epithelial_cells_per_side=cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, n_nk=n_nk, n_cd8=n_cd8,
        margin_sites=margin_sites, separation_sites=separation_sites, seed=seed)

    # Sheet-wide infection scatter (run.py:1776-1787, verbatim: `seed+7`,
    # replaces the single central lesion `build_cytotoxic_scenario_spec`
    # places with `n_infected` cells scattered across the whole patch).
    if n_infected > 0:
        epi_idx = [i for i, c in enumerate(spec["cells"])
                   if c["type"] in (inf_types.H, inf_types.I)]
        for i in epi_idx:
            spec["cells"][i]["type"] = inf_types.H
        if epi_idx:
            scatter_rng = np.random.default_rng(seed + 7)
            chosen = scatter_rng.choice(np.array(epi_idx),
                                        size=min(n_infected, len(epi_idx)),
                                        replace=False)
            for i in chosen:
                spec["cells"][int(i)]["type"] = inf_types.I

    # Initial off-lattice immune_agents from the scenario's M/K/E cluster
    # cell blocks (see docstring "departures"), BEFORE those cells are
    # dropped from the epithelium's own spec below.
    immune_agents = []
    next_id = 1
    for c in spec["cells"]:
        if c["type"] in (inf_types.M, inf_types.K, inf_types.E):
            x0, y0, _z0, x1, y1, _z1 = c["seed_block"]
            immune_agents.append({
                "id": next_id, "type": int(c["type"]),
                "x": (x0 + x1) / 2.0, "y": (y0 + y1) / 2.0,
            })
            next_id += 1

    epithelium_spec = dict(spec)
    epithelium_spec["cells"] = [c for c in spec["cells"]
                                if c["type"] in (inf_types.H, inf_types.I)]

    nx, ny, _nz = spec["potts"]["dims"]
    epithelial_side_sites = int(cells_per_side) * sheet.CELL_SIDE_SITES
    # margin_box: the EXTERIOR-to-the-epithelial-patch region -- everything
    # to the right of the patch (where build_cytotoxic_scenario_spec already
    # places the macrophage + NK/CD8 clusters), inside the domain's own
    # margin_sites wall clearance. Does not overlap the epithelial content
    # region (x < margin_sites + epithelial_side_sites).
    margin_box = [margin_sites + epithelial_side_sites, margin_sites,
                 nx - margin_sites, ny - margin_sites]

    consts = price_ode.resolve_constants(params["price_ode"], num_epithelial=tot_cell)

    return {
        # --- pre-initialized top-level stores (existing-composite convention) ---
        "fates": {},
        "field_deposit": [],
        "ode_state": {},
        "immune_kills": [],
        "immune_agents": immune_agents,
        "chemo_field": [],
        "virus_field": [],
        "il10_field": [],
        "dims": [],
        "positions": [],
        "types": [],
        "ode_inputs": {},
        "immune_counts": {},
        # epithelial_positions/infected_ids: the immune process needs the
        # INFECTED cells' positions specifically, which requires filtering
        # epithelium's flat `positions`/`types` list outputs by type -- not
        # expressible as declarative store wiring (no process-bigraph
        # transform/Step for it yet). Seeded empty and left UNWIRED to any
        # process output here; deferred to Task 3.3's driver (see report).
        "epithelial_positions": {},
        "infected_ids": [],
        "recruit_drivers": {},
        "sig_1": 0.0,
        "margin_box": margin_box,
        # Recruitment cadence gating (matching run_full_model's per-MCS
        # recruitment step) is a Task 3.3 driver concern; always-on here so
        # the single-step wiring is exercised.
        "apply_recruitment": True,

        "epithelium": {
            "_type": "process",
            "address": EPITHELIUM_ADDR,
            "config": {
                "spec": epithelium_spec,
                "mcs_per_update": mcs_per_step,
                "mcs_per_step": mcs_per_step,
                "secretory_types": [inf_types.I],
                # Source-faithfulness fix: matches run_full_model's default
                # epithelial subsystem set INCLUDING "chemokine" (see
                # FULL_MODEL_EPITHELIUM_ENABLE docstring above) so the
                # uninfected-H IL-10 gate is active, same as run_full_model.
                "enable": list(FULL_MODEL_EPITHELIUM_ENABLE),
                # fig5/fig7 viral-load IC (0.0 = the init_infection_frac path).
                "init_viral_load": float(init_viral_load or 0.0),
            },
            "inputs": {
                "fates": ["fates"],
                "field_deposit": ["field_deposit"],
                "ode_state": ["ode_state"],
                "immune_kills": ["immune_kills"],
            },
            "outputs": {
                "types": ["types"],
                "positions": ["positions"],
                "field_at_cell_all": ["field_at_cell_all"],
                "chemo_field": ["chemo_field"],
                "virus_field": ["virus_field"],
                "il10_field": ["il10_field"],
                "dims": ["dims"],
                "counts": ["counts"],
                "ode_inputs": ["ode_inputs"],
            },
        },
        "immune": {
            "_type": "process",
            "address": IMMUNE_ADDR,
            "config": {"seed": seed},
            "inputs": {
                "chemo_field": ["chemo_field"],
                "virus_field": ["virus_field"],
                "il10_field": ["il10_field"],
                "dims": ["dims"],
                "immune_agents": ["immune_agents"],
                "epithelial_positions": ["epithelial_positions"],
                "infected_ids": ["infected_ids"],
                "sig_1": ["sig_1"],
                "recruit_drivers": ["recruit_drivers"],
                "margin_box": ["margin_box"],
                "apply_recruitment": ["apply_recruitment"],
            },
            "outputs": {
                "immune_agents": ["immune_agents"],
                "immune_kills": ["immune_kills"],
                "field_deposit": ["field_deposit"],
                "immune_counts": ["immune_counts"],
            },
        },
        "ode": {
            "_type": "process",
            "address": ODE_ADDR,
            "config": {
                "consts": consts,
                "num_epithelial": tot_cell,
                # one ODE step per composite update = mcs_per_step raw MCS
                # of epithelium time (matches EpitheliumProcess's own
                # mcs_per_step granularity; run_full_model instead steps the
                # ODE every single MCS -- a coarser cadence here, deferred/
                # noted for Task 3.3, see report).
                "dt_seconds": mcs_per_step * s_per_mcs,
            },
            "inputs": {
                # H/I/DH/V/F/C/L/B_ei/G_ki: sub-path wiring into
                # epithelium's `ode_inputs` map output (Task 3.2's key
                # wiring decision -- verified to work, see report).
                "H": ["ode_inputs", "H"],
                "I": ["ode_inputs", "I"],
                "DH": ["ode_inputs", "DH"],
                "V": ["ode_inputs", "V"],
                "F": ["ode_inputs", "F"],
                "C": ["ode_inputs", "C"],
                "L": ["ode_inputs", "L"],
                "B_ei": ["ode_inputs", "B_ei"],
                "G_ki": ["ode_inputs", "G_ki"],
                # M/K/E: sub-path wiring into immune's `immune_counts` map
                # output (Task 3.2 controller ruling R3 addendum).
                "M": ["immune_counts", "M"],
                "K": ["immune_counts", "K"],
                "E": ["immune_counts", "E"],
            },
            "outputs": {
                "ode_state": ["ode_state"],
                "recruit_drivers": ["recruit_drivers"],
                "sig_1": ["sig_1"],
            },
        },
    }


@composite_generator(
    name="full_model", default_n_steps=60,
    description=(
        "Influenza-sego2022 -- Task 3.2 (influenza-immune-process-composite "
        "epic): the process-bigraph Composite wiring EpitheliumProcess + "
        "ImmuneProcess (off-lattice M/K/E agents) + SystemicODEProcess (the "
        "Price-2015 global ODE) into one runnable document -- the declarative "
        "'swappable' alternative to run.run_full_model's custom per-MCS driver "
        "loop, which remains the source-faithful reference the capstone "
        "reproduction studies (repro-fig3b/5/7) use."
    ),
    parameters={
        "cells_per_side": {"type": "int", "default": DEMO_FULL_MODEL_CELLS_PER_SIDE,
                           "description": "epithelial patch side length in cells"},
        "seed": {"type": "int", "default": DEMO_SEED,
                 "description": "RNG seed for the scenario build + every process RNG stream"},
        "init_infection_frac": {"type": "float", "default": DEMO_INIT_INFECTED_FRAC,
                               "description": "fraction of epithelial cells seeded as I, scattered sheet-wide"},
    },
)
def full_model(core=None, cells_per_side: int = DEMO_FULL_MODEL_CELLS_PER_SIDE,
               seed: int = DEMO_SEED,
               init_infection_frac: float = DEMO_INIT_INFECTED_FRAC) -> dict:
    return full_model_composite_document(
        cells_per_side=cells_per_side, seed=seed, init_infection_frac=init_infection_frac)
