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

from ..influenza import fields, immune, sheet
from ..influenza import types as inf_types
from ..influenza.params import load_params

CPM_ADDR = "local:!cpm.processes.cpm_process.CPMProcess"
INFECTION_ADDR = "local:!pbg_cpm_studies.influenza.infection_process.InfectionProcess"

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
# The complete multiscale scene: epithelium + viral infection + innate + cytotoxic
# immune cells + the systemic ODE compartment. The full per-MCS mechanism
# (infection/death/Allee transitions, per-cell IL-10-Hill secretion, the hybrid
# Price-2015 ODE + dynamic sig_1, ODE-driven recruitment, NK/CD8 contact+nearby
# killing) runs in run.run_full_model -- which the capstone reproduction studies
# (repro-fig3b / repro-fig5 / repro-fig7) drive and measure against the digitized
# Fig 3B/5/7 acceptance bands. This composite carries the full spatial scene.


def full_model_composite_document(*, epithelial_cells_per_side=DEMO_MACROPHAGE_CELLS_PER_SIDE,
                                  n_infected=DEMO_N_INFECTED, n_macrophages=DEMO_N_MACROPHAGES,
                                  n_nk=DEMO_N_NK, n_cd8=DEMO_N_CD8,
                                  margin_sites=DEMO_CYTOTOXIC_MARGIN_SITES,
                                  separation_sites=DEMO_CYTOTOXIC_SEPARATION_SITES, seed=DEMO_SEED,
                                  chemotaxis_v_macro=DEMO_CHEMOTAXIS_V_MACRO,
                                  chemotaxis_v_nk=DEMO_CHEMOTAXIS_V_NK,
                                  chemotaxis_v_cd8=DEMO_CHEMOTAXIS_V_CD8) -> dict:
    """The complete spatial scene (all immune subsystems); the full mechanism +
    the systemic ODE run in run.run_full_model (see module note). Reuses the
    ``cytotoxic_immunity`` scene as the spatial substrate for all subsystems."""
    return cytotoxic_immunity_composite_document(
        epithelial_cells_per_side=epithelial_cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, n_nk=n_nk, n_cd8=n_cd8,
        margin_sites=margin_sites, separation_sites=separation_sites, seed=seed,
        chemotaxis_v_macro=chemotaxis_v_macro, chemotaxis_v_nk=chemotaxis_v_nk,
        chemotaxis_v_cd8=chemotaxis_v_cd8)


@composite_generator(
    name="full_model", default_n_steps=60,
    description=(
        "Influenza-sego2022 -- the full multiscale model: the composition of all "
        "five biological subsystems (epithelium + viral_infection + innate_immunity "
        "+ cytotoxic_immunity + systemic_ode). This composite carries the complete "
        "spatial scene; the full per-MCS mechanism (infection/death/Allee, per-cell "
        "IL-10-Hill secretion, the hybrid Price-2015 ODE + dynamic sig_1, ODE-driven "
        "recruitment, NK/CD8 contact+nearby killing) runs in run.run_full_model, "
        "which the capstone reproduction studies (repro-fig3b/5/7) drive and measure "
        "against the digitized Fig 3B/5/7 acceptance bands."
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
def full_model(core=None, chemotaxis_v_macro: float = DEMO_CHEMOTAXIS_V_MACRO,
               chemotaxis_v_nk: float = DEMO_CHEMOTAXIS_V_NK,
               chemotaxis_v_cd8: float = DEMO_CHEMOTAXIS_V_CD8,
               seed: int = DEMO_SEED) -> dict:
    return full_model_composite_document(
        seed=seed, chemotaxis_v_macro=chemotaxis_v_macro,
        chemotaxis_v_nk=chemotaxis_v_nk, chemotaxis_v_cd8=chemotaxis_v_cd8)
