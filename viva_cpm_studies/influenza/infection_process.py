"""``InfectionProcess`` -- a population-level process-bigraph `Process` that
wraps the pure `transitions.infection_step` function so the stochastic H -> I
infection transition can be wired into a composite alongside `CPMProcess`.

Unlike the per-cell subcell pattern in `cpm.subcellular.receptor`
(one `Process` per cell, deterministic), infection is a single stochastic
draw per H cell that needs one shared `numpy.random.Generator` -- so this is
ONE process for the whole population, reading the CPM's `types` and
`field_at_cell` outputs and writing the *entire* `fates` map back each
update (see `cpm.processes.cpm_process.CPMProcess.update`: a fate of 0 --
including a cell id simply absent from the dict -- means "no override").
"""
from __future__ import annotations

import numpy as np
from process_bigraph import Process

from bigraph_schema.contract import ProcessContract

from . import transitions, types
from .params import load_params

# Resolved from params.yaml at import time (single source of truth for the
# rate coefficient) rather than duplicated as a hardcoded literal here; the
# `virus_infection` composite always passes `g_hv` explicitly, so this only
# matters as the schema's documented default.
_DEFAULT_G_HV = float(load_params()["virus"]["infection_g_hv"])


class InfectionProcess(Process):
    """Stochastic H -> I viral-infection transition for the influenza-Sego2022 CPM.

    One population-level `Process` that, each MCS, reads the CPM's per-cell
    type codes and the per-cell mean local virus-field concentration, draws an
    infection event for every healthy (H) epithelial cell, and writes back the
    set of cells that just became infected (I) as a `fates` override map. A
    single shared RNG is used for the whole population so the draws are
    reproducible from `seed`.
    """

    description = (
        "Influenza H -> I infection transition (Sego 2022 cellularization).\n\n"
        "For each healthy epithelial cell, infection in one MCS is a Poisson "
        "event whose rate is proportional to the local virus concentration:\n"
        "    rate = g_hv * v_bar\n"
        "    Pr(infect) = 1 - exp(-rate)\n"
        "The cell flips to I iff rng.random() < Pr(infect). Cells that are "
        "already I/D, immune, or the index-0 medium sentinel pass through "
        "unchanged, and only newly-infected cells appear in the fates map "
        "(an absent cell id means 'no override')."
    )

    contract = ProcessContract(
        summary=(
            "Stochastic healthy->infected (H->I) epithelial transition for the "
            "influenza-Sego2022 CPM. Once per MCS, each H cell is infected with "
            "probability 1 - exp(-g_hv * v_bar), where v_bar is the mean virus-field "
            "concentration over that cell and g_hv is the per-MCS infection rate "
            "coefficient. Draws use one shared seeded RNG for the whole population, "
            "and the process emits only the cell ids that just became infected as a "
            "fates override map consumed by CPMProcess."
        ),
        description=(
            "Wraps the pure `transitions.infection_step` function so the stochastic "
            "H->I transition can be wired into a composite alongside CPMProcess. Infection "
            "is modeled as a single process for the whole epithelial population (rather than "
            "one Process per cell) because each cell's draw needs the same shared RNG and the "
            "entire fates map is rewritten each update. The virus field itself (diffusion, "
            "secretion by infected cells, decay) is owned by other components; this process "
            "only reads the already-sampled per-cell virus concentration and converts it into "
            "infection events."
        ),
        inputs={
            "types": (
                "list[integer] of per-cell type codes indexed by cell id "
                "(0=medium sentinel, 1=H healthy, 2=I infected, 3=D dead, 4=M, 5=K, 6=E). "
                "Index 0 is always left unchanged; only cells currently equal to H (1) are "
                "candidates for infection this step."
            ),
            "field_at_cell": (
                "map[float] from cell-id string to the mean local virus-field concentration "
                "seen by that cell (field amount / cell volume, lattice units). Entries are "
                "matched to `types` by integer cell id; cells with no entry default to 0 virus "
                "(no infection pressure)."
            ),
        },
        outputs={
            "fates": (
                "overwrite[map[integer]] from cell-id string to the new type code (always "
                "types.I = 2) for the cells that transitioned H->I this step. The map is an "
                "OVERWRITE replacing the prior fates each update; a cell id absent from the map "
                "means 'no fate override' (fate 0), so already-infected and non-H cells are "
                "never listed."
            ),
        },
        config={
            "g_hv": (
                "Per-MCS infection rate coefficient g_hv (default from params.yaml "
                "virus.infection_g_hv). Multiplies the per-cell mean virus concentration to "
                "give the infection hazard rate for one MCS."
            ),
            "seed": (
                "Integer seed for the single numpy Generator shared across all per-cell "
                "infection draws, making the stochastic transition reproducible."
            ),
        },
        math=[
            "rate_i = g_hv * v_bar_i",
            "Pr(infect_i) = 1 - exp(-rate_i)",
            "cell i (type H) -> I  iff  U_i < Pr(infect_i),  U_i ~ Uniform(0, 1)",
        ],
        symbols={
            "g_hv": "per-MCS infection rate coefficient (1/MCS per unit virus concentration)",
            "v_bar_i": "mean local virus-field concentration over cell i (virus amount / cell volume, lattice units)",
            "rate_i": "per-MCS infection hazard for cell i, g_hv * v_bar_i (1/MCS)",
            "Pr(infect_i)": "probability cell i is infected during one MCS (dimensionless, 0..1)",
            "U_i": "uniform random draw for cell i on [0, 1) (dimensionless)",
            "interval": "update window passed by the scheduler; one call advances one MCS (60 s/MCS in Sego 2022)",
        },
        assumptions=[
            "Infection is a single stochastic draw per healthy (H) cell per MCS; one shared RNG serves the whole population, so this is one Process for all cells rather than one per cell.",
            "Only cells currently of type H (1) can be infected; already-infected (I), dead (D), immune (M/K/E), and the index-0 medium sentinel are copied through unchanged.",
            "The per-MCS infection probability is the exponential survival form Pr = 1 - exp(-g_hv * v_bar), equivalent to a constant-hazard Poisson event over the step.",
            "v_bar is supplied read-only via field_at_cell; this process does not modify the virus field (no deposition or decay here). Cells missing from field_at_cell are treated as seeing zero virus.",
            "fates is an overwrite map rewritten in full each update; a fate of 0 (or an absent cell id) means CPMProcess applies no type override to that cell.",
        ],
        references=[
            "Sego, T.J. et al. (2022) A modular framework for multiscale, multicellular, "
            "spatiotemporal modeling of acute primary viral infection and immune response in "
            "epithelial tissues. J. Theor. Biol. 532:110918 (CC3D ViralInfectionVTM / ImmuneModel).",
        ],
    )

    config_schema = {
        "g_hv": {"_type": "float", "_default": _DEFAULT_G_HV},
        "seed": {"_type": "integer", "_default": 17},
    }

    def initialize(self, config):
        self.g_hv = float(config["g_hv"])
        self.rng = np.random.default_rng(int(config["seed"]))

    def inputs(self):
        return {"types": "list[integer]", "field_at_cell": "map[float]"}

    def outputs(self):
        return {"fates": "overwrite[map[integer]]"}

    def update(self, state, interval):
        state = state or {}
        types_list = list(state.get("types") or [])
        field_at_cell = state.get("field_at_cell") or {}
        if not types_list:
            return {"fates": {}}

        virus_at_cell = [0.0] * len(types_list)
        for cid_key, v in field_at_cell.items():
            cid = int(cid_key)
            if 0 <= cid < len(types_list):
                virus_at_cell[cid] = float(v)

        new_types = transitions.infection_step(types_list, virus_at_cell, self.g_hv, self.rng)

        fates = {
            str(cid): types.I
            for cid in range(1, len(types_list))
            if new_types[cid] == types.I and types_list[cid] != types.I
        }
        return {"fates": fates}
