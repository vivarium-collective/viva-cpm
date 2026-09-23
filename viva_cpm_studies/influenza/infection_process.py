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

from . import transitions, types
from .params import load_params

# Resolved from params.yaml at import time (single source of truth for the
# rate coefficient) rather than duplicated as a hardcoded literal here; the
# `virus_infection` composite always passes `g_hv` explicitly, so this only
# matters as the schema's documented default.
_DEFAULT_G_HV = float(load_params()["virus"]["infection_g_hv"])


class InfectionProcess(Process):
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
