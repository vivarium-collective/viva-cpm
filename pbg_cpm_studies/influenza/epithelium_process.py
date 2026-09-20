"""EpitheliumProcess (Task 1.2): `CPMProcess` extended with an all-field
per-cell readout, raw-field outputs, and a point-deposit input -- the
spatial-substrate half of the influenza immune-process composite. A later
agent-based immune process (Task 1.3+) samples/deposits these fields through
`field_deposit` (input) and `field_at_cell_all`/`virus_field`/`chemo_field`/
`dims` (outputs).

Field construction (controller ruling R2, task-1.2-brief.md): `load_world`
(`cpm/schema.py`, used by `CPMProcess.initialize`) builds `fields:` entries
declaratively from the spec, but only 2 of the 4 fields this process needs
have a declarative `*_field_spec_entry` helper (`fields.virus_field_spec_entry`,
`fields.chemokine_field_spec_entry` -- no ifn/il10 entries exist). The
reliable, already-tested path is IMPERATIVE construction, the same one every
`run.py` driver uses for this exact field set (see `run_full_model` etc.:
`build.world_from_spec(spec, finalize=False)` -> `fields.add_virus_field` ->
`add_ifn_field` -> `add_chemokine_field` -> `add_il10_field` ->
`world.finalize(seed)`).

Deviation from the brief's literal `initialize` (documented in
task-1.2-report.md): the brief says to call `super().initialize(config)`
first and add fields after. That is not possible -- `CPMProcess.initialize`
calls `cpm.schema.load_world`, which unconditionally calls `world.finalize()`
at the end, and `World.add_field` raises `PyRuntimeError("add_field after
finalize")` once finalized (`crates/cpm-py/src/lib.rs`; confirmed by direct
repro and by `build.world_from_spec`'s own `finalize=False` docstring: "pass
finalize=False when a caller needs to add fields first"). So
`EpitheliumProcess.initialize` does not call `super().initialize(config)`; it
builds the world itself via the same pre-finalize pattern `run.py` uses, then
sets the same `mcs`/`n_fields`/`secretory`/`dims` attributes
`CPMProcess.initialize` would have set. The composite's spec must NOT carry a
`spec["fields"]` key (this process adds fields itself, imperatively).
"""
from cpm.processes.cpm_process import CPMProcess

from . import build, fields


class EpitheliumProcess(CPMProcess):
    """CPMProcess + all-field per-cell readout, raw fields, point deposit.

    Does NOT move the epithelial cell-fate transitions (Task 1.3) -- this
    task only extends the ports/outputs around the existing `fates` ->
    `world.step` -> `volumes/types/positions/field_at_cell/neighbor_secretory`
    behavior `CPMProcess.update` already provides.
    """

    def initialize(self, config):
        spec = self.config["spec"]
        # Imperative field construction, pre-finalize (see module docstring):
        # mirrors run.py's tested `build.world_from_spec(spec, finalize=False)`
        # -> add_*_field -> `world.finalize(seed)` pattern exactly, giving
        # indices 0=virus, 1=ifn, 2=chemokine, 3=il10.
        self.world = build.world_from_spec(spec, finalize=False)
        fields.add_virus_field(self.world)
        fields.add_ifn_field(self.world)
        fields.add_chemokine_field(self.world)
        fields.add_il10_field(self.world)
        self.world.finalize(int(spec["potts"]["seed"]))

        self.mcs = int(self.config["mcs_per_update"])
        # Always exactly the 4 fields added above, regardless of the
        # `n_fields` config value (this process owns field construction).
        self.n_fields = 4
        self.secretory = set(int(t) for t in self.config.get("secretory_types") or [])
        self.dims = self.world.dims()

    def inputs(self):
        base = dict(super().inputs())
        base["field_deposit"] = "list"  # each item [field_idx, x, y, z, amount]
        return base

    def outputs(self):
        out = dict(super().outputs())
        out["field_at_cell_all"] = "overwrite[map[list]]"
        # raw flat fields + dims so a later ImmuneProcess can sample gradients
        # in numpy (process-bigraph has no within-update request/response, so
        # the whole field crosses the store boundary once per update).
        out["chemo_field"] = "overwrite[list]"
        out["virus_field"] = "overwrite[list]"
        out["dims"] = "overwrite[list]"
        return out

    def update(self, state, interval):
        state = state or {}
        for item in state.get("field_deposit") or []:
            fidx, x, y, z, amt = item
            self.world.field_add_source_at(int(fidx), int(x), int(y), int(z), float(amt))

        base = super().update(state, interval)  # applies fates, steps, base outputs

        types = base["types"]; n = len(types)
        base["field_at_cell_all"] = {
            str(cid): [self.world.field_mean_at_cell(f, cid) for f in range(self.n_fields)]
            for cid in range(1, n)
        }
        # field indices (per fields.py order): 0=virus 1=ifn 2=chemokine 3=il10
        base["chemo_field"] = list(self.world.field_conc(2)) if self.n_fields > 2 else []
        base["virus_field"] = list(self.world.field_conc(0)) if self.n_fields > 0 else []
        base["dims"] = list(self.world.dims())
        return base
