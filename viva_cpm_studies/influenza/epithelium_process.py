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

Task 1.3 (this revision) moves `run.py`'s `run_full_model` per-MCS epithelial
cell-fate pipeline (`_one_mcs` steps 2-8.5, `run.py:1975`-`2097`) INTO
`update()`, so every raw MCS advanced by this process runs infection,
resistance, secretion, Allee and ROS fate transitions in the exact source
order -- see `_epithelial_fates` below and task-1.3-report.md for the full
port/RNG/ode_inputs accounting.
"""
import math

import numpy as np

from bigraph_schema.contract import ProcessContract

from viva_cpm.processes.cpm_process import CPMProcess

from . import allee, build, fields, price_ode, signaling, transitions, types
from .params import load_params
from .resistance import cell_resistance
from .run import _EPITHELIAL_TYPES, _epithelial_contact_totals

# Mirrors run_full_model's `enable` set (task-1.3-brief.md): this task only
# implements the EPITHELIAL-side tokens by default. "macrophage"/"nk_cd8"/
# "recruitment"/"ode" (chemotaxis wiring, ODE stepping) are entirely out of
# scope for Task 1.3 -- deferred to the later ImmuneProcess/ODE-process
# tasks. "chemokine" IS understood (it gates the macrophage-independent
# uninfected-H IL-10 secretion loop, matching run.py:1999-2010's nesting --
# see `_epithelial_fates`), but is NOT in the default set, matching the
# brief's literal 6-token default; pass it explicitly in `enable` to turn
# that IL-10 secretion on (see task-1.3-report.md "Deviations"/fix log).
_DEFAULT_ENABLE = ("infection", "ifn", "death", "ros", "allee", "killing")


class EpitheliumProcess(CPMProcess):
    """CPMProcess + all-field per-cell readout, raw fields, point deposit,
    and the ordered epithelial cell-fate pipeline (Task 1.3).

    Every `update()` call advances `mcs_per_step` raw MCS; each raw MCS runs
    `world.step(1)` followed by the full source-ordered fate pipeline
    (`_epithelial_fates`), matching `run_full_model`'s `_one_mcs` cadence
    exactly (fates fire EVERY MCS, not once per `update()` call). This is why
    `update()` no longer delegates to `CPMProcess.update` (which does one
    big `world.step(mcs_per_update)` with no interleaved fates) -- the base
    per-cell/neighbor outputs are recomputed here instead, after the loop.
    """

    contract = ProcessContract(
        summary=(
            "Epithelial-sheet half of the influenza-Sego2022 immune composite: a Cellular "
            "Potts (Potts/CPM) lattice of epithelial cells carrying a viral-infection state "
            "(healthy H -> infected I -> dead D, with D->H Allee recovery) coupled to four "
            "diffusing fields (virus, IFN, chemokine, IL-10). Each update() advances mcs_per_step "
            "raw Monte-Carlo steps; every MCS runs world.step(1) then a source-ordered cell-fate "
            "pipeline (infection, immune kills, IFN resistance, field secretion, Allee death/recovery, "
            "infected apoptosis, ROS death) so cell state and the virus/signalling fields co-evolve "
            "in lockstep. It also exposes the raw fields and spatial->ODE integrals a downstream "
            "immune/global-ODE process consumes."
        ),
        description=(
            "Reproduces run.py run_full_model's per-MCS epithelial pipeline (_one_mcs steps 2-8.5) "
            "inside a process-bigraph Process, extending CPMProcess with all-field per-cell readouts, "
            "raw field outputs, a point-deposit input, and spatial->ODE coupling outputs. Which fate "
            "tokens fire is gated by the `enable` config set."
        ),
        inputs={
            "fates": (
                "Inherited from CPMProcess. Map keyed by string cell id -> target cell type; an "
                "external (e.g. subcellular) process can override a single cell's fate. Applied once "
                "per update() call before the MCS loop (cid>0 only)."
            ),
            "field_deposit": (
                "List of point deposits, each [field_idx, x, y, z, amount], applied once per update() "
                "via world.field_add_source_at so a downstream agent-based immune process can inject "
                "into any field (e.g. cytokine release at a macrophage site)."
            ),
            "ode_state": (
                "Well-mixed global-ODE state map pushed in each update; only the ROS oxidant scalar "
                "X is consumed here (drives step 8.5 ROS death). The remaining entries are owned by a "
                "separate global-ODE process."
            ),
            "immune_kills": (
                "List of infected cell ids the immune layer decided to kill this record; each id still "
                "in infected_ids is switched I->D (step gated by the 'killing' enable token), replacing "
                "run.py's in-place contact/nearby-kill computation."
            ),
        },
        outputs={
            "volumes": "Inherited. Per-cell pixel volumes (list indexed by cell id) after the MCS loop.",
            "types": "Inherited. Per-cell integer type after the MCS loop (index 0 = medium sentinel).",
            "positions": "Inherited. Per-cell centre-of-mass [x,y,z] coordinates.",
            "field_at_cell": (
                "Inherited. Map str(cid) -> mean VIRUS (field 0) concentration over the cell's pixels."
            ),
            "neighbor_secretory": (
                "Inherited. Map str(cid) -> count of face-adjacent cells whose type is in secretory_types."
            ),
            "field_at_cell_all": (
                "Map str(cid) -> [mean of each of the 4 fields over the cell] (order virus, IFN, "
                "chemokine, IL-10), so a downstream immune process can read the full local environment."
            ),
            "virus_field": "Raw flat virus field (index 0) as a list, for numpy gradient sampling downstream.",
            "chemo_field": "Raw flat chemokine field (index 2) as a list.",
            "il10_field": (
                "Raw flat IL-10 field (index 3) as a list, so an immune process can read local IL-10 "
                "for self-limiting macrophage feedback (field_mean_at_cell equivalent)."
            ),
            "dims": "Lattice dimensions [nx, ny, nz].",
            "counts": "Map of epithelial population counts {H, I, D} after the loop.",
            "ode_inputs": (
                "Spatial->ODE coupling map a downstream global-ODE process needs: H, I epithelial counts, "
                "DH (dead-from-healthy), per-z-slice field integrals V/F/C/L of the virus/IFN/chemokine/"
                "IL-10 fields, and B_ei/G_ki (sum of infected-cell resistance scaled by b_ei / g_ki)."
            ),
        },
        config={
            "spec": (
                "Inherited. Declarative world spec (potts block with seed, cells, lattice). This process "
                "owns field construction, so spec must NOT carry a 'fields' key."
            ),
            "mcs_per_update": "Inherited. Historical single-knob MCS count; the default for mcs_per_step.",
            "n_fields": (
                "Inherited from CPMProcess but overridden: this process always builds exactly 4 fields "
                "(virus, IFN, chemokine, IL-10) regardless of the config value."
            ),
            "secretory_types": "Inherited. Cell types counted as secretory for the neighbor_secretory readout.",
            "mcs_per_step": (
                "Raw MCS advanced per update() call; each MCS runs world.step(1) then the full fate "
                "pipeline. Defaults to mcs_per_update."
            ),
            "enable": (
                "Set of fate tokens that fire each MCS: 'infection', 'ifn', 'death', 'ros', 'allee', "
                "'killing' (default), plus optional 'chemokine' to turn on the uninfected-H IL-10 "
                "secretion gate. Mirrors run_full_model's enable set."
            ),
            "init_viral_load": (
                "If > 0, seeds a ~uniform virus-field initial condition at construction (the fig5/fig7 "
                "'viral load' scenario, no pre-infected cells); 0.0 uses the init_infection_frac path."
            ),
        },
        math=[
            "Infection H->I:        rate = g_hv * v_bar;   Pr(infect) = 1 - exp(-rate)",
            "Resistance:            rho = f_bar / (a_rf + f_bar)",
            "Virus release scale:   s_virus(cid) = 1 - rho(cid)   (infected cells)",
            "Uninfected IL-10 gate: s_il10(cid)  = clip(1 - rho(cid), 0, 1)   (H cells, 'chemokine')",
            "Infected apoptosis I->D: rate = mu_i * (1 - rho);  Pr(die) = 1 - exp(-rate)",
            "Allee death H->D:      rate = b_h*(1-rho)*srf_D*(srf_thr - srf_u)/srf_tot^2,  srf_thr = theta*srf_tot",
            "Allee recovery D->H:   rate = b_h*(1-rho)*srf_u*(srf_u - srf_thr)/srf_tot^2",
            "ROS death:             Pr_I = 1 - exp(-g_ix*hill(X,a_ix,h_x)),  Pr_H = 1 - exp(-g_hx*hill(X,a_hx,h_x))",
            "Hill:                  hill(X,a,h) = 1 / (1 + (a/X)^h)   (0 for X<=0)",
            "ODE integrals:         V,F,C,L = sum(field_conc) / nz;   B_ei = b_ei*sum_I(rho),  G_ki = g_ki*sum_I(rho)",
        ],
        symbols={
            "v_bar": "mean local virus-field concentration over a cell's pixels (virus-field units)",
            "f_bar": "mean local IFN-field concentration over a cell's pixels (IFN-field units)",
            "rho": "per-cell viral resistance from local IFN, in [0,1) (dimensionless)",
            "g_hv": "infection rate coefficient, per MCS per unit virus concentration (1/MCS per conc)",
            "a_rf": "IFN half-saturation constant for resistance (IFN-field units)",
            "mu_i": "base infected-cell apoptosis rate (1/MCS)",
            "b_h": "Allee death/recovery rate coefficient (1/MCS, per contact-geometry factor)",
            "theta": "Allee surface-fraction threshold srf_thr/srf_tot (dimensionless)",
            "srf_u": "cell surface area in contact with uninfected neighbours (lattice sites)",
            "srf_D": "cell surface area in contact with dead neighbours (lattice sites)",
            "srf_tot": "total cell surface area (lattice sites)",
            "X": "well-mixed ROS oxidant level from ode_state (ODE concentration units)",
            "g_ix": "ROS kill-rate coefficient for infected cells (1/MCS)",
            "g_hx": "ROS kill-rate coefficient for healthy cells (1/MCS)",
            "a_ix": "ROS Hill half-max for infected-cell death (X units)",
            "a_hx": "ROS Hill half-max for healthy-cell death (X units)",
            "h_x": "ROS Hill exponent (dimensionless)",
            "b_ei": "scaling from summed infected resistance to the B_ei ODE input (dimensionless)",
            "g_ki": "scaling from summed infected resistance to the G_ki ODE input (dimensionless)",
            "nz": "lattice z-depth, the divisor turning field sums into per-slice integrals (lattice sites)",
            "V, F, C, L": "per-slice integrals of the virus, IFN, chemokine, IL-10 fields (conc * area)",
        },
        assumptions=[
            "Fates fire EVERY raw MCS (per world.step(1)), not once per update() call, matching "
            "run_full_model's _one_mcs cadence; update() does not delegate to CPMProcess.update.",
            "The four fields are built imperatively pre-finalize (virus, IFN, chemokine, IL-10 = "
            "indices 0,1,2,3); initialize() does NOT call super().initialize because World.add_field "
            "raises after finalize.",
            "Infected-cell death uses the CC3D source's simplified flat form rate = mu_i*(1-rho), with "
            "no viral-load/ROS Hill saturation term (deviates from the paper's printed Table 2 row).",
            "Resistance uses the source's saturating form rho = f_bar/(a_rf+f_bar), not the paper's "
            "printed rho = 1 - f_bar'/(theta*a_f + f_bar'); consumers apply (1-rho) protectively.",
            "Immune killing (I->D) consumes a precomputed immune_kills list rather than computing "
            "contact/nearby kills in place (the future agent-based ImmuneProcess owns that decision).",
            "Macrophage-driven dynamic chemokine/IL-10 secretion (sig_1) is out of scope; only the "
            "macrophage-independent uninfected-H IL-10 gate is implemented, under the 'chemokine' token.",
            "init_viral_load seeds a documented ~uniform virus IC (one transient H-cell source advance), "
            "an IC approximation, not a secretion-rate tune.",
        ],
        references=[
            "Sego et al. (2022), 'A modular framework for multiscale, multicellular, spatiotemporal "
            "modeling of acute primary viral infection and immune response in epithelial tissues' "
            "(CC3D ViralInfectionVTM), reproduced here natively.",
        ],
    )

    config_schema = dict(CPMProcess.config_schema)
    config_schema.update({
        # Number of raw MCS advanced (Potts + fate pipeline) per `update()`
        # call. Defaults to `mcs_per_update` (Task 1.2's knob) so existing
        # single-knob callers are unaffected; set independently to decouple
        # "how many MCS per update" from any historical `mcs_per_update` use.
        "mcs_per_step": {"_type": "integer", "_default": 1},
        "enable": {"_type": "list", "_default": list(_DEFAULT_ENABLE)},
        # init_viral_load > 0 seeds a ~uniform virus-field IC at construction
        # (the fig5/fig7 "viral load" scenario: no pre-infected cells, infection
        # emerges from the field). 0.0 = not used (init_infection_frac path).
        "init_viral_load": {"_type": "float", "_default": 0.0},
    })

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
        seed = int(spec["potts"]["seed"])
        self.world.finalize(seed)

        # init_viral_load: seed a ~uniform virus-field IC (no pre-infected
        # cells), mirroring run.py's `_seed_uniform_virus` -- every H cell is a
        # transient virus source for exactly one field advance (the engine has
        # no field-write primitive), then reset to 0 so only infected cells
        # secrete during the run proper. Documented IC approximation, NOT a
        # rate tune (deposit scales linearly with the load). Virus field = idx 0.
        init_viral_load = float(self.config.get("init_viral_load") or 0.0)
        if init_viral_load > 0.0:
            _d, _decay, _dt, _sub, _sec = fields._virus_field_params()
            # v0 is TOTAL virus per epithelial cell -> v0/cell_sites per pixel
            # (source line 92-93); depositing per-pixel `load` over-seeds
            # ~cell_sites-fold. Mirrors run.py `_seed_uniform_virus`.
            cell_sites = int(load_params()["cpm"]["cell_sites"])
            per_pixel = init_viral_load / cell_sites
            self.world.set_secretion(0, types.H, per_pixel / _dt)
            self.world.advance_fields(1)
            self.world.set_secretion(0, types.H, 0.0)

        self.mcs = int(self.config["mcs_per_update"])
        self.mcs_per_step = int(self.config.get("mcs_per_step") or self.mcs)
        # Always exactly the 4 fields added above, regardless of the
        # `n_fields` config value (this process owns field construction).
        self.n_fields = 4
        self.secretory = set(int(t) for t in self.config.get("secretory_types") or [])
        self.dims = self.world.dims()

        self.enable = set(self.config.get("enable") or _DEFAULT_ENABLE)

        # --- source-faithful constants (run.py run_full_model setup,
        # run.py:~1722-1849) resolved once here, never re-derived per MCS. ---
        cells = spec["cells"]
        self.n_cells = len(cells)
        num_epithelial = sum(1 for c in cells if c["type"] in (types.H, types.I))

        params = load_params()
        self.g_hv = float(params["virus"]["infection_g_hv"])
        self.a_rf = float(params["resistance"]["a_rf"])
        self.mu_i = float(params["cell_death"]["mu_i_per_mcs"])
        self.b_h = float(params["allee"]["b_h"])
        self.theta = float(params["allee"]["srf_threshold"])

        consts = price_ode.resolve_constants(params["price_ode"], num_epithelial=num_epithelial)
        self.b_ei = consts["b_ei"]
        self.g_ki = consts["g_ki"]
        self.g_ix = consts["g_ix"]
        self.a_ix = consts["a_ix"]
        self.g_hx = consts["g_hx"]
        self.a_hx = consts["a_hx"]
        self.h_x = int(params["price_ode"]["hill_exponents"]["h_x"])

        # --- independent RNG streams (offsets IDENTICAL to run_full_model,
        # run.py:1852-1855, so Task 1.4 parity can match seed-for-seed). ---
        self.infect_rng = np.random.default_rng(seed + 1)
        self.death_rng = np.random.default_rng(seed + 2)
        self.allee_rng = np.random.default_rng(seed + 3)
        self.ros_rng = np.random.default_rng(seed + 4)

        # --- cell-fate state, persisted across updates (closure vars in
        # run.py's run_full_model -> instance attributes here). ---
        self.infected_ids = set(
            i + 1 for i, c in enumerate(cells) if c["type"] == types.I)
        self.dead_from_healthy_ids: set[int] = set()
        self.dead_from_infected_ids: set[int] = set()
        self.resist_at_cell = [0.0] * (self.n_cells + 1)

    def inputs(self):
        base = dict(super().inputs())
        base["field_deposit"] = "list"  # each item [field_idx, x, y, z, amount]
        # ODE state pushed in each update (only "X", the ROS oxidant scalar,
        # is consumed here; a later ODE process owns/produces the rest).
        base["ode_state"] = "map[float]"
        # Infected cell ids the immune layer killed this record (I->D).
        base["immune_kills"] = "list[integer]"
        return base

    def outputs(self):
        out = dict(super().outputs())
        out["field_at_cell_all"] = "overwrite[map[list]]"
        # raw flat fields + dims so a later ImmuneProcess can sample gradients
        # in numpy (process-bigraph has no within-update request/response, so
        # the whole field crosses the store boundary once per update).
        out["chemo_field"] = "overwrite[list]"
        out["virus_field"] = "overwrite[list]"
        # IL-10 field (index 3) -- so ImmuneProcess can sample the local IL-10
        # concentration at each macrophage agent, matching run_full_model's
        # `world.field_mean_at_cell(il10_fi, cid)` self-limiting feedback
        # (source-faithfulness fix; see immune_process.py's Secretion section).
        out["il10_field"] = "overwrite[list]"
        out["dims"] = "overwrite[list]"
        out["counts"] = "overwrite[map[integer]]"
        # Controller ruling R3 (task-1.3-brief.md): the spatial->ODE inputs a
        # later GlobalODE process needs each record -- H/I epithelial counts,
        # DH (dead-from-healthy), V/F/C/L field integrals, B_ei/G_ki (summed
        # infected resistance * b_ei/g_ki). See `_ode_inputs`.
        out["ode_inputs"] = "overwrite[map[float]]"
        return out

    def update(self, state, interval):
        state = state or {}
        # External fate overrides (e.g. a subcellular process), applied once
        # per update call -- same semantics as CPMProcess.update.
        fates = state.get("fates") or {}
        for cid_key, t in fates.items():
            cid = int(cid_key)
            if t and cid > 0:
                self.world.set_cell_type(cid, int(t))

        for item in state.get("field_deposit") or []:
            fidx, x, y, z, amt = item
            self.world.field_add_source_at(int(fidx), int(x), int(y), int(z), float(amt))

        ode_state = state.get("ode_state") or {}
        immune_kills = state.get("immune_kills") or []

        # Per-MCS interleave: world.step(1) then the full fate pipeline,
        # `mcs_per_step` times -- matches run_full_model's
        # `for _ in range(mcs_per_step): _one_mcs(...)` cadence exactly (fates
        # fire EVERY MCS, not once after all MCS have been stepped).
        for _ in range(self.mcs_per_step):
            self.world.step(1)
            self._epithelial_fates(ode_state, immune_kills)

        types_now = list(self.world.cell_types())
        n = len(types_now)
        field_at = {}
        if self.n_fields > 0:
            for cid in range(1, n):
                field_at[str(cid)] = self.world.field_mean_at_cell(0, cid)
        neigh = self._neighbor_secretory_counts(types_now)
        base = {
            "volumes": list(self.world.cell_volumes()),
            "types": types_now,
            "positions": [list(c) for c in self.world.cell_coms()],
            "field_at_cell": field_at,
            "neighbor_secretory": {str(cid): neigh[cid] for cid in range(1, n)},
        }

        base["field_at_cell_all"] = {
            str(cid): [self.world.field_mean_at_cell(f, cid) for f in range(self.n_fields)]
            for cid in range(1, n)
        }
        # field indices (per fields.py order): 0=virus 1=ifn 2=chemokine 3=il10
        base["chemo_field"] = list(self.world.field_conc(2)) if self.n_fields > 2 else []
        base["virus_field"] = list(self.world.field_conc(0)) if self.n_fields > 0 else []
        base["il10_field"] = list(self.world.field_conc(3)) if self.n_fields > 3 else []
        base["dims"] = list(self.world.dims())

        H = sum(1 for t in types_now[1:] if t == types.H)
        I = sum(1 for t in types_now[1:] if t == types.I)
        D = sum(1 for t in types_now[1:] if t == types.D)
        base["counts"] = {"H": H, "I": I, "D": D}
        base["ode_inputs"] = self._ode_inputs(types_now, H, I)
        return base

    def _epithelial_fates(self, ode_state, immune_kills):
        """One source-faithful MCS of the epithelial cell-fate pipeline:
        `run.py` `_one_mcs` steps 2-8.5 (`run.py:1975`-`2097`), moved
        verbatim (exact order, exact helper calls, exact rates). The only
        structural change from the source is step 6 (contact + nearby
        killing): that in-place `kill_rng` computation is replaced by
        consuming a precomputed `immune_kills` list (the future
        agent-based ImmuneProcess's decision), merged right after infection
        and before resistance/apoptosis (task-1.3-brief.md Interfaces).
        The macrophage secretion-scale HALF of the source's "chemokine"
        block (sig_1-dependent chemokine/IL-10 secretion by macrophage
        cells) is out of scope for this task -- deferred, see module
        docstring and task-1.3-report.md. The uninfected-H IL-10 gate half
        of that same block IS implemented, gated by "chemokine" in
        `self.enable` exactly like the source (run.py:1999-2010).
        """
        world = self.world
        n_cells = self.n_cells
        current_types = list(world.cell_types())

        # (2) infection H->I (local virus * g_hv).
        if "infection" in self.enable:
            virus_at_cell = [world.field_mean_at_cell(0, cid)
                             for cid in range(n_cells + 1)]
            after_inf = transitions.infection_step(
                current_types, virus_at_cell, self.g_hv, self.infect_rng)
            for cid in range(1, n_cells + 1):
                if after_inf[cid] == types.I and current_types[cid] != types.I:
                    world.set_cell_type(cid, types.I)
                    self.infected_ids.add(cid)
            current_types = after_inf

        # Immune-layer kills I->D, merged right after infection, before
        # resistance/apoptosis (replaces run.py step 6's in-place contact +
        # nearby kill computation; gated like the source's "killing" token).
        if "killing" in self.enable:
            for cid in immune_kills:
                cid = int(cid)
                if cid in self.infected_ids:
                    world.set_cell_type(cid, types.D)
                    self.infected_ids.discard(cid)
                    self.dead_from_infected_ids.add(cid)

        # (3) per-cell resistance rho from local IFN (every epithelial cell).
        resist_at_cell = [0.0] * (n_cells + 1)
        if "ifn" in self.enable:
            for cid in range(1, n_cells + 1):
                if current_types[cid] in _EPITHELIAL_TYPES:
                    f_bar = world.field_mean_at_cell(1, cid)
                    resist_at_cell[cid] = cell_resistance(f_bar, self.a_rf)
        self.resist_at_cell = resist_at_cell

        # (4) field-secretion scales. Virus release by infected gated (1-rho),
        # UNCONDITIONAL (matches run.py: not nested under "chemokine").
        # Uninfected IL-10 gate (1-rho) IS nested under "chemokine" in the
        # source (run.py:1999-2010, same block as the macrophage secretion
        # loop -- verified by indentation) -- gated here identically. The
        # macrophage-driven chemokine/IL-10 scale itself (dynamic sig_1) is
        # still NOT implemented (out of scope, see class docstring); only
        # the macrophage-independent H-cell IL-10 gate is, under the same
        # "chemokine" token the source uses for that whole block.
        for cid in self.infected_ids:
            world.set_cell_secretion_scale(0, cid, 1.0 - resist_at_cell[cid])
        if "chemokine" in self.enable:
            for cid in range(1, n_cells + 1):
                if current_types[cid] == types.H:
                    world.set_cell_secretion_scale(
                        3, cid, signaling.uninfected_il10_scale(resist_at_cell[cid]))

        # (5) chemotaxis: persistent engine setting, applied inside world.step
        # (not wired by this task -- see module docstring).

        # (7) Allee death H->D / recovery D->H (contact-geometry rule).
        if "allee" in self.enable:
            snapshot = list(world.cell_types())
            allee_writes = []
            for cid in range(1, n_cells + 1):
                t = snapshot[cid]
                if t not in (types.H, types.D):
                    continue
                srf_uninfected, srf_D, srf_total = _epithelial_contact_totals(world, cid)
                resist = resist_at_cell[cid]
                if t == types.H:
                    rate = allee.allee_death_rate(srf_D, srf_uninfected, srf_total,
                                                  self.b_h, resist, self.theta)
                    if self.allee_rng.random() < 1.0 - math.exp(-rate):
                        allee_writes.append((cid, types.D, "death"))
                else:  # types.D
                    rate = allee.allee_recovery_rate(srf_uninfected, srf_total,
                                                     self.b_h, resist, self.theta)
                    if self.allee_rng.random() < 1.0 - math.exp(-rate):
                        allee_writes.append((cid, types.H, "recover"))
            for cid, new_t, kind in allee_writes:
                world.set_cell_type(cid, new_t)
                if kind == "death":  # H->D: dead-from-HEALTHY (source DH)
                    self.dead_from_healthy_ids.add(cid)
                else:                # D->H: recovered, no longer dead
                    self.dead_from_healthy_ids.discard(cid)
                    self.dead_from_infected_ids.discard(cid)

        # (8) infected apoptosis I->D (mu_i*(1-rho)).
        if "death" in self.enable and self.infected_ids:
            current_types = list(world.cell_types())
            after_apop = transitions.infected_death_step(
                current_types, resist_at_cell, self.mu_i, self.death_rng)
            for cid in range(1, n_cells + 1):
                if after_apop[cid] == types.D and current_types[cid] == types.I:
                    world.set_cell_type(cid, types.D)
                    self.infected_ids.discard(cid)
                    self.dead_from_infected_ids.add(cid)

        # (8.5) ROS-driven death (source OxidationAgentModelSteppable) --
        # the well-mixed oxidant X (from `ode_state`) kills infected
        # (I->D, dead-from-INFECTED) and uninfected (H->D, dead-from-HEALTHY)
        # cells this MCS.
        if "ros" in self.enable:
            X = float(ode_state.get("X", 0.0))
            if X > 0.0:
                pr_I = 1.0 - math.exp(-self.g_ix * price_ode.hill(X, self.a_ix, self.h_x))
                pr_H = 1.0 - math.exp(-self.g_hx * price_ode.hill(X, self.a_hx, self.h_x))
                if pr_I > 0.0 and self.infected_ids:
                    for cid in list(self.infected_ids):
                        if self.ros_rng.random() < pr_I:
                            world.set_cell_type(cid, types.D)
                            self.infected_ids.discard(cid)
                            self.dead_from_infected_ids.add(cid)
                if pr_H > 0.0:
                    types_ros = list(world.cell_types())
                    for cid in range(1, n_cells + 1):
                        if types_ros[cid] == types.H and self.ros_rng.random() < pr_H:
                            world.set_cell_type(cid, types.D)
                            self.dead_from_healthy_ids.add(cid)

    def _ode_inputs(self, types_now, H, I):
        """Controller ruling R3 (task-1.3-brief.md): the spatial->ODE inputs
        a later GlobalODE process needs, reusing `run.py` `_ode_couple`'s
        exact expressions (run.py:2098-2128) for DH/V/F/C/L/B_ei/G_ki."""
        DH = sum(1 for cid in self.dead_from_healthy_ids if types_now[cid] == types.D)
        dim_z = self.dims[2]
        V = float(sum(self.world.field_conc(0))) / dim_z
        F = float(sum(self.world.field_conc(1))) / dim_z
        C = float(sum(self.world.field_conc(2))) / dim_z
        L = float(sum(self.world.field_conc(3))) / dim_z
        sum_resist = sum(self.resist_at_cell[cid] for cid in self.infected_ids)
        B_ei = sum_resist * self.b_ei
        G_ki = sum_resist * self.g_ki
        return {
            "H": float(H), "I": float(I), "DH": float(DH),
            "V": V, "F": F, "C": C, "L": L,
            "B_ei": B_ei, "G_ki": G_ki,
        }
