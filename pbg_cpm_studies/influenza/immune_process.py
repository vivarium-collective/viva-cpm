"""``ImmuneProcess`` -- a process-bigraph `Process` holding OFF-LATTICE immune
agents (macrophages, NK cells, CD8+ T cells) that chemotax up a diffusive
field gradient, kill nearby infected epithelial cells (NK/CD8), and secrete
chemokine/IL-10 (macrophages).

Task 2.1 shipped agent state + chemotaxis only. Task 2.2 added the other two
off-lattice mechanisms (killing, secretion). This revision (Task 2.3) adds
ODE-DRIVEN RECRUITMENT (spawning/removing agents) -- see "Recruitment" below.

Why the process re-samples the whole field itself (rather than requesting a
point sample from the CPM engine, as `EpitheliumProcess`'s `field_at_cell`
does): process-bigraph has no within-update request/response channel, so an
off-lattice agent can't ask the lattice "what's the gradient at my
position?" mid-step. Instead the WHOLE flat field crosses the store boundary
each update (`chemo_field`, `virus_field`, `dims`) and this process samples
a local central-difference gradient itself, in numpy, at each agent's
nearest lattice site.

## Killing (Task 2.2)

`killing.contact_kill_rate` is the on-lattice LOCAL kill-rate term (see that
module's docstring / `run.py`'s `_apply_killing`):

    rate = g_i * tot_ec * srf_immune * cell_resist / cell_volume

`srf_immune` (CPM contact-surface area, lattice sites) and `cell_volume`
(the infected cell's CPM volume) are both CPM-geometry quantities this
off-lattice process has no access to (agents aren't wired to a `World`).
Per the brief, PROXIMITY (agent within `kill_radius` of the infected cell's
centroid) REPLACES the contact-area gate: an agent within range contributes
exactly a UNIT contact area (`srf_immune = 1.0`) -- "in range" is a boolean
gate, not a graded overlap measure. `cell_volume` is likewise fixed at
`1.0` (a unit-volume proxy) and `cell_resist` at `1.0` (neutral/no local-IFN
resistance, since this process doesn't yet see per-epithelial-cell IFN
readings) -- i.e. everything except `g_i` (`nk.g_ik`/`cd8.g_ie`, unchanged
source constants) and `tot_ec` (`scaling.ode_epithelial_population`,
unchanged) is held at its identity/neutral value. This reduces
`contact_kill_rate` to exactly params.yaml's already-recorded
`pr_cf_nk_loc`/`pr_cf_cd8_loc` constants (`g_i * tot_ec`) -- a deliberate,
non-fabricated choice (no new number is invented; the "proxy" is the
identity element of the existing formula), NOT a fit to any target outcome.
`kill_radius` defaults to `sqrt(cpm.cell_sites)` (~5.0 lattice sites) -- the
side length of a single CPM epithelial cell's footprint (`cell_sites` = 25
sites, a 5x5 square) -- i.e. "adjacent" means within about one cell-width,
NOT a fitted/tuned distance.

Per-agent-per-infected-cell kill probability is `1 - exp(-rate)`, drawn once
per K/E agent x in-range infected-cell pair per `update()` call (one MCS);
an infected cell already killed earlier in the SAME update (by an earlier
agent) is skipped for the rest of that update (killed at most once/update).

## Secretion (Task 2.2)

Each macrophage (`types.M`) agent deposits chemokine (field idx 2) and IL-10
(field idx 3) point sources at its own (x, y) via the `field_deposit`
output (consumed by `EpitheliumProcess.field_add_source_at`, which truncates
to integer lattice coordinates):

    amount_chemo = b_c_per_site * cell_sites * scale
    amount_il10  = b_l_per_site * cell_sites * scale
    scale = signaling.macrophage_secretion_scale(L_local, sig_1, g_1, g_2, d_2)

`b_c_per_site`/`b_l_per_site` derive `params.yaml`'s `chemokine.b_c`/
`il10.b_l` (the CC3D-source per-CELL, dim.z=2-baked coefficients) down to a
per-SITE amount the SAME way `fields._per_pixel_secretion_rate` does for the
continuous-secretion path -- halve for dim.z (z=1 field here) and spread
over `cpm.cell_sites` pixels -- MINUS that helper's `/dt` term, which exists
only to cancel `world.step`'s internal `rate*dt`-per-MCS re-multiplication
in the CONTINUOUS `world.set_secretion` path; a `field_deposit` point source
is added once, undivided, via `field_add_source_at`, so no `dt` term
applies here. `g_1`, `g_2`, `d_2` come from `fields.il10_hill_constants()`
(unchanged); `sig_1` is NOT the stub in that tuple -- it is read fresh each
update from the `sig_1` input (the dynamic ODE value, per the brief). Local
IL-10 (`L_local`) is fixed at `0.0`: this process does not yet receive an
`il10_field` input (out of Task 2.2 scope per the brief -- macrophage agents
are themselves an IL-10 SOURCE, not yet a consumer of the diffused field
here), so the Hill self-regulation term sees no local IL-10 buildup.

**Whole-cell flux, not one pixel (controller ruling R4, fix round 1):** an
off-lattice macrophage AGENT represents one WHOLE CPM cell. On-lattice
(`run_full_model`/`fields.add_chemokine_field`/`add_il10_field`), that whole
cell secretes its per-site rate (`b_c_per_site`/`b_l_per_site`) at EVERY one
of its `cell_sites` (~25) occupied pixels each MCS -- its total per-MCS mass
release is `cell_sites` times the single-pixel rate, not the single-pixel
rate itself. Depositing only `b_c_per_site * scale` at the agent's one point
would understate a macrophage's true secreted mass by ~25x, weakening the
chemokine gradient the NK/CD8 recruitment loop chemotaxes up -- so both
deposited amounts are multiplied by `cell_sites` here. This is
source-faithful (it matches the on-lattice per-macrophage total mass), NOT a
tuning change -- no constant changes value, only the point-vs-whole-cell
accounting is corrected.

## Recruitment (Task 2.3)

`run_full_model` (`run.py`'s `_recruit_step`) activates new CPM cells from a
fixed-size RESERVE POOL per type (`pool_t.pop(0)`): a Poisson-drawn inflow
count in excess of the remaining pool simply fails to seed (a hard cap). Off-
lattice agents here don't occupy CPM lattice sites, so there is no pool/seat
limit to draw from -- every `recruitment.poisson_inflow_count` draw becomes a
newly minted agent. This is the intended difference from `run_full_model`
and the core fix this refactor exists to make: it lets the immune
population actually reach the ODE's predicted equilibrium instead of
saturating a small fixed pool.

Per type (M/K/E), each `update()` call for which `apply_recruitment` is
truthy:
  - **Outflow (remove):** each existing agent of that type is removed with
    probability `recruitment.ul_rate_to_prob(local_ratio * outflow)`
    (`state["recruit_drivers"]`'s `*_outflow` rate).
  - **Inflow (spawn):** `recruitment.poisson_inflow_count(local_ratio *
    inflow, self.recruit_rng)` new agents of that type are minted at
    uniform-random positions inside `state["margin_box"]` (`[x0, y0, x1,
    y1]`), each with a fresh id from a monotonic counter (`self.
    _next_agent_id`) held on the process. No cap.
  - Outflow is applied before inflow for a given type, per-type in M, K, E
    order -- matching `_recruit_step`'s own within-type and across-type
    ordering (so the two implementations' RNG-consumption *shape* agrees,
    even though the two processes draw from independent RNG streams and are
    not expected to reproduce identical outcomes).

**`local_ratio` = 1.0 for all three types (design decision, not a rate
tuning):** `params["coupling"]["recruitment"]["local_ratios"]` (macro=1.0,
nk=0.75, cd8=0.75) splits each type's ODE inflow into a LOCAL fraction
(seeded as an actual CPM cell) and a NEARBY fraction (accrued into a
well-mixed ODE surrogate population, `M_nb`/`K_nb`/`E_nb` in `run.py`). That
split exists *because* `run_full_model`'s local fraction is capped by a
small reserve pool -- the nearby surrogate is a bookkeeping device so
inflow the pool has no room for isn't silently discarded. Uncapped agents
remove the reason for that split entirely, so per the spec ("with agents
uncapped, local_ratio -> 1.0") this process ignores the nk=0.75/cd8=0.75
entries and uses `local_ratio = 1.0` for M, K, and E alike, sending the
FULL ODE-predicted inflow/outflow to spatial agents with no nearby
surrogate. This changes how much of a fixed rate becomes a spatial agent
vs. an ODE bookkeeping number; it does NOT change any inflow/outflow RATE
constant or formula (`recruitment.macrophage_inflow` etc. are called
unmodified, on drivers computed upstream and passed in via
`recruit_drivers`).

`self.recruit_rng = np.random.default_rng(seed + 900)` is a dedicated RNG
stream (offset convention matches `run_full_model`'s `recruit_rng`) so
recruitment draws never share/correlate with the movement-noise (`self.rng`)
or kill-draw (`self.kill_rng`) streams. `self._next_agent_id` is lazily
seeded on the FIRST `update()` call, above the max id present in that
call's `immune_agents` (so ids picked up from upstream composite-provided
initial agents are never collided with), then increments monotonically
across all subsequent calls on this same process instance.
"""
from __future__ import annotations

import math

import numpy as np
from process_bigraph import Process

from . import recruitment, types
from .fields import il10_hill_constants
from .killing import contact_kill_rate
from .params import load_params
from .signaling import macrophage_secretion_scale

# Field routing: which field each agent type chemotaxes toward. Macrophages
# home in on the VIRUS field (the infection itself); NK and CD8+ T cells
# home in on the CHEMOKINE field (macrophage-secreted recruitment signal) --
# same routing as the on-lattice `immune.set_macrophage_chemotaxis` /
# `immune.set_nk_cd8_chemotaxis` engine wiring, reproduced here for the
# off-lattice agent layer.
_VIRUS_FIELD_TYPES = {types.M}
_CHEMO_FIELD_TYPES = {types.K, types.E}

# Cytotoxic agent types (proximity killing, Task 2.2). Coincides with
# `_CHEMO_FIELD_TYPES` above (NK/CD8 both chemotax toward chemokine AND both
# kill) -- kept as a separate name since the two sets mean different things.
_CYTOTOXIC_TYPES = {types.K, types.E}

# `kill_radius` default: the side length of one CPM epithelial cell's
# footprint (`cpm.cell_sites` = 25 sites, a 5x5 square) -- see this module's
# docstring "Killing" section for why this is a documented CPM-cell-footprint
# derivation, not a fitted outcome.
_params_at_import = load_params()
_CELL_SITES = int(_params_at_import["cpm"]["cell_sites"])
_DEFAULT_KILL_RADIUS = math.sqrt(_CELL_SITES)


def _gradient_at(field2d: np.ndarray, ix: int, iy: int) -> tuple[float, float]:
    """Central-difference gradient of ``field2d`` (shape (ny, nx)) at the
    integer site ``(ix, iy)``, clamped at the domain edges (falls back to a
    one-sided difference there instead of reading out of bounds)."""
    ny, nx = field2d.shape
    ix_lo, ix_hi = max(ix - 1, 0), min(ix + 1, nx - 1)
    iy_lo, iy_hi = max(iy - 1, 0), min(iy + 1, ny - 1)
    dx = ix_hi - ix_lo
    dy = iy_hi - iy_lo
    gx = (field2d[iy, ix_hi] - field2d[iy, ix_lo]) / dx if dx > 0 else 0.0
    gy = (field2d[iy_hi, ix] - field2d[iy_lo, ix]) / dy if dy > 0 else 0.0
    return float(gx), float(gy)


class ImmuneProcess(Process):
    """Off-lattice immune agents that chemotax up a field gradient.

    Agents are plain dicts ``{"id": int, "type": int, "x": float, "y": float}``
    held in the ``immune_agents`` store. Each update samples the local
    gradient of the RELEVANT field (per ``types``-based routing above) at the
    agent's nearest lattice site and moves the agent ``step_len`` sites
    up-gradient, plus a small random component.
    """

    config_schema = {
        "seed": {"_type": "integer", "_default": 17},
        "step_len": {"_type": "float", "_default": 1.0},
        "kill_radius": {"_type": "float", "_default": _DEFAULT_KILL_RADIUS},
    }

    def initialize(self, config):
        self.step_len = float(config["step_len"])
        self.rng = np.random.default_rng(int(config["seed"]))
        # Dedicated kill-draw RNG stream (offset convention matches run.py's
        # `kill_rng = np.random.default_rng(seed + 500)`) so kill rolls never
        # share/correlate with the movement-noise draws above.
        self.kill_rng = np.random.default_rng(int(config["seed"]) + 500)
        self.kill_radius = float(config["kill_radius"])
        # Dedicated recruitment RNG stream (Task 2.3; offset convention
        # matches run_full_model's `recruit_rng = np.random.default_rng(seed
        # + 900)`) -- see module docstring "Recruitment" section.
        self.recruit_rng = np.random.default_rng(int(config["seed"]) + 900)
        # Monotonic new-agent-id counter; lazily seeded above the max id seen
        # on the FIRST update() call (see module docstring).
        self._next_agent_id: int | None = None

        params = load_params()
        # Killing-rate coefficients (source constants, unchanged) -- see this
        # module's docstring "Killing" section for the proximity/unit-proxy
        # convention `contact_kill_rate` is called with below.
        self.g_ik = float(params["nk"]["g_ik"])
        self.g_ie = float(params["cd8"]["g_ie"])
        self.tot_ec = float(params["scaling"]["ode_epithelial_population"])

        # Secretion: IL-10-Hill constants (g_1/g_2/d_2; sig_1 here is the
        # tuple's STUB and is NOT used -- the dynamic `sig_1` input is used
        # instead each update, see this module's docstring "Secretion"
        # section) plus the per-site chemokine/IL-10 base rates.
        _, self.g_1, self.g_2, self.d_2 = il10_hill_constants()
        self.cell_sites = int(params["cpm"]["cell_sites"])
        self.b_c_per_site = (float(params["chemokine"]["b_c"]) / 2.0) / self.cell_sites
        self.b_l_per_site = (float(params["il10"]["b_l"]) / 2.0) / self.cell_sites

    def inputs(self):
        return {
            "chemo_field": "list",
            "virus_field": "list",
            "dims": "list",
            "immune_agents": "list",
            # Task 2.2 inputs (killing + secretion):
            "epithelial_positions": "map[list]",  # cell-id-str -> [x, y]
            "infected_ids": "list[integer]",
            "sig_1": "float",
            # Task 2.3 inputs (recruitment):
            "recruit_drivers": "map[float]",  # macro/nk/cd8 in/outflow rates
            "margin_box": "list",  # [x0, y0, x1, y1] spawn region
            "apply_recruitment": "boolean",  # cadence gate set by Composite
        }

    def outputs(self):
        return {
            "immune_agents": "overwrite[list]",
            # Infected cell ids killed THIS update (Task 2.2).
            "immune_kills": "overwrite[list[integer]]",
            # [field_idx, x, y, z, amount] point sources (Task 2.2).
            "field_deposit": "overwrite[list]",
        }

    def update(self, state, interval):
        state = state or {}
        agents = list(state.get("immune_agents") or [])
        epithelial_positions = state.get("epithelial_positions") or {}
        infected_ids = list(state.get("infected_ids") or [])
        sig_1 = float(state.get("sig_1", 0.0))

        # Task 2.3: lazily seed the new-agent-id counter above the max id
        # present on the FIRST update() call this process instance sees (see
        # module docstring "Recruitment" section) -- before any early return,
        # so recruitment can spawn even from a zero-agent start.
        if self._next_agent_id is None:
            existing_ids = [int(a.get("id", 0)) for a in agents]
            self._next_agent_id = (max(existing_ids) + 1) if existing_ids else 1

        if not agents:
            updated = self._apply_recruitment([], state)
            return {"immune_agents": updated, "immune_kills": [], "field_deposit": []}

        dims = state.get("dims") or [0, 0, 1]
        nx, ny = int(dims[0]), int(dims[1])
        if nx <= 0 or ny <= 0:
            updated = self._apply_recruitment(agents, state)
            return {"immune_agents": updated, "immune_kills": [], "field_deposit": []}

        chemo = np.asarray(state.get("chemo_field") or [], dtype=float).reshape(ny, nx)
        virus = np.asarray(state.get("virus_field") or [], dtype=float).reshape(ny, nx)

        updated = []
        killed: set[int] = set()
        field_deposit: list = []

        for agent in agents:
            agent_type = int(agent.get("type", 0))
            field2d = virus if agent_type in _VIRUS_FIELD_TYPES else chemo

            x, y = float(agent["x"]), float(agent["y"])
            ix = int(round(min(max(x, 0.0), nx - 1)))
            iy = int(round(min(max(y, 0.0), ny - 1)))

            gx, gy = _gradient_at(field2d, ix, iy)
            grad = np.array([gx, gy])
            norm = float(np.linalg.norm(grad))
            direction = grad / norm if norm > 1e-12 else np.zeros(2)

            noise = self.rng.normal(scale=0.1 * self.step_len, size=2)
            move = self.step_len * direction + noise

            new_x = min(max(x + float(move[0]), 0.0), nx - 1)
            new_y = min(max(y + float(move[1]), 0.0), ny - 1)

            updated.append({**agent, "x": new_x, "y": new_y})

            # --- Task 2.2: proximity killing (K/E agents only) ---
            if agent_type in _CYTOTOXIC_TYPES:
                g_i = self.g_ik if agent_type == int(types.K) else self.g_ie
                for cid in infected_ids:
                    if cid in killed:
                        continue
                    pos = epithelial_positions.get(str(cid))
                    if pos is None:
                        continue
                    dist = math.hypot(x - float(pos[0]), y - float(pos[1]))
                    if dist > self.kill_radius:
                        continue
                    rate = contact_kill_rate(
                        srf_immune=1.0, cell_resist=1.0, g_i=g_i,
                        tot_ec=self.tot_ec, cell_volume=1.0)
                    if self.kill_rng.random() < 1.0 - math.exp(-rate):
                        killed.add(cid)

            # --- Task 2.2: macrophage chemokine/IL-10 secretion ---
            if agent_type == int(types.M):
                il10_local = 0.0  # no il10_field input yet, see module docstring
                scale = macrophage_secretion_scale(
                    il10_local, sig_1, self.g_1, self.g_2, self.d_2)
                # Whole-cell flux (controller ruling R4, fix round 1): see this
                # module's docstring "Secretion" section, final paragraph.
                amount_chemo = self.b_c_per_site * self.cell_sites * scale
                amount_il10 = self.b_l_per_site * self.cell_sites * scale
                if amount_chemo > 0.0:
                    field_deposit.append([2, x, y, 0, amount_chemo])
                if amount_il10 > 0.0:
                    field_deposit.append([3, x, y, 0, amount_il10])

        # --- Task 2.3: ODE-driven recruitment (uncapped spawn/remove) ---
        updated = self._apply_recruitment(updated, state)

        return {
            "immune_agents": updated,
            "immune_kills": sorted(killed),
            "field_deposit": field_deposit,
        }

    def _apply_recruitment(self, agents: list, state: dict) -> list:
        """Task 2.3: mint/remove immune agents directly, UNCAPPED, driven by
        the ODE's per-MCS inflow/outflow rates. No-op unless
        ``state["apply_recruitment"]`` is truthy (the Composite sets this so
        recruitment runs at the right cadence, see Task 3.3) or
        ``state["margin_box"]`` is missing. See module docstring
        "Recruitment" section for the full rationale, including why
        ``local_ratio`` is 1.0 for every type here."""
        if not state.get("apply_recruitment"):
            return agents
        margin_box = state.get("margin_box")
        if not margin_box:
            return agents

        drivers = state.get("recruit_drivers") or {}
        x0, y0, x1, y1 = (float(v) for v in margin_box)
        lo_x, hi_x = min(x0, x1), max(x0, x1)
        lo_y, hi_y = min(y0, y1), max(y0, y1)

        # local_ratio=1.0 for all types -- see module docstring "Recruitment"
        # section (uncapped agents remove the reason for the
        # local/nearby-surrogate split; this is a design decision, not a
        # tuning change to any rate).
        local_ratio = 1.0

        # (type, inflow_rate, outflow_rate) in M, K, E order -- matches
        # run.py's `_recruit_step` per-type ordering (see module docstring).
        by_type = [
            (types.M, float(drivers.get("macro_inflow", 0.0)),
             float(drivers.get("macro_outflow", 0.0))),
            (types.K, float(drivers.get("nk_inflow", 0.0)),
             float(drivers.get("nk_outflow", 0.0))),
            (types.E, float(drivers.get("cd8_inflow", 0.0)),
             float(drivers.get("cd8_outflow", 0.0))),
        ]

        updated = list(agents)
        for agent_type, inflow, outflow in by_type:
            t = int(agent_type)

            # --- Outflow: remove each existing agent of this type w.p. ---
            pr_out = recruitment.ul_rate_to_prob(local_ratio * outflow)
            if pr_out > 0.0:
                updated = [
                    a for a in updated
                    if int(a.get("type", 0)) != t or self.recruit_rng.random() >= pr_out
                ]

            # --- Inflow: spawn fresh agents at random margin-box positions ---
            n_new = recruitment.poisson_inflow_count(local_ratio * inflow, self.recruit_rng)
            for _ in range(n_new):
                new_x = float(self.recruit_rng.uniform(lo_x, hi_x))
                new_y = float(self.recruit_rng.uniform(lo_y, hi_y))
                updated.append({"id": self._next_agent_id, "type": t, "x": new_x, "y": new_y})
                self._next_agent_id += 1

        return updated
