"""``ImmuneProcess`` -- a process-bigraph `Process` holding OFF-LATTICE immune
agents (macrophages, NK cells, CD8+ T cells) that chemotax up a diffusive
field gradient, kill nearby infected epithelial cells (NK/CD8), and secrete
chemokine/IL-10 (macrophages).

Task 2.1 shipped agent state + chemotaxis only. This revision (Task 2.2)
adds the other two off-lattice mechanisms; agent RECRUITMENT (spawning new
agents) is still deferred to Task 2.3 and is NOT implemented here.

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

    amount_chemo = b_c_per_site * scale
    amount_il10  = b_l_per_site * scale
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
"""
from __future__ import annotations

import math

import numpy as np
from process_bigraph import Process

from . import types
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
        cell_sites = int(params["cpm"]["cell_sites"])
        self.b_c_per_site = (float(params["chemokine"]["b_c"]) / 2.0) / cell_sites
        self.b_l_per_site = (float(params["il10"]["b_l"]) / 2.0) / cell_sites

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

        if not agents:
            return {"immune_agents": [], "immune_kills": [], "field_deposit": []}

        dims = state.get("dims") or [0, 0, 1]
        nx, ny = int(dims[0]), int(dims[1])
        if nx <= 0 or ny <= 0:
            return {"immune_agents": agents, "immune_kills": [], "field_deposit": []}

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
                amount_chemo = self.b_c_per_site * scale
                amount_il10 = self.b_l_per_site * scale
                if amount_chemo > 0.0:
                    field_deposit.append([2, x, y, 0, amount_chemo])
                if amount_il10 > 0.0:
                    field_deposit.append([3, x, y, 0, amount_il10])

        return {
            "immune_agents": updated,
            "immune_kills": sorted(killed),
            "field_deposit": field_deposit,
        }
