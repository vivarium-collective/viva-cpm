"""``ImmuneProcess`` -- a process-bigraph `Process` holding OFF-LATTICE immune
agents (macrophages, NK cells, CD8+ T cells) that chemotax up a diffusive
field gradient.

This is the FIRST slice (Task 2.1) of the swappable immune layer: agent
state + chemotaxis only. Killing, secretion, and recruitment are later
increments (Tasks 2.2/2.3) and are deliberately NOT implemented here.

Why the process re-samples the whole field itself (rather than requesting a
point sample from the CPM engine, as `EpitheliumProcess`'s `field_at_cell`
does): process-bigraph has no within-update request/response channel, so an
off-lattice agent can't ask the lattice "what's the gradient at my
position?" mid-step. Instead the WHOLE flat field crosses the store boundary
each update (`chemo_field`, `virus_field`, `dims`) and this process samples
a local central-difference gradient itself, in numpy, at each agent's
nearest lattice site.
"""
from __future__ import annotations

import numpy as np
from process_bigraph import Process

from . import types

# Field routing: which field each agent type chemotaxes toward. Macrophages
# home in on the VIRUS field (the infection itself); NK and CD8+ T cells
# home in on the CHEMOKINE field (macrophage-secreted recruitment signal) --
# same routing as the on-lattice `immune.set_macrophage_chemotaxis` /
# `immune.set_nk_cd8_chemotaxis` engine wiring, reproduced here for the
# off-lattice agent layer.
_VIRUS_FIELD_TYPES = {types.M}
_CHEMO_FIELD_TYPES = {types.K, types.E}


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
    }

    def initialize(self, config):
        self.step_len = float(config["step_len"])
        self.rng = np.random.default_rng(int(config["seed"]))

    def inputs(self):
        return {
            "chemo_field": "list",
            "virus_field": "list",
            "dims": "list",
            "immune_agents": "list",
        }

    def outputs(self):
        return {"immune_agents": "overwrite[list]"}

    def update(self, state, interval):
        state = state or {}
        agents = list(state.get("immune_agents") or [])
        if not agents:
            return {"immune_agents": []}

        dims = state.get("dims") or [0, 0, 1]
        nx, ny = int(dims[0]), int(dims[1])
        if nx <= 0 or ny <= 0:
            return {"immune_agents": agents}

        chemo = np.asarray(state.get("chemo_field") or [], dtype=float).reshape(ny, nx)
        virus = np.asarray(state.get("virus_field") or [], dtype=float).reshape(ny, nx)

        updated = []
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

        return {"immune_agents": updated}
