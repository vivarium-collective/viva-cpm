"""``ProtrusionEnclosureProcess`` — the S3 enclosure model as a process-bigraph
Process (campaign composite B, the active-membrane / stochastic-encounter model).

It owns a two-cell host/symbiont Cellular Potts world (``enclosure.build_world``)
and, each process step, advances the Rust engine by ``mcs_per_update`` Monte-Carlo
sweeps and emits the wrapping observables the S3 gate (G3 Enclosure) reads:
``enc/enclosure_fraction`` (host coverage of the symbiont boundary, 1.0 = fully
enclosed), ``enc/enclosed`` (past the enclosure threshold), and the raw boundary
counts. Enclosure is therefore an *emergent, measured* outcome of the mechanics —
not a threshold imposed by hand.
"""
from __future__ import annotations

from process_bigraph import Process

from bigraph_schema.contract import ProcessContract

from .enclosure import (
    ENCLOSED_THRESH,
    HOST,
    MEDIUM,
    SYMBIONT,
    EnclosureParams,
    build_world,
    type_grid,
    wrapping_fraction,
)


class ProtrusionEnclosureProcess(Process):
    """Host/symbiont differential-adhesion wrapping on the Rust CPM engine."""

    contract = ProcessContract(
        summary=(
            "Advances a two-cell Cellular Potts world — a deformable archaeal-like HOST and a "
            "compact bacterial SYMBIONT in a shared medium — by `mcs_per_update` Monte-Carlo "
            "sweeps each step, and emits how far the host has wrapped the symbiont. Differential "
            "adhesion (Glazier-Graner Fig. 18 engulfment) drives wrapping when the host-symbiont "
            "adhesion J_hs beats the combined surface-tension cost J_sM - J_hM; enclosure is the "
            "emergent plateau of the host's coverage of the symbiont boundary, not an imposed "
            "threshold."
        ),
        description=(
            "Campaign composite B for eukaryogenesis S3 (protrusions to enclosure). The energies "
            "J_hM, J_sM, J_hs set host and symbiont surface tensions against the medium and the "
            "host-symbiont adhesion; the volume constraint sets cell sizes. Sweeping J_hs across "
            "J_sM - J_hM moves the system through the wrapping transition (stable enclosure -> "
            "partial cup -> free symbiont), the G3 discriminator."
        ),
        inputs={},
        outputs={
            "enclosure_fraction": "host coverage of the symbiont boundary (0..1; 1 = fully enclosed)",
            "enclosed": "enclosure_fraction has passed the enclosure threshold",
            "host_contact_faces": "symbiont boundary faces touching the host",
            "medium_contact_faces": "symbiont boundary faces exposed to the medium",
            "sweeps": "cumulative Monte-Carlo sweeps advanced",
        },
    )

    config_schema = {
        "nx": {"_type": "integer", "_default": 90},
        "temperature": {"_type": "float", "_default": 10.0},
        "lambda_volume": {"_type": "float", "_default": 2.0},
        "j_sM": {"_type": "float", "_default": 16.0},
        "j_hM": {"_type": "float", "_default": 4.0},
        "j_hs": {"_type": "float", "_default": 6.0},
        "symbiont_radius": {"_type": "float", "_default": 8.0},
        "host_radius": {"_type": "float", "_default": 12.0},
        "symbiont_target": {"_type": "float", "_default": 200.0},
        "host_target": {"_type": "float", "_default": 400.0},
        "seed": {"_type": "integer", "_default": 7},
        "mcs_per_update": {"_type": "integer", "_default": 60},
    }

    def initialize(self, config):
        self.params = EnclosureParams(
            nx=int(config["nx"]), temperature=float(config["temperature"]),
            lambda_volume=float(config["lambda_volume"]), j_sM=float(config["j_sM"]),
            j_hM=float(config["j_hM"]), j_hs=float(config["j_hs"]),
            symbiont_radius=float(config["symbiont_radius"]),
            host_radius=float(config["host_radius"]),
            symbiont_target=float(config["symbiont_target"]),
            host_target=float(config["host_target"]), seed=int(config["seed"]))
        self.nx = self.params.nx
        self.mcs = int(config["mcs_per_update"])
        self.world = build_world(self.params)
        self._sweeps = 0

    def inputs(self):
        return {}

    def outputs(self):
        # overwrite semantics: each step's reading REPLACES the store (a bare scalar
        # apply accumulates/adds, which would sum the per-step fractions).
        return {
            "enclosure_fraction": "overwrite[float]",
            "enclosed": "overwrite[boolean]",
            "host_contact_faces": "overwrite[integer]",
            "medium_contact_faces": "overwrite[integer]",
            "sweeps": "overwrite[integer]",
        }

    def update(self, state, interval):
        self.world.step_parallel(self.mcs, 16)
        self._sweeps += self.mcs
        g = type_grid(self.world, self.nx)
        sym = g == SYMBIONT
        import numpy as np
        host_faces = med_faces = 0
        for ax, sh in [(0, 1), (0, -1), (1, 1), (1, -1)]:
            ng = np.roll(g, sh, axis=ax)
            host_faces += int(np.sum(sym & (ng == HOST)))
            med_faces += int(np.sum(sym & (ng == MEDIUM)))
        frac = host_faces / (host_faces + med_faces) if (host_faces + med_faces) else 0.0
        return {
            "enclosure_fraction": float(frac),
            "enclosed": bool(frac >= ENCLOSED_THRESH),
            "host_contact_faces": int(host_faces),
            "medium_contact_faces": int(med_faces),
            "sweeps": int(self._sweeps),
        }
