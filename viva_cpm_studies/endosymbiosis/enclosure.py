"""S3 — protrusions to enclosure: the host/symbiont wrapping model on the real
Rust Cellular Potts engine.

A deformable HOST cell and a compact SYMBIONT cell share a medium, with actual
spatial boundaries. Differential adhesion (Glazier & Graner 1993, Fig. 18
engulfment) drives the host to wrap the symbiont when the energy gained by
host-symbiont contact beats the cost of host membrane deformation plus the lost
host-medium and symbiont-medium interface:

    engulfment when  gamma(symbiont, medium) > gamma(host, medium) + gamma(host, symbiont)
    with self-energies ~0 this is       J_sM  >  J_hM + J_hs

so the wrapping transition in the host-symbiont adhesion J_hs sits at
``J_sM - J_hM``.  Lower J_hs = stronger adhesion = enclosure.

This module is the shared kernel used by both the process-bigraph composite
(``ProtrusionEnclosureProcess``) and the phase-diagram driver.  The Hamiltonian,
Metropolis dynamics and lattice live in the compiled ``cpm`` engine; here we only
build the two-cell world and read the wrapping fraction back out.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from cpm import cpm_core

MEDIUM, SYMBIONT, HOST = 0, 1, 2  # medium is always owner/type 0


@dataclass
class EnclosureParams:
    nx: int = 90
    temperature: float = 10.0
    lambda_volume: float = 2.0
    j_sM: float = 16.0            # symbiont-medium surface tension (high -> wants covering)
    j_hM: float = 4.0            # host-medium surface tension (low -> fluid, spreads)
    j_hs: float = 6.0            # host-symbiont adhesion (SWEEP; low -> strong adhesion)
    symbiont_radius: float = 8.0
    host_radius: float = 12.0
    symbiont_target: float = 200.0
    host_target: float = 400.0
    seed: int = 7

    @property
    def wrapping_threshold(self) -> float:
        """Predicted J_hs below which the host engulfs the symbiont."""
        return self.j_sM - self.j_hM


def _disc(nx, cx, cy, r):
    yy, xx = np.mgrid[0:nx, 0:nx]
    return ((xx - cx) ** 2 + (yy - cy) ** 2) <= r * r


def build_world(p: EnclosureParams):
    """Two cells in contact: symbiont disc above, host disc just below it."""
    w = cpm_core.World((p.nx, p.nx, 1), "noflux", 2, p.temperature)
    for (a, b), j in {
        (MEDIUM, MEDIUM): 0.0,
        (HOST, MEDIUM): p.j_hM,
        (SYMBIONT, MEDIUM): p.j_sM,
        (HOST, SYMBIONT): p.j_hs,
        (HOST, HOST): 0.0,
        (SYMBIONT, SYMBIONT): 0.0,
    }.items():
        w.set_contact(a, b, float(j))

    labels = np.zeros((p.nx, p.nx), dtype=np.uint32)
    cx = p.nx // 2
    labels[_disc(p.nx, cx, int(p.nx * 0.42), p.symbiont_radius)] = 1
    host = _disc(p.nx, cx, int(p.nx * 0.62), p.host_radius) & (labels == 0)
    labels[host] = 2
    types_map = {1: SYMBIONT, 2: HOST}
    label_to_cid = w.seed_from_labels(
        labels.ravel().tolist(), types_map, HOST, p.symbiont_target, p.lambda_volume)
    targets = {SYMBIONT: p.symbiont_target, HOST: p.host_target}
    for label, cid in label_to_cid.items():
        w.set_target_volume(int(cid), targets[types_map[int(label)]])
    w.finalize(p.seed)
    return w


def type_grid(w, nx):
    owners = np.asarray(w.snapshot(), dtype=np.int64).reshape(nx, nx)
    types_of = np.asarray(w.cell_types(), dtype=np.int64)
    return types_of[owners]


def wrapping_fraction(g) -> float:
    """Fraction of the symbiont's exposed boundary that the host covers.

    1.0 = fully enclosed (no symbiont-medium interface left); 0.0 = free symbiont.
    """
    sym = g == SYMBIONT
    host_faces = med_faces = 0
    for ax, sh in [(0, 1), (0, -1), (1, 1), (1, -1)]:
        ng = np.roll(g, sh, axis=ax)
        host_faces += int(np.sum(sym & (ng == HOST)))
        med_faces += int(np.sum(sym & (ng == MEDIUM)))
    tot = host_faces + med_faces
    return host_faces / tot if tot else 0.0


@dataclass
class Trajectory:
    sweeps: np.ndarray
    wrapping: np.ndarray
    final_grid: np.ndarray
    params: EnclosureParams = field(default=None)


def run_trajectory(p: EnclosureParams, total_sweeps: int = 2400, block: int = 60) -> Trajectory:
    """Run one realization, sampling the wrapping fraction every ``block`` sweeps."""
    w = build_world(p)
    xs = [0]
    ys = [wrapping_fraction(type_grid(w, p.nx))]
    for s in range(block, total_sweeps + 1, block):
        w.step_parallel(block, 16)
        xs.append(s)
        ys.append(wrapping_fraction(type_grid(w, p.nx)))
    return Trajectory(np.array(xs), np.array(ys), type_grid(w, p.nx), p)


# Enclosure-regime classification for the phase diagram. A regime is read off the
# wrapping trajectory: the stable plateau (last fifth mean) and whether it held or
# decayed from an earlier peak (transient wrapping that did not stick).
ENCLOSED_THRESH = 0.85
PARTIAL_THRESH = 0.45
TRANSIENT_DROP = 0.20


def classify_regime(traj: Trajectory) -> str:
    y = traj.wrapping
    plateau = float(np.mean(y[-max(1, len(y) // 5):]))
    peak = float(np.max(y))
    if plateau >= ENCLOSED_THRESH:
        return "stable enclosure"
    if peak >= ENCLOSED_THRESH and (peak - plateau) >= TRANSIENT_DROP:
        return "transient internalization"
    if plateau >= PARTIAL_THRESH:
        return "partial wrapping"
    return "no enclosure"
