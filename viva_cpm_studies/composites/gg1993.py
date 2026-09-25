"""process-bigraph composite factories for the 11 Glazier & Graner (1993)
simulation examples.

Each factory returns a composite *document* wrapping the real Rust CPM engine
(`cpm.processes.cpm_process.CPMProcess`) configured with that study's exact
energies, temperature, area constraint and initial condition.  These are the
pbg-native, registry-resolvable encodings referenced by each study's
`baseline[].composite`.  The published figures are produced by the study driver
(`viva_cpm_studies.gg1993.driver`) at full scale; these composites run a modest
aggregate so the dashboard's "run baseline" is fast and live.

Factory names match the study slugs so the study `composite:` field reads
`viva_cpm_studies.composites.gg1993.<slug>`.
"""

from __future__ import annotations

import numpy as np

from process_bigraph.composite_generator import composite_generator

from ..gg1993.params import STUDIES
from ..gg1993 import engine, run as runmod
from ..gg1993.types import MEDIUM, DARK, LIGHT

CPM_ADDR = "local:!viva_cpm.processes.cpm_process.CPMProcess"

# modest live-demo aggregate (full-scale runs come from the driver)
_DEMO_NX = 90
_DEMO_RECT = (12, 12, 78, 78)


def _type_assignment(slug, owner):
    s = STUDIES[slug]
    ic = s["ic"]
    if ic in ("equilibrated_light", "brick_equilibrate"):
        return runmod.assign_all_light(owner)
    if ic == "half_split":
        return runmod.assign_half_split(owner, light_on_top=True)
    return runmod.assign_random(owner, s["frac_light"], seed=runmod.__dict__.get("SEED", 0) or 0)


def build_spec(slug, nx=_DEMO_NX, rect=_DEMO_RECT):
    """Return a load_world spec dict for a modest live demo of `slug`."""
    s = STUDIES[slug]
    labels = engine.brick_tiling_labels(nx, nx, rect, 5, 8)
    owner = labels.reshape(nx, nx)
    types_map = _type_assignment(slug, owner)
    en = engine.energies_from_paper(**s["J"])
    contact = [{"a": int(a), "b": int(b), "j": float(j)} for (a, b), j in en.items()]
    tgt = s["target"]
    # per-label target volume via seed_labels default + explicit cells is not
    # available; encode each cell explicitly so per-type targets apply.
    ny = nx
    cells = []
    for l in np.unique(labels):
        if l == 0:
            continue
        ys, xs = np.where(owner == l)
        t = int(types_map.get(int(l), LIGHT))
        cells.append({
            "type": t,
            "target_volume": float(tgt[t]),
            "lambda_volume": float(s["lam"]),
            "target_surface": 0.0, "lambda_surface": 0.0,
            "seed_block": [int(xs.min()), int(ys.min()), 0,
                           int(xs.max()) + 1, int(ys.max()) + 1, 1],
        })
    return {
        "potts": {"dims": [nx, ny, 1], "boundary": "noflux",
                  "neighbor_order": 2, "temperature": float(s["temperature"]),
                  "seed": 17},
        "cells": cells,
        "contact": contact,
    }


def composite_document(slug):
    """A process-bigraph composite document embedding the CPM engine."""
    spec = build_spec(slug)
    return {
        "cpm": {
            "_type": "process",
            "address": CPM_ADDR,
            "config": {"spec": spec, "mcs_per_update": 16},  # 1 paper-MCS/update
            "inputs": {"fates": ["fates"]},
            "outputs": {
                "volumes": ["volumes"],
                "types": ["types"],
                "positions": ["positions"],
                "field_at_cell": ["field_at_cell"],
                "neighbor_secretory": ["neighbor_secretory"],
            },
        },
        "fates": {},
    }


_STUDY_SLUGS = list(STUDIES.keys())


# ONE parameterized composite for all 11 examples, not one-per-study. Every
# example wraps the SAME single CPMProcess; the 11 differ only in the config
# (energies / temperature / area constraint / initial condition) pulled from
# STUDIES[study]. The `study` parameter selects which one — each study passes
# it via baseline[].params, e.g. ``params: {study: annealing}`` — and the
# dashboard renders it as a dropdown of the allowed slugs.
@composite_generator(
    name="gg1993", default_n_steps=16,
    parameters={
        "study": {
            "type": "string",
            "default": _STUDY_SLUGS[0],
            "choices": _STUDY_SLUGS,
            "description": "which Glazier & Graner (1993) paper example to encode",
        },
    },
    description=("Glazier & Graner (1993) CPM example (single CPMProcess). The `study` "
                 "parameter selects which of the 11 paper examples to encode; its exact "
                 "energies, temperature, area constraint and initial condition come from "
                 "viva_cpm_studies.gg1993.params.STUDIES[study]."))
def gg1993(core=None, study=_STUDY_SLUGS[0]):
    if study not in STUDIES:
        raise ValueError(
            f"unknown gg1993 study {study!r}; choices: {', '.join(_STUDY_SLUGS)}")
    return composite_document(study)
