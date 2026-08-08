"""process-bigraph composite factories for the chemotactic-recruitment slice.

Three registry-resolvable conditions form the experimental logic — recruitment
requires BOTH a secreted cue AND a competent chemotactic response:

    recruitment_baseline     cue ON  (rate 10), response ON  (lambda 14)  -> recruited
    recruitment_inhibited    cue ON  (rate 10), response OFF (lambda  0)  -> not recruited
                             (the intervention: the response is blocked)
    recruitment_adversarial  cue OFF (rate  0), response ON  (lambda 14)  -> not recruited
                             (negative control: competent responders, no cue)

Each factory returns a composite *document* embedding the real Rust CPM engine
(``cpm.processes.cpm_process.CPMProcess``). Source cells (type 1) form a slab on
the left; responder cells (type 2) start in the gradient zone on the right. The
parameters below are the de-risked recipe (responders sweep from x~43 to x~4
under the baseline, and stay at x~43 when cue or response is removed).

Referenced from a study's ``baseline[].composite`` as
``pbg_cpm_studies.composites.chemotaxis.<name>``.
"""
from __future__ import annotations

import itertools

from process_bigraph.composite_generator import composite_generator

CPM_ADDR = "local:!cpm.processes.cpm_process.CPMProcess"

SOURCE_TYPE = 1
RESPONDER_TYPE = 2
NX, NY = 80, 40

# de-risked recipe (see pbg_cpm_studies/chemotaxis validation)
CUE_RATE = 10.0
CHEMO_LAMBDA = 14.0
FIELD_D = 0.22
FIELD_DECAY = 0.001
TEMPERATURE = 10.0
SEED = 17
RECRUITMENT_RADIUS = 15.0


def _cell(t, x0, y0, x1, y1, *, target_volume=30, lambda_volume=1.2):
    return {
        "type": t,
        "target_volume": float(target_volume),
        "lambda_volume": float(lambda_volume),
        "target_surface": 0.0,
        "lambda_surface": 0.0,
        "seed_block": [x0, y0, 0, x1, y1, 1],
    }


def build_spec(*, cue_rate=CUE_RATE, chemo_lambda=CHEMO_LAMBDA, seed=SEED):
    """A load_world spec for the recruitment world at one (cue, response) setting."""
    cells = [_cell(SOURCE_TYPE, 5, 10, 15, 30, target_volume=200, lambda_volume=1.0)]
    # six responder cells spread through the gradient zone (right of source)
    for (x0, y0) in itertools.product((30, 40, 50), (8, 26)):
        cells.append(_cell(RESPONDER_TYPE, x0, y0, x0 + 7, y0 + 7))
    return {
        "potts": {"dims": [NX, NY, 1], "boundary": "noflux",
                  "neighbor_order": 2, "temperature": TEMPERATURE, "seed": seed},
        "cells": cells,
        "contact": [
            {"a": 0, "b": SOURCE_TYPE, "j": 12.0},
            {"a": 0, "b": RESPONDER_TYPE, "j": 12.0},
            {"a": SOURCE_TYPE, "b": SOURCE_TYPE, "j": 4.0},
            {"a": RESPONDER_TYPE, "b": RESPONDER_TYPE, "j": 6.0},
            {"a": SOURCE_TYPE, "b": RESPONDER_TYPE, "j": 9.0},
        ],
        "fields": [{
            "name": "C", "d": FIELD_D, "decay": FIELD_DECAY,
            "secretion": [{"type": SOURCE_TYPE, "rate": float(cue_rate)}],
            "chemotaxis": [{"type": RESPONDER_TYPE, "lambda": float(chemo_lambda)}],
        }],
    }


def composite_document(*, cue_rate=CUE_RATE, chemo_lambda=CHEMO_LAMBDA, seed=SEED):
    spec = build_spec(cue_rate=cue_rate, chemo_lambda=chemo_lambda, seed=seed)
    return {
        "cpm": {
            "_type": "process",
            "address": CPM_ADDR,
            "config": {"spec": spec, "mcs_per_update": 16, "n_fields": 1},
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


def meta():
    """Static description a viewer / metric / claim binding can key on."""
    return {
        "dims": [NX, NY, 1],
        "source_type": SOURCE_TYPE,
        "responder_type": RESPONDER_TYPE,
        "recruitment_radius": RECRUITMENT_RADIUS,
        "type_names": ["Medium", "Source", "Responder"],
    }


_VIZ = [{"name": "Recruitment over time", "address": "local:ChemotaxisRecruitment"}]


# ONE parameterized composite for the whole slice, not one-per-condition. It wraps
# a single CPMProcess; the conditions (baseline / inhibited / adversarial) are just
# different (cue_rate, chemo_lambda) PARAMS the studies pass — not separate
# composites. This is the idiomatic pattern for a family of configs of one process.
@composite_generator(
    name="recruitment", default_n_steps=32, visualizations=_VIZ,
    description=("Chemotactic recruitment (single CPMProcess). Source cells secrete a "
                 "cue; responders chemotax up it. cue_rate=0 removes the cue; "
                 "chemo_lambda=0 blocks the response."),
    parameters={
        "cue_rate": {"type": "float", "default": CUE_RATE,
                     "description": "source secretion rate (0 = no cue)"},
        "chemo_lambda": {"type": "float", "default": CHEMO_LAMBDA,
                         "description": "responder chemotaxis strength (0 = response blocked)"},
    })
def recruitment(core=None, cue_rate=CUE_RATE, chemo_lambda=CHEMO_LAMBDA):
    return composite_document(cue_rate=cue_rate, chemo_lambda=chemo_lambda)
