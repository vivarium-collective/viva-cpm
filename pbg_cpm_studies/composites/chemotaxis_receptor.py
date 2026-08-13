"""``recruitment_receptor`` — a receptor-level realization of chemotactic
recruitment.

Where ``chemotaxis.recruitment`` gives every responder a hardcoded, always-on
chemotactic bias, this composite makes the response *emergent*: each
responder cell carries a per-cell ``ReceptorSubcell`` (Task 1) that reads the
local cue concentration and computes fractional receptor occupancy via a Hill
relation. Only once a responder's occupancy crosses threshold does it flip
from ``NAIVE_TYPE`` (2) to ``ACTIVATED_TYPE`` (3) — and only type 3 carries a
nonzero chemotactic ``lambda`` in the field's ``chemotaxis`` list. Naive
responders sit adjacent to the gradient exactly like activated ones (contact
J is identical, see ``build_receptor_spec``); the only thing that changes
with activation is whether the cell chemotaxes.

The per-cell wiring (subcell address, ports, and the ``fates`` map pre-init)
comes from ``cpm.coupling.receptor_coupling`` (Task 2); this module supplies
the CPM ``spec`` (built on ``chemotaxis.build_spec``, with the activated
sub-type's contact rows and the field's ``chemotaxis`` list) and assembles
the full composite document.

Referenced from a study's ``baseline[].composite`` as
``pbg_cpm_studies.composites.chemotaxis_receptor.recruitment_receptor``.
"""
from __future__ import annotations

from process_bigraph.composite_generator import composite_generator

from cpm.coupling import receptor_coupling
from pbg_cpm_studies.composites import chemotaxis as CT

NAIVE_TYPE = 2
ACTIVATED_TYPE = 3

# Hill-occupancy defaults for the per-cell ReceptorSubcell. conc_scale brings
# the field's raw amplitude (cue_rate / decay, same units as chemotaxis.py's
# field) into the Hill relation's sensitive window around KD_DEFAULT; hill=1
# and activate_occupancy=0.5 put the activation threshold exactly at
# ligand*conc_scale == kd (see ReceptorSubcell docstring).
KD_DEFAULT = 1.0
HILL_DEFAULT = 2.0
CONC_SCALE_DEFAULT = 0.02
ACTIVATE_OCCUPANCY_DEFAULT = 0.5

CUE_RATE = CT.CUE_RATE
CHEMO_LAMBDA = CT.CHEMO_LAMBDA
SEED = CT.SEED


def _responder_ids(spec):
    """CPM cell ids (1-based, in ``cells`` list order) of the naive-seeded
    responders in ``spec`` -- the cells eligible to activate."""
    return [i + 1 for i, cell in enumerate(spec["cells"]) if cell["type"] == CT.RESPONDER_TYPE]


def build_receptor_spec(*, cue_rate=CUE_RATE, chemo_lambda=CHEMO_LAMBDA, blocked=False, seed=SEED):
    """A ``load_world`` spec for the receptor-recruitment world: naive(2)/
    activated(3) responder sub-types with matched contact-J, and a field
    whose ``chemotaxis`` list only biases the activated type."""
    spec = CT.build_spec(cue_rate=cue_rate, chemo_lambda=chemo_lambda, seed=seed)

    # Add ACTIVATED_TYPE contact rows by copying NAIVE_TYPE's rows pair-by-pair,
    # so J(3, other) == J(2, other) for every other type -- adhesion is not
    # confounded with activation, only chemotaxis differs.
    existing = {(min(r["a"], r["b"]), max(r["a"], r["b"])): r["j"] for r in spec["contact"]}
    contact = list(spec["contact"])
    for other in (0, CT.SOURCE_TYPE, NAIVE_TYPE):
        key = (min(NAIVE_TYPE, other), max(NAIVE_TYPE, other))
        a, b = min(ACTIVATED_TYPE, other), max(ACTIVATED_TYPE, other)
        contact.append({"a": a, "b": b, "j": existing[key]})
    contact.append({"a": ACTIVATED_TYPE, "b": ACTIVATED_TYPE, "j": existing[(NAIVE_TYPE, NAIVE_TYPE)]})
    spec["contact"] = contact

    # Only the activated type chemotaxes; naive (2) is omitted -> lambda 0.
    # ``blocked`` forces the activated lambda to 0 too (response disabled
    # even for cells that do activate).
    eff_lambda = 0.0 if blocked else float(chemo_lambda)
    spec["fields"][0]["chemotaxis"] = (
        [{"type": ACTIVATED_TYPE, "lambda": eff_lambda}] if eff_lambda else []
    )
    return spec


def composite_document(*, cue_rate=CUE_RATE, chemo_lambda=CHEMO_LAMBDA, kd=KD_DEFAULT,
                       seed=SEED, blocked=False):
    spec = build_receptor_spec(cue_rate=cue_rate, chemo_lambda=chemo_lambda, blocked=blocked, seed=seed)
    responder_ids = _responder_ids(spec)
    receptor_config = {
        "kd": float(kd),
        "hill": HILL_DEFAULT,
        "conc_scale": CONC_SCALE_DEFAULT,
        "activate_occupancy": ACTIVATE_OCCUPANCY_DEFAULT,
        "naive_type": NAIVE_TYPE,
        "activated_type": ACTIVATED_TYPE,
    }
    doc = {
        "cpm": {
            "_type": "process",
            "address": CT.CPM_ADDR,
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
    }
    # receptor_<cid> per-responder subcells + the pre-initialized fates store.
    doc.update(receptor_coupling(responder_ids, receptor_config=receptor_config))
    return doc


def meta():
    m = CT.meta()
    m.update({"naive_type": NAIVE_TYPE, "activated_type": ACTIVATED_TYPE,
              "type_names": ["Medium", "Source", "Responder (naive)", "Responder (activated)"]})
    return m


_VIZ = [{"name": "Recruitment over time", "address": "local:ChemotaxisRecruitment"}]


@composite_generator(
    name="recruitment_receptor", default_n_steps=40, visualizations=_VIZ,
    description=("Chemotactic recruitment via a receptor-level response (single CPMProcess "
                 "+ one ReceptorSubcell per responder). Responders activate (naive(2) -> "
                 "activated(3)) once local cue occupancy crosses threshold, and only "
                 "activated cells chemotax. blocked=True zeroes the activated lambda even "
                 "for cells that do activate."),
    parameters={
        "cue_rate": {"type": "float", "default": CUE_RATE,
                     "description": "source secretion rate (0 = no cue)"},
        "chemo_lambda": {"type": "float", "default": CHEMO_LAMBDA,
                         "description": "activated-responder chemotaxis strength"},
        "kd": {"type": "float", "default": KD_DEFAULT,
               "description": "receptor half-occupancy concentration"},
        "blocked": {"type": "boolean", "default": False,
                    "description": "force the activated response off (intervention)"},
    })
def recruitment_receptor(core=None, cue_rate=CUE_RATE, chemo_lambda=CHEMO_LAMBDA,
                         kd=KD_DEFAULT, seed=SEED, blocked=False):
    return composite_document(cue_rate=cue_rate, chemo_lambda=chemo_lambda, kd=kd,
                              seed=seed, blocked=blocked)
