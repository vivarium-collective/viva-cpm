"""``recruitment_adaptive`` — fold-change-detecting recruitment.

Where ``chemotaxis_receptor.recruitment_receptor`` gates activation on raw Hill
occupancy of the local cue (a cell that started in a high background would
already read "active"), this composite swaps in ``AdaptiveReceptorSubcell``
(Task 1): each responder carries a slow per-cell setpoint ``m`` that adapts
toward the local occupancy, so the signal that actually gates activation is
the *change* in occupancy (``response = max(0, theta - m)``), not its
absolute level. A responder sitting in a high but STATIC background
eventually re-adapts (``m -> theta``) and stops activating; only a rising cue
(a source gradient) keeps ``response`` above threshold. ``epsilon=0``
disables adaptation entirely (m frozen at 0), reducing the mechanism back to
pure Hill occupancy -- the knockout condition.

As in ``chemotaxis_receptor``, only ``ACTIVATED_TYPE`` (3) carries a nonzero
chemotactic ``lambda``; naive (2) and activated (3) share identical contact J
so adhesion is not confounded with activation. The per-cell wiring (subcell
address, ports, the pre-inited ``fates`` AND ``adaptation`` map stores) comes
from ``cpm.coupling.adaptation_coupling`` (Task 2); this module supplies the
CPM ``spec`` (built on ``chemotaxis.build_spec``) and assembles the full
composite document plus the five background-cue conditions used to probe
fold-change detection (``CONDITIONS``).

Referenced from a study's ``baseline[].composite`` as
``pbg_cpm_studies.composites.chemotaxis_adaptive.recruitment_adaptive``.
"""
from __future__ import annotations

from process_bigraph.composite_generator import composite_generator

from cpm.coupling import adaptation_coupling
from pbg_cpm_studies.composites import chemotaxis as CT

NAIVE_TYPE = 2
ACTIVATED_TYPE = 3

# Hill-occupancy defaults for the per-cell AdaptiveReceptorSubcell -- mirrors
# chemotaxis_receptor's KD/hill/conc_scale/activate_occupancy defaults; kd here
# matches AdaptiveReceptorSubcell.config_schema's own default (2.9).
KD_DEFAULT = 2.9
HILL_DEFAULT = 2.0
CONC_SCALE_DEFAULT = 0.02
ACTIVATE_OCCUPANCY_DEFAULT = 0.5
EPSILON_DEFAULT = 0.1
BACKGROUND_DEFAULT = 0.0

CUE_RATE = CT.CUE_RATE
CHEMO_LAMBDA = CT.CHEMO_LAMBDA
SEED = CT.SEED


def _responder_ids(spec):
    """CPM cell ids (1-based, in ``cells`` list order) of the naive-seeded
    responders in ``spec`` -- the cells eligible to activate."""
    return [i + 1 for i, cell in enumerate(spec["cells"]) if cell["type"] == CT.RESPONDER_TYPE]


def build_adaptive_spec(*, cue_rate=CUE_RATE, chemo_lambda=CHEMO_LAMBDA, blocked=False, seed=SEED):
    """A ``load_world`` spec for the adaptive-recruitment world: naive(2)/
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
    # even for cells that do activate) -- the receptor_gating intervention.
    spec["fields"][0]["chemotaxis"] = (
        [] if blocked else [{"type": ACTIVATED_TYPE, "lambda": float(chemo_lambda)}]
    )
    return spec


def composite_document(*, cue_rate=CUE_RATE, chemo_lambda=CHEMO_LAMBDA, kd=KD_DEFAULT,
                       epsilon=EPSILON_DEFAULT, background=BACKGROUND_DEFAULT,
                       seed=SEED, blocked=False):
    spec = build_adaptive_spec(cue_rate=cue_rate, chemo_lambda=chemo_lambda, blocked=blocked, seed=seed)
    responder_ids = _responder_ids(spec)
    receptor_config = {
        "kd": float(kd),
        "hill": HILL_DEFAULT,
        "conc_scale": CONC_SCALE_DEFAULT,
        "activate_occupancy": ACTIVATE_OCCUPANCY_DEFAULT,
        "epsilon": float(epsilon),
        "background": float(background),
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
    # receptor_<cid> per-responder subcells + the pre-initialized fates/adaptation stores.
    doc.update(adaptation_coupling(responder_ids, receptor_config=receptor_config))
    return doc


def meta():
    m = CT.meta()
    m.update({"naive_type": NAIVE_TYPE, "activated_type": ACTIVATED_TYPE,
              "type_names": ["Medium", "Source", "Responder (naive)", "Responder (activated)"]})
    return m


_VIZ = [{"name": "Recruitment over time", "address": "local:ChemotaxisRecruitment"}]

# Five background-cue conditions probing fold-change detection: recruitment
# should track a RISING cue regardless of background level (low_bg/mid_bg/
# high_bg all recruit), while receptor_gating (high_bg_blocked, lambda forced
# to 0) and the adaptation knockout (high_bg_knockout, epsilon=0 -> no
# fold-change, pure Hill occupancy re-saturates in high background) do not.
# NOTE: the background values below are PLACEHOLDER starting points for
# Task 4's empirical calibration, not tuned final values.
CONDITIONS = {
    "low_bg":           {"background": 0.0,   "epsilon": 0.1, "blocked": False},
    "mid_bg":           {"background": 80.0,  "epsilon": 0.1, "blocked": False},
    "high_bg":          {"background": 200.0, "epsilon": 0.1, "blocked": False},
    "high_bg_blocked":  {"background": 200.0, "epsilon": 0.1, "blocked": True},   # receptor_gating
    "high_bg_knockout": {"background": 200.0, "epsilon": 0.0, "blocked": False},  # adaptation off
}


@composite_generator(
    name="recruitment_adaptive", default_n_steps=40, visualizations=_VIZ,
    description=("Adaptive receptor recruitment (fold-change detection): single CPMProcess "
                 "+ one AdaptiveReceptorSubcell per responder. Responders activate (naive(2) "
                 "-> activated(3)) once response = max(0, theta - m) crosses threshold, where "
                 "m is a slow per-cell setpoint adapting toward occupancy theta -- so a rising "
                 "cue recruits regardless of background level, while a static background does "
                 "not. blocked=True zeroes the activated lambda (receptor_gating); epsilon=0 "
                 "disables adaptation (knockout, pure Hill occupancy)."),
    parameters={
        "cue_rate": {"type": "float", "default": CUE_RATE,
                     "description": "source secretion rate (0 = no cue)"},
        "chemo_lambda": {"type": "float", "default": CHEMO_LAMBDA,
                         "description": "activated-responder chemotaxis strength"},
        "kd": {"type": "float", "default": KD_DEFAULT,
               "description": "receptor half-occupancy concentration"},
        "epsilon": {"type": "float", "default": EPSILON_DEFAULT,
                    "description": "adaptation rate (0 = knockout, pure Hill occupancy)"},
        "background": {"type": "float", "default": BACKGROUND_DEFAULT,
                       "description": "static background cue level added to the local field"},
        "blocked": {"type": "boolean", "default": False,
                    "description": "force the activated response off (intervention)"},
    })
def recruitment_adaptive(core=None, cue_rate=CUE_RATE, chemo_lambda=CHEMO_LAMBDA,
                         kd=KD_DEFAULT, epsilon=EPSILON_DEFAULT, background=BACKGROUND_DEFAULT,
                         seed=SEED, blocked=False):
    return composite_document(cue_rate=cue_rate, chemo_lambda=chemo_lambda, kd=kd,
                              epsilon=epsilon, background=background, seed=seed, blocked=blocked)
