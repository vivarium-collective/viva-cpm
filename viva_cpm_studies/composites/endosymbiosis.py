"""Composite generator for eukaryogenesis S3 — protrusions to enclosure.

``protrusion_enclosure`` embeds the real CPM engine
(``ProtrusionEnclosureProcess``) as campaign composite B: a deformable host and a
compact symbiont whose host-symbiont adhesion ``j_hs`` is the sweepable gate
parameter. Running it emits the wrapping observables at ``enc/*``; sweeping
``j_hs`` across ``j_sM - j_hM`` drives the enclosure transition that the S3 gate
(G3 Enclosure) grades.
"""
from __future__ import annotations

from process_bigraph.composite_generator import composite_generator

ENCLOSURE_ADDR = "local:!viva_cpm_studies.endosymbiosis.process.ProtrusionEnclosureProcess"


@composite_generator(
    name="protrusion_enclosure", default_n_steps=40,
    parameters={
        "j_hs": {
            "type": "number", "default": 6.0,
            "description": "host-symbiont adhesion energy (low = strong adhesion = enclosure; "
                           "the wrapping transition sits at j_sM - j_hM = 12)",
        },
        "j_sM": {"type": "number", "default": 16.0,
                 "description": "symbiont-medium surface tension"},
        "j_hM": {"type": "number", "default": 4.0,
                 "description": "host-medium surface tension (host cortical tension)"},
        "seed": {"type": "integer", "default": 7, "description": "RNG seed for the realization"},
    },
    description=("Eukaryogenesis S3 (protrusions to enclosure), campaign composite B. A host and a "
                 "symbiont on the Rust CPM engine; sweeping host-symbiont adhesion j_hs across "
                 "j_sM - j_hM moves the system through the wrapping transition. Emits "
                 "enc/enclosure_fraction (1 = fully enclosed) and enc/enclosed."))
def protrusion_enclosure(core=None, j_hs=6.0, j_sM=16.0, j_hM=4.0, seed=7):
    return {
        "enclosure": {
            "_type": "process",
            "address": ENCLOSURE_ADDR,
            "config": {"j_hs": float(j_hs), "j_sM": float(j_sM), "j_hM": float(j_hM),
                       "seed": int(seed), "mcs_per_update": 60},
            "inputs": {},
            "outputs": {
                "enclosure_fraction": ["enc", "enclosure_fraction"],
                "enclosed": ["enc", "enclosed"],
                "host_contact_faces": ["enc", "host_contact_faces"],
                "medium_contact_faces": ["enc", "medium_contact_faces"],
                "sweeps": ["enc", "sweeps"],
            },
        },
        "enc": {},
    }
