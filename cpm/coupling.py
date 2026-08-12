"""Down-scale coupling helper: wire per-cell subcellular processes to the CPM
`fates` port. Encapsulates the map-key contract (a `fates` write to an absent
key is dropped on apply, so the store is pre-initialized per cell). When the
bigraph-schema keyed-map affordance (H5) lands, the pre-init here is removed and
the `fates` port declared with the affordance instead.
"""
from __future__ import annotations

RECEPTOR_ADDR = "local:!cpm.subcellular.receptor.ReceptorSubcell"


def receptor_coupling(cell_ids, *, receptor_config, receptor_addr=RECEPTOR_ADDR):
    frag = {}
    naive = int(receptor_config["naive_type"])
    fates = {}
    for cid in cell_ids:
        key = str(int(cid))
        frag[f"receptor_{cid}"] = {
            "_type": "process",
            "address": receptor_addr,
            "config": dict(receptor_config),
            "inputs": {"ligand": ["field_at_cell", key]},
            "outputs": {"fate": ["fates", key]},
        }
        fates[key] = naive  # H4 pre-init fallback
    frag["fates"] = fates
    return frag
