"""Instantiate a finalized cpm_core.World from a load_world sheet spec."""
from __future__ import annotations
from cpm import cpm_core


def world_from_spec(spec: dict):
    po = spec["potts"]
    w = cpm_core.World(tuple(po["dims"]), po["boundary"],
                       int(po["neighbor_order"]), float(po["temperature"]))
    for c in spec["contact"]:
        w.set_contact(int(c["a"]), int(c["b"]), float(c["j"]))
    for cell in spec["cells"]:
        cid = w.add_cell(int(cell["type"]), float(cell["target_volume"]),
                         float(cell["lambda_volume"]),
                         float(cell["target_surface"]), float(cell["lambda_surface"]))
        x0, y0, z0, x1, y1, z1 = cell["seed_block"]
        w.seed_block(cid, x0, y0, z0, x1, y1, z1)
    w.finalize(int(po.get("seed", 0)))
    return w
