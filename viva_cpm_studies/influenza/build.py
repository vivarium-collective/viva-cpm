"""Instantiate a finalized cpm_core.World from a load_world sheet spec."""
from __future__ import annotations
from viva_cpm import cpm_core


def world_from_spec(spec: dict, finalize: bool = True):
    """Build a `cpm_core.World` from a sheet spec.

    `finalize=True` (default, matches every existing caller): the returned
    world is finalized and ready to `step()`. Pass `finalize=False` when a
    caller needs to add fields first (`World.add_field` is only valid
    pre-finalize, see `fields.add_virus_field`); the caller is then
    responsible for calling `world.finalize(seed)` itself before stepping.
    """
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
    if finalize:
        w.finalize(int(po.get("seed", 0)))
    return w
