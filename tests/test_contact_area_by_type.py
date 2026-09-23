"""Task 4.1: per-cell contact-area-by-type (Allee surface-fraction substrate)."""
from viva_cpm_studies.influenza import build, sheet, types


def _cells_per_side(spec):
    n = len(spec["cells"])
    return round(n ** 0.5)


def test_contact_area_sums_to_cell_surface_for_interior_cell():
    spec = sheet.build_sheet_spec(0.3)
    w = build.world_from_spec(spec)
    cells_per_side = _cells_per_side(spec)  # 30
    gx = gy = cells_per_side // 2  # interior, far from every sheet edge
    cid = 1 + gy * cells_per_side + gx

    contact = w.cell_contact_area_by_type(cid)
    surfaces = w.cell_surfaces()
    assert sum(contact.values()) == surfaces[cid]
    assert surfaces[cid] > 0


def test_interior_all_h_cell_contacts_mostly_type_h():
    spec = sheet.build_sheet_spec(0.3)
    w = build.world_from_spec(spec)
    cells_per_side = _cells_per_side(spec)
    gx = gy = cells_per_side // 2
    cid = 1 + gy * cells_per_side + gx

    contact = w.cell_contact_area_by_type(cid)
    total = sum(contact.values())
    assert contact.get(types.H, 0) > 0
    # Fully interior confluent cell: contact is (near-)entirely with other H
    # cells; medium (type 0), if present at all, is a small remainder.
    assert contact.get(types.H, 0) >= 0.9 * total


def test_retyping_a_neighbor_shows_up_as_nonzero_new_type_contact():
    spec = sheet.build_sheet_spec(0.3)
    w = build.world_from_spec(spec)
    cells_per_side = _cells_per_side(spec)
    gx = gy = cells_per_side // 2
    cid = 1 + gy * cells_per_side + gx
    neighbor_cid = 1 + (gy - 1) * cells_per_side + gx  # cell directly above

    before = w.cell_contact_area_by_type(cid)
    assert before.get(types.I, 0) == 0

    w.set_cell_type(neighbor_cid, types.I)

    after = w.cell_contact_area_by_type(cid)
    assert after.get(types.I, 0) > 0
    # Retyping doesn't move pixels: total contact area is unchanged.
    assert sum(after.values()) == sum(before.values())
    surfaces = w.cell_surfaces()
    assert sum(after.values()) == surfaces[cid]
