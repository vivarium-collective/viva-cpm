from cpm import cpm_core


def test_default_scale_matches_unscaled_secretion():
    w = cpm_core.World((10, 10, 1), "noflux", 2, 10.0)
    a = w.add_cell(1, 9.0, 1.0, 0.0, 0.0)
    w.set_contact(0, 1, 8.0)
    w.seed_block(a, 3, 3, 0, 6, 6, 1)

    fi = w.add_field("virus", 0.0, 0.0)
    w.set_secretion(fi, 1, 10.0)

    w.finalize(11)
    w.advance_fields(5)

    conc = w.field_conc(fi)
    assert max(conc) > 0.0, "unscaled secretion should raise concentration"
    assert max(conc) == 50.0  # 10.0 * dt(1.0) * 5 steps, no diffusion/decay


def test_scale_half_halves_secretion():
    w = cpm_core.World((10, 10, 1), "noflux", 2, 10.0)
    a = w.add_cell(1, 9.0, 1.0, 0.0, 0.0)
    w.set_contact(0, 1, 8.0)
    w.seed_block(a, 3, 3, 0, 6, 6, 1)

    fi = w.add_field("virus", 0.0, 0.0)
    w.set_secretion(fi, 1, 10.0)
    w.set_cell_secretion_scale(fi, a, 0.5)

    w.finalize(11)
    w.advance_fields(5)

    conc = w.field_conc(fi)
    assert abs(max(conc) - 25.0) < 1e-4  # half of the unscaled 50.0


def test_scale_zero_suppresses_secretion():
    w = cpm_core.World((10, 10, 1), "noflux", 2, 10.0)
    a = w.add_cell(1, 9.0, 1.0, 0.0, 0.0)
    w.set_contact(0, 1, 8.0)
    w.seed_block(a, 3, 3, 0, 6, 6, 1)

    fi = w.add_field("virus", 0.0, 0.0)
    w.set_secretion(fi, 1, 10.0)
    w.set_cell_secretion_scale(fi, a, 0.0)

    w.finalize(11)
    w.advance_fields(5)

    conc = w.field_conc(fi)
    assert max(conc) == 0.0


def test_unrelated_cell_unaffected_by_another_cells_scale():
    w = cpm_core.World((16, 10, 1), "noflux", 2, 10.0)
    a = w.add_cell(1, 9.0, 1.0, 0.0, 0.0)
    b = w.add_cell(1, 9.0, 1.0, 0.0, 0.0)
    w.set_contact(0, 1, 8.0)
    w.seed_block(a, 1, 1, 0, 4, 4, 1)
    w.seed_block(b, 11, 1, 0, 14, 4, 1)

    fi = w.add_field("virus", 0.0, 0.0)
    w.set_secretion(fi, 1, 10.0)
    w.set_cell_secretion_scale(fi, a, 0.0)

    w.finalize(11)
    w.advance_fields(1)

    mean_a = w.field_mean_at_cell(fi, a)
    mean_b = w.field_mean_at_cell(fi, b)
    assert mean_a == 0.0, "a's secretion should be suppressed"
    assert abs(mean_b - 10.0) < 1e-4, "b's secretion should be untouched"
