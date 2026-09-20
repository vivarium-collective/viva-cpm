import cpm
from cpm.cpm_core import World


def test_field_value_and_deposit_at_point():
    w = World((10, 10, 1), "noflux", 2, 10.0)
    idx = w.add_field("c", 0.0, 0.0)  # no diffusion, no decay
    w.finalize(1)
    assert w.field_value_at(idx, 3, 4, 0) == 0.0
    w.field_add_source_at(idx, 3, 4, 0, 2.5)
    assert abs(w.field_value_at(idx, 3, 4, 0) - 2.5) < 1e-9
    # out-of-range coordinates are clamped, not a panic
    w.field_add_source_at(idx, 99, 99, 0, 1.0)
    assert w.field_value_at(idx, 9, 9, 0) >= 1.0
