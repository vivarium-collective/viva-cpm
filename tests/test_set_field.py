"""Bulk ``World.set_field`` — overwrite a whole field grid from an external engine.

viva-multicell-blueprint runs ONE authoritative reaction-diffusion field engine
and wants CPM to only READ that field for chemotaxis. ``set_field`` lets the
external engine overwrite a field's grid each step, matching the exact
layout/ordering ``field_conc`` returns, so native diffusion/decay/secretion can
be zeroed and the externally-set field is used verbatim.
"""

import numpy as np

import viva_cpm
from viva_cpm.cpm_core import World


def test_set_field_roundtrips():
    # set_field(grid) then field_conc() returns exactly that grid (same layout).
    w = World((4, 3, 1), "noflux", 2, 10.0)
    fi = w.add_field("c", 0.1, 1e-3)  # native dynamics irrelevant to the roundtrip
    w.finalize(1)

    n = 4 * 3 * 1
    grid = np.arange(n, dtype=np.float32) * 1.25
    w.set_field(fi, grid)

    out = w.field_conc(fi)
    assert len(out) == n
    assert np.allclose(np.asarray(out, dtype=np.float32), grid, atol=1e-6)

    # Also addressable by field name, not just index.
    grid2 = (np.arange(n, dtype=np.float32) + 7.0)
    w.set_field("c", grid2)
    assert np.allclose(np.asarray(w.field_conc(fi), dtype=np.float32), grid2, atol=1e-6)


def test_native_dynamics_zeroed_keeps_field():
    # d=0, decay=0, no secretion: native dynamics must not perturb an
    # externally-set field. set_field each step, then step() (which runs
    # advance_fields internally); the field equals exactly what was set.
    w = World((5, 5, 1), "noflux", 2, 10.0)
    fi = w.add_field("c", 0.0, 0.0)  # d=0, decay=0; no set_secretion -> no source
    w.finalize(1)

    n = 5 * 5
    grid = [float((i * 7) % 11) for i in range(n)]

    for _ in range(5):
        w.set_field(fi, grid)
        w.step(1)

    out = w.field_conc(fi)
    assert np.allclose(np.asarray(out, dtype=np.float32), np.asarray(grid, dtype=np.float32), atol=1e-6), out
