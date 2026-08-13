"""Regression + extension tests for pbg_cpm_studies.chemotaxis.metrics."""
from pbg_cpm_studies.chemotaxis import metrics as M


def test_single_responder_type_unchanged():
    """The existing chemotaxis studies call with a single responder_type
    (source at origin, radius covers only the near responder). This must
    stay byte-for-byte identical after generalizing to multi-type."""
    coms = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (5.0, 0.0, 0.0), (50.0, 0.0, 0.0)]
    types = [0, 1, 2, 2]
    idx = M.recruitment_index_from_coms(coms, types, source_type=1,
                                        responder_type=2, radius=10.0)
    assert idx == 0.5  # 1 of 2 type-2 responders within radius


def test_responder_types_counts_multiple_types():
    """A composite with naive(2)+activated(3) responders: responder_types={2,3}
    counts both, unlike the single-type default."""
    coms = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (5.0, 0.0, 0.0),
            (5.0, 0.0, 0.0), (50.0, 0.0, 0.0)]
    types = [0, 1, 2, 3, 3]
    single = M.recruitment_index_from_coms(coms, types, source_type=1,
                                           responder_type=2, radius=10.0)
    both = M.recruitment_index_from_coms(coms, types, source_type=1,
                                         responder_types={2, 3}, radius=10.0)
    assert single == 1.0  # only the one type-2 cell, and it's within radius
    assert both == 2 / 3  # 2 of 3 type-{2,3} cells within radius


def test_recruitment_index_all_responders_wrapper():
    coms = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (5.0, 0.0, 0.0),
            (5.0, 0.0, 0.0), (50.0, 0.0, 0.0)]
    types = [0, 1, 2, 3, 3]
    assert M.recruitment_index_all_responders(
        coms, types, source_type=1, radius=10.0) == 2 / 3
