import numpy as np
import process_bigraph as pb
from pbg_cpm_studies.influenza.epithelium_process import EpitheliumProcess
from pbg_cpm_studies.influenza import build, immune


def _spec():
    return immune.build_cytotoxic_scenario_spec(
        epithelial_cells_per_side=4, n_infected=1, n_macrophages=2,
        n_nk=2, n_cd8=2, margin_sites=10, separation_sites=8, seed=1)


def test_epithelium_process_exposes_all_fields_and_deposits():
    # NOTE (controller ruling R2, task-1.2-brief.md): unlike the brief's
    # original `spec["fields"] = None` line, the spec here carries NO
    # "fields" key at all -- EpitheliumProcess.initialize adds the four
    # fields itself, imperatively, before finalize.
    spec = _spec()
    core = pb.allocate_core()
    p = EpitheliumProcess({"spec": spec, "mcs_per_update": 1, "n_fields": 4}, core=core)
    out = p.update({"fates": {}, "field_deposit": []}, 1.0)

    # per-cell readout is a list of 4 field means per cell
    any_cell = next(iter(out["field_at_cell_all"]))
    assert len(out["field_at_cell_all"][any_cell]) == 4

    # raw field + dims outputs
    assert len(out["dims"]) == 3
    assert len(out["virus_field"]) == out["dims"][0] * out["dims"][1] * out["dims"][2]

    # depositing chemokine (field idx 2) at a point then reading it back rises
    p.update({"fates": {}, "field_deposit": [[2, 20, 20, 0, 5.0]]}, 1.0)
    assert p.world.field_value_at(2, 20, 20, 0) > 0.0


def test_epithelium_process_fates_deplete_living_epithelium():
    spec = _spec()
    core = pb.allocate_core()
    p = EpitheliumProcess({"spec": spec, "mcs_per_update": 1, "n_fields": 4,
                           "enable": ["infection", "ifn", "death", "ros", "allee"]},
                          core=core)
    # feed a large ROS X so ROS death fires
    living0 = None
    D = 0
    H = I = 0
    for _ in range(25):
        out = p.update({"fates": {}, "field_deposit": [],
                        "ode_state": {"X": 50.0}, "immune_kills": []}, 1.0)
        H, I, D = out["counts"]["H"], out["counts"]["I"], out["counts"]["D"]
        living0 = living0 if living0 is not None else H + I
    assert D > 0 and (H + I) < living0


def test_epithelium_process_emits_ode_inputs():
    spec = _spec()
    core = pb.allocate_core()
    p = EpitheliumProcess({"spec": spec, "mcs_per_update": 1, "n_fields": 4}, core=core)
    out = p.update({"fates": {}, "field_deposit": [[0, 20, 20, 0, 50.0]],
                    "ode_state": {"X": 0.0}, "immune_kills": []}, 1.0)
    ode_inputs = out["ode_inputs"]
    assert set(ode_inputs.keys()) == {
        "H", "I", "DH", "V", "F", "C", "L", "B_ei", "G_ki"}
    assert ode_inputs["V"] > 0.0


def _il10_scenario_spec():
    # No macrophages/NK/CD8: isolates IL-10 secretion to H cells only, so the
    # "chemokine" gate's effect on the uninfected-H IL-10 loop (run.py:
    # 1999-2010) isn't swamped by the ungated macrophage IL-10 contribution
    # (macrophage secretion-scale regulation is out of scope for this
    # process either way -- see epithelium_process.py module docstring).
    return immune.build_cytotoxic_scenario_spec(
        epithelial_cells_per_side=4, n_infected=0, n_macrophages=0,
        n_nk=0, n_cd8=0, margin_sites=10, separation_sites=8, seed=1)


def _drive_il10(enable):
    core = pb.allocate_core()
    spec = _il10_scenario_spec()
    p = EpitheliumProcess({"spec": spec, "mcs_per_update": 1, "n_fields": 4,
                           "enable": enable}, core=core)
    dims = p.world.dims()
    # Strong local IFN under the whole patch so resist -> ~1 for every H
    # cell (cell_resistance = ifn/(a_rf+ifn) -> 1 as ifn grows), so the
    # (1-resist) IL-10 gate (applied only when "chemokine" is enabled)
    # visibly suppresses secretion relative to the unregulated default.
    deposits = [[1, x, y, 0, 500.0] for x in range(0, dims[0], 2)
                for y in range(0, dims[1], 2)]
    out = None
    for _ in range(10):
        out = p.update({"fates": {}, "field_deposit": deposits,
                        "ode_state": {"X": 0.0}, "immune_kills": []}, 1.0)
    return out["ode_inputs"]["L"]


def test_epithelium_process_chemokine_gates_uninfected_il10_secretion():
    # run.py:1999-2010: the uninfected-H IL-10 secretion-scale gate
    # (signaling.uninfected_il10_scale(1-resist)) is nested INSIDE
    # `if "chemokine" in enable:`, the same block as the macrophage
    # secretion loop -- verified by indentation. `world.set_secretion`
    # wires the IL-10 field's H-cell base rate UNCONDITIONALLY at
    # construction (fields.add_il10_field), and the engine's per-cell
    # secretion scale defaults to 1.0 (UNREGULATED, i.e. NOT zero) until
    # explicitly set -- so with "chemokine" off the source does NOT
    # secrete zero IL-10, it secretes at the un-gated (scale=1.0) rate.
    # The observable, source-faithful effect of the gate is therefore:
    # "chemokine" ON (high local IFN -> resist~1 -> scale~0) accumulates
    # LESS IL-10 than "chemokine" OFF (scale stays at the ungated
    # default) over the same MCS/IFN schedule.
    base_enable = ["infection", "ifn", "death", "allee"]
    L_off = _drive_il10(base_enable)
    L_on = _drive_il10(base_enable + ["chemokine"])
    assert L_off > 0.0
    assert L_on >= 0.0
    assert L_on < L_off
