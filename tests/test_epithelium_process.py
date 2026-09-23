import numpy as np
import process_bigraph as pb
from viva_cpm_studies.influenza.epithelium_process import EpitheliumProcess
from viva_cpm_studies.influenza import build, immune, types


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


# --- Task 1.4: Phase-1 parity gate --------------------------------------


def _matching_scene(*, cells_per_side, seed, init_infection_frac,
                    n_macrophages=4, n_nk=4, n_cd8=4,
                    margin_sites=10, separation_sites=8):
    """Build the EXACT scene `run.run_full_model` builds for the
    `init_infection_frac` scenario branch (run.py:1726-1787): the same
    `immune.build_cytotoxic_scenario_spec` call (matching run_full_model's
    own defaults for the immune-cluster kwargs it doesn't vary in this
    parity test), then the SAME init_infection_frac re-scatter -- reset
    every epithelial cell to H, then re-infect `n_infected` of them chosen
    by `np.random.default_rng(seed + 7).choice(...)` -- so an
    EpitheliumProcess built from this spec sees the identical initial
    lesion (and, via `spec["potts"]["seed"] == seed`, the identical Potts
    RNG stream) as run_full_model's own world.
    """
    tot_cell = cells_per_side * cells_per_side
    n_infected = int(round(float(init_infection_frac) * tot_cell))
    n_infected = max(0, min(n_infected, tot_cell))

    spec = immune.build_cytotoxic_scenario_spec(
        epithelial_cells_per_side=cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, n_nk=n_nk, n_cd8=n_cd8,
        margin_sites=margin_sites, separation_sites=separation_sites, seed=seed)

    if init_infection_frac is not None and n_infected > 0:
        epi_idx = [i for i, c in enumerate(spec["cells"])
                   if c["type"] in (types.H, types.I)]
        for i in epi_idx:
            spec["cells"][i]["type"] = types.H
        if epi_idx:
            scatter_rng = np.random.default_rng(seed + 7)
            chosen = scatter_rng.choice(np.array(epi_idx),
                                        size=min(n_infected, len(epi_idx)),
                                        replace=False)
            for i in chosen:
                spec["cells"][int(i)]["type"] = types.I
    return spec


def _drive_epithelium_process(*, cells_per_side, steps, seed, init_infection_frac,
                              enable, mcs_per_step=7):
    """Build the same scene `run_full_model` builds (`_matching_scene`),
    construct an `EpitheliumProcess` with matching config, and step it
    `mcs_per_step` MCS per record for `steps - 1` records -- mirroring
    `run_full_model`'s own outer loop (`_record(0)` is the seeded initial
    state with no MCS advanced, then `for i in range(1, steps): step
    mcs_per_step MCS; _record(i)`, run.py:2148-2154), so the LAST record
    here (index `steps - 1`) lines up MCS-for-MCS with `r["counts"][...
    ][-1]`. `immune_kills=[]` and `ode_state={"X": 0.0}` keep the
    immune-kill merge and ROS death pipeline stages inert (both are also
    absent from `enable` in the Task 1.4 parity test, so this is belt and
    suspenders) for a clean epithelial-only comparison.

    Returns `(dead_final, H_series, I_series, D_series)` -- the H/I/D
    counts recorded after each of the `steps - 1` update() calls.
    """
    core = pb.allocate_core()
    spec = _matching_scene(cells_per_side=cells_per_side, seed=seed,
                           init_infection_frac=init_infection_frac)
    p = EpitheliumProcess({
        "spec": spec, "mcs_per_update": 1, "mcs_per_step": mcs_per_step,
        "n_fields": 4, "enable": list(enable),
    }, core=core)

    H_series, I_series, D_series = [], [], []
    for _ in range(1, steps):
        out = p.update({"fates": {}, "field_deposit": [],
                        "ode_state": {"X": 0.0}, "immune_kills": []}, 1.0)
        H_series.append(out["counts"]["H"])
        I_series.append(out["counts"]["I"])
        D_series.append(out["counts"]["D"])
    dead_final = D_series[-1] if D_series else 0
    return dead_final, H_series, I_series, D_series


def test_epithelium_process_parity_with_run_full_model_epithelial_only():
    # ACHIEVED PARITY (2026-09-20, this seed + 3 additional seed/scale
    # combos probed during development -- see task-1.4-report.md):
    # dead_proc == dead_ref EXACTLY (|delta| = 0, tolerance = 2) -- the
    # per-record dead-count series matches run_full_model's `counts.dead`
    # series element-for-element, not just the final value. With
    # "killing"/"ros"/"ode"/"recruitment"/"macrophage"/"nk_cd8" all off and
    # immune_kills=[]/X=0.0, the two drivers consume IDENTICAL RNG streams
    # (infect_rng/death_rng/allee_rng seeded seed+1/+2/+3, world Potts RNG
    # seeded from the shared spec["potts"]["seed"] == seed) over the
    # IDENTICAL scene and MCS cadence, so exact agreement is expected, not
    # just "within tolerance".
    from viva_cpm_studies.influenza import run
    enable = ["infection", "ifn", "death", "allee"]   # ros/ode/immune off for a clean parity
    r = run.run_full_model(cells_per_side=8, steps=20, seed=3,
                           init_infection_frac=0.1, enable=enable)
    dead_ref = r["counts"]["dead"][-1]
    dead_proc, H, I, D = _drive_epithelium_process(
        cells_per_side=8, steps=20, seed=3, init_infection_frac=0.1, enable=enable)

    tol = max(2, int(0.1 * max(dead_ref, 1)))
    assert abs(dead_proc - dead_ref) <= tol
    # the achieved delta, logged for visibility when this test is run -v
    assert abs(dead_proc - dead_ref) == 0

    # Element-for-element series parity (the whole point of the gate --
    # not just the final scalar). run_full_model's counts series index 0 is
    # the seeded initial state (before any MCS); indices 1..steps-1 are
    # after each record's mcs_per_step MCS -- exactly the records
    # `_drive_epithelium_process`'s H/I/D series holds one-to-one, hence the
    # `[1:]` slice.
    assert H == r["counts"]["uninfected"][1:]
    assert I == r["counts"]["infected"][1:]
    assert D == r["counts"]["dead"][1:]
