from pbg_cpm_studies.composites.influenza import full_model_composite_document
from pbg_cpm_studies.core import build_core
from process_bigraph import Composite


# Task 3.4 (influenza-immune-process-composite epic, final task): the end-to-
# end gate proving `run_full_model_composite` (the pb-Composite driver, Task
# 3.3) is faithful on the EPITHELIAL fate cascade AND that it removes the
# fixed-pool immune-count cap `run_full_model` carries (`recruit_pool_per_
# type`, Task 8.4's reserve-pool device -- see run.py's `RECRUIT_RESERVE_TYPE`
# module note). The composite's immune layer is OFF-LATTICE agents with
# `local_ratio=1.0` and no nearby-surrogate (see immune_process.py's module
# docstring "Recruitment" section) -- a deliberate, documented divergence from
# `run_full_model`'s NK/CD8 `local_ratio` split and nearby-kill term, so
# immune-count/kill-rate PARITY with `run_full_model` is not expected or
# tested here; only the epithelial fate cascade (H/I/D) is.


def test_full_model_composite_builds_and_steps():
    doc = full_model_composite_document(cells_per_side=6, seed=1, init_infection_frac=0.1)
    assert set(doc).issuperset({"epithelium", "immune", "ode"})
    core = build_core()
    comp = Composite({"state": doc}, core=core)
    comp.run(1)   # one interval, no exception


def test_run_full_model_composite_shape_and_decline():
    from pbg_cpm_studies.influenza import run
    r = run.run_full_model_composite(cells_per_side=8, steps=12, seed=1, init_infection_frac=0.1)
    for k in ("counts", "fields", "ode", "t_days", "params"):
        assert k in r
    u = r["counts"]["uninfected"]
    assert u[-1] <= u[0]


def test_run_full_model_composite_immune_killing_wired():
    """Task 3.3 hard requirement: `epithelial_positions`/`infected_ids` must be
    computed + wired every MCS so ImmuneProcess's proximity-killing loop is
    not a silent no-op. Drives the SAME three processes with the SAME cadence
    `run_full_model_composite` uses (epithelium+immune every MCS, ODE once
    per `mcs_per_step`-MCS record -- the ODE coupling matters here: NK/CD8
    chemotax up the CHEMOKINE field, which macrophages only secrete once the
    ODE's dynamic `sig_1` rises above 0, so skipping the ODE leaves `sig_1`
    pinned at its 0.0 IC and the immune cluster never closes in), but reads
    `immune_kills` directly each MCS (`run_full_model_composite`'s own result
    dict doesn't expose the per-MCS kill list, only aggregate counts/fields/
    ode series) -- proof the wiring is live, not just that the driver runs."""
    from pbg_cpm_studies.composites.influenza import full_model_composite_document
    from pbg_cpm_studies.core import build_core
    from pbg_cpm_studies.influenza import types as inf_types

    mcs_per_step = 7
    # separation_sites/margin_sites shrunk from the defaults (8/10) so the
    # NK/CD8 cluster starts closer to the epithelial patch -- purely a test
    # speed-up (fewer MCS needed for chemotaxis to close the gap), not a
    # change to any kill-rate constant.
    doc = full_model_composite_document(cells_per_side=6, seed=3, init_infection_frac=0.3,
                                        separation_sites=2, margin_sites=6,
                                        mcs_per_step=mcs_per_step)
    doc["epithelium"]["config"]["mcs_per_step"] = 1
    core = build_core()
    comp = Composite({"state": doc}, core=core)
    epi_proc = comp.state["epithelium"]["instance"]
    imm_proc = comp.state["immune"]["instance"]
    ode_proc = comp.state["ode"]["instance"]

    store = {"fates": {}, "field_deposit": [], "ode_state": {}, "immune_kills": [],
             "immune_agents": list(doc["immune_agents"]),
             "sig_1": 0.0, "recruit_drivers": {}, "margin_box": doc["margin_box"],
             "ode_inputs": {}, "immune_counts": {}}

    any_kill = False
    for _record in range(1, 9):  # 8 records * 7 MCS/record = 56 raw MCS
        for _mcs in range(mcs_per_step):
            epi_out = epi_proc.update({
                "fates": store["fates"], "field_deposit": store["field_deposit"],
                "ode_state": store["ode_state"], "immune_kills": store["immune_kills"],
            }, 1)
            types_now = epi_out["types"]
            positions_now = epi_out["positions"]
            infected_ids = [cid for cid in range(1, len(types_now)) if types_now[cid] == inf_types.I]
            epithelial_positions = {str(cid): [float(positions_now[cid][0]), float(positions_now[cid][1])]
                                    for cid in infected_ids}
            store["ode_inputs"] = epi_out["ode_inputs"]

            imm_out = imm_proc.update({
                "chemo_field": epi_out["chemo_field"], "virus_field": epi_out["virus_field"],
                "dims": epi_out["dims"], "immune_agents": store["immune_agents"],
                "epithelial_positions": epithelial_positions, "infected_ids": infected_ids,
                "sig_1": store["sig_1"], "recruit_drivers": store["recruit_drivers"],
                "margin_box": store["margin_box"], "apply_recruitment": True,
            }, 1)
            store["immune_agents"] = imm_out["immune_agents"]
            store["immune_kills"] = imm_out["immune_kills"]
            store["field_deposit"] = imm_out["field_deposit"]
            store["immune_counts"] = imm_out["immune_counts"]
            if imm_out["immune_kills"]:
                any_kill = True

        oi = store["ode_inputs"]
        ode_out = ode_proc.update({
            "H": oi.get("H", 0.0), "I": oi.get("I", 0.0), "DH": oi.get("DH", 0.0),
            "V": oi.get("V", 0.0), "F": oi.get("F", 0.0), "C": oi.get("C", 0.0), "L": oi.get("L", 0.0),
            "B_ei": oi.get("B_ei", 0.0), "G_ki": oi.get("G_ki", 0.0),
            "M": store["immune_counts"].get("M", 0), "K": store["immune_counts"].get("K", 0),
            "E": store["immune_counts"].get("E", 0),
        }, mcs_per_step)
        store["ode_state"] = ode_out["ode_state"]
        store["recruit_drivers"] = ode_out["recruit_drivers"]
        store["sig_1"] = ode_out["sig_1"]

    assert any_kill, "expected immune killing to fire at least once over the run"
    # And the ODE-input plumbing this hinges on is non-trivial (sig_1 rose off
    # its 0.0 IC, immune agents present) -- not just an accidental early kill.
    assert store["sig_1"] > 0.0
    assert sum(store["immune_counts"].values()) > 0


def test_composite_epithelial_parity_reduced_scale():
    """Task 3.4 -- the end-to-end parity gate: at reduced scale, the
    composite's EPITHELIAL fate cascade (infection -> Allee/ROS/apoptosis/
    immune-kill death) tracks `run_full_model`'s hand-loop reference within a
    documented tolerance. Only `dead[-1]` (the terminal DEAD count) is
    compared -- NOT any immune count -- because the two drivers' immune
    layers are deliberately different (see this module's header comment):
    `run_full_model`'s NK/CD8 split 75% local CPM cells / 25% into a
    well-mixed "nearby" surrogate consumed by `killing.nearby_kill_rate`,
    while the composite's off-lattice `ImmuneProcess` sends the FULL
    ODE-predicted inflow to local agents (`local_ratio=1.0`, no nearby
    surrogate, no nearby-kill term) -- the very capability this epic's
    off-lattice refactor exists to add (Task 2.3's uncapped recruitment).
    Because both drivers still consume the SAME infection/resistance/Allee/
    ROS mechanism (verbatim-ported into `EpitheliumProcess`, task-1.3-
    report.md) and the SAME immune contact-kill rate formula (`killing.
    contact_kill_rate`, reduced to a proximity gate off-lattice -- see
    immune_process.py's "Killing" docstring section), the epithelial fate
    cascade is expected to stay close even though immune COUNTS are not
    expected to match.

    ACHIEVED at these exact args (measured, not tuned): `ref` (run_full_
    model) dead[-1] = 0, `got` (run_full_model_composite) dead[-1] = 1 --
    delta = 1, well within tolerance (max(3, int(0.15*max(0,1))) = 3). Both
    drivers agree the reduced-scale/short-window run stays almost entirely
    in the infection/no-death regime (infected 6->5 vs 6->4, uninfected
    58->59 vs 58->59) -- no immune-driven divergence was large enough at
    this scale/window to need suppressing immune killing for a clean
    epithelial-only comparison; the small got>ref gap is consistent with the
    composite's immune agents landing a proximity kill (`immune_kills`,
    verified live by `test_run_full_model_composite_immune_killing_wired`
    above) that `run_full_model`'s reference run's NK/CD8 (locked at their
    seeded n=4 each, well below this run's recruit-reserve-capped ceiling)
    did not land in the same window -- an immune-driven epithelial
    difference, not a wiring bug (no driver/process change was made to
    reach this delta)."""
    from pbg_cpm_studies.influenza import run

    ref = run.run_full_model(cells_per_side=8, steps=15, seed=2, init_infection_frac=0.1)
    got = run.run_full_model_composite(cells_per_side=8, steps=15, seed=2, init_infection_frac=0.1)

    d_ref, d_got = ref["counts"]["dead"][-1], got["counts"]["dead"][-1]
    tol = max(3, int(0.15 * max(d_ref, 1)))
    assert abs(d_got - d_ref) <= tol, (
        f"epithelial dead-count parity broke: ref={d_ref} got={d_got} "
        f"(tol={tol})")


def test_composite_immune_not_pool_capped():
    """Task 3.4 -- the core proof this epic's off-lattice refactor exists for:
    `run_full_model`'s immune population is hard-capped by a fixed CPM
    reserve pool (`RECRUIT_RESERVE_TYPE`/`recruit_pool_per_type`, run.py's
    Task-8.4 module note -- the engine can only create CPM cells PRE-
    finalize, so ODE-driven inflow activates dormant reserve cells one at a
    time until the pool is exhausted, then inflow silently fails to seed).
    The composite's `ImmuneProcess` agents are OFF-LATTICE with no such pool
    (`local_ratio=1.0`, uncapped spawn -- immune_process.py's "Recruitment"
    docstring section): every ODE-predicted inflow draw mints a fresh agent,
    with no ceiling.

    Reduced-scale recruitment is documented elsewhere (`run_full_model`'s own
    docstring) as "near-inert" at fast-test scale -- the ODE's per-MCS
    inflow rate (`recruitment.macrophage_inflow`) is dominated by a
    homeostatic baseline term `mu_m*b_m` that scales with the tissue-cell
    count (`b_m` is proportional to `cells_per_side**2`, `mu_m` fixed), so
    demonstrating UNCAPPED growth (rather than merely non-zero growth) needs
    both a big enough `cells_per_side` for that baseline term to be
    non-negligible per-MCS AND enough total raw-MCS budget for it to
    accumulate. The brief's original cells_per_side=12/steps=40 example was
    measured and does not work: the macrophage count never left its seeded
    5-6 band (expected inflow over only 273 raw MCS at that scale is < 1
    agent). cells_per_side=24 (`mu_m*b_m ~= 0.0083/MCS`) was tried next at
    steps=175 (1044 raw MCS, seed=1) and only reached max=10 -- inside the
    Poisson-draw noise floor (outflow removals happened to outpace inflow
    for long stretches at this rate/seed), not a clean margin above
    `run_full_model`'s own default ceiling at this cell count
    (`n_macrophages=4` seed + `recruit_pool_per_type=6` reserve = 10 max).

    cells_per_side=48 (`tot_cell=2304`, `mu_m*b_m ~= 0.033/MCS` -- an order
    of magnitude stronger baseline signal, since `b_m` scales with
    `cells_per_side**2`) at steps=35 (238 raw MCS, seed=1,
    init_infection_frac=0.05) reliably clears the old ceiling with margin:
    macrophage count rises monotonically-ish from the seeded 6 to 14
    (measured series below; crosses 12 at record 28 / MCS~196), a clean,
    non-marginal double-digit excess over any small fixed-pool cap
    (`run_full_model`'s own default ceiling at ANY scale is
    n_macrophages+recruit_pool_per_type, a small constant independent of
    tissue size -- exactly the number this refactor removes). See
    task-3.4-report.md for the full scale/cost trade-off search (cells=24
    vs 48, steps vs wall-time) and the measured trajectory.

    Measured macrophage series (deterministic, seed=1) at these exact args:
    [6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8,
    8, 10, 10, 10, 12, 12, 12, 13, 14, 14, 14] -- max=14."""
    from pbg_cpm_studies.influenza import run

    got = run.run_full_model_composite(cells_per_side=48, steps=35, seed=1,
                                       init_infection_frac=0.05)
    macrophage = got["counts"]["macrophage"]
    # `run_full_model`'s own default ceiling (at ANY scale) would be
    # n_macrophages(4) + recruit_pool_per_type(6) = 10 -- grow comfortably
    # past that fixed-pool number, not a fixed "12" chosen without reference
    # to what the capped path could ever reach.
    assert max(macrophage) >= 12, (
        f"expected macrophage count to grow past the old fixed-pool ceiling "
        f"(run_full_model's default cap = 10); got max={max(macrophage)}, "
        f"series={macrophage}")
