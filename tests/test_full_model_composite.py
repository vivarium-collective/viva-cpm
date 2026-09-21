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


def test_full_model_composite_epithelium_enables_chemokine_il10_gate():
    """Source-faithfulness fix (final review round): the composite's
    epithelium node previously carried no `enable` config at all, silently
    falling to `EpitheliumProcess._DEFAULT_ENABLE`, which OMITS "chemokine"
    -- so the uninfected-H IL-10 gate (`EpitheliumProcess._epithelial_fates`,
    gated by `"chemokine" in self.enable`) was always off, unlike
    `run_full_model`'s default (`run.py`'s `_FULL_MODEL_SUBSYSTEMS`, which
    includes "chemokine"). Proves the fix at both the document level (the
    config the composite ships) and the constructed-instance level (what
    `EpitheliumProcess.initialize` actually resolves it to)."""
    from pbg_cpm_studies.composites.influenza import FULL_MODEL_EPITHELIUM_ENABLE

    assert "chemokine" in FULL_MODEL_EPITHELIUM_ENABLE

    doc = full_model_composite_document(cells_per_side=6, seed=1, init_infection_frac=0.1)
    assert "chemokine" in doc["epithelium"]["config"]["enable"]

    core = build_core()
    comp = Composite({"state": doc}, core=core)
    epi_proc = comp.state["epithelium"]["instance"]
    assert "chemokine" in epi_proc.enable
    # And the rest of the epithelial-fate pipeline stays active alongside it.
    for token in ("infection", "ifn", "death", "ros", "allee", "killing"):
        assert token in epi_proc.enable


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
    """Task 3.4 (fix round 1, review finding): the ORIGINAL version of this
    test (`cells_per_side=8, steps=15, seed=2, init_infection_frac=0.1`) drove
    `dead_ref` to exactly 0, collapsing the tolerance to its floor
    (`max(3, int(0.15*max(0,1)))=3`) and passing for ANY `dead_composite` in
    [0,3] -- not a real parity check. Fixed by moving to
    `cells_per_side=8, seed=2, init_infection_frac=0.3, steps=30` (reviewer-
    verified: gives `run_full_model` a NON-TRIVIAL, fluctuating dead series
    peaking ~4), comparing dead AND uninfected ELEMENT-WISE across the full
    30-record series (not just the terminal value), each record's tolerance
    computed from the plan's formula (`max(3, int(0.15*max(dead_ref[i],1)))`
    applied per record for dead; same form for uninfected).

    **Debugging the divergence (per the fix instructions):** at this
    regime's FULL default `enable` (both drivers' own defaults --
    `run.run_full_model`'s 11-token `_FULL_MODEL_SUBSYSTEMS`, the composite's
    always-active pipeline), `ref` dead[-1]=1 (series peaks at 4, matching
    the reviewer's report) but `got` (`run_full_model_composite`, unmodified)
    dead[-1]=8 -- delta 7, WAY beyond tolerance (3). Traced the cause (not a
    composite wiring bug):

    `run_full_model`'s `_one_mcs` only calls `world.set_cell_secretion_scale`
    for macrophage cells INSIDE its "chemokine" block, itself inside the
    per-MCS loop -- meaning the CPM engine's per-cell secretion-scale default
    (`Field::cell_secretion_scale`, `crates/cpm-core/src/field.rs`: `1.0`
    when never explicitly set for a cell) is still in effect for the VERY
    FIRST `world.step(1)` call of the whole run, since that call happens
    BEFORE `_one_mcs`'s step-4 secretion-scale block ever runs once. So
    macrophages secrete chemokine at the FULL, UNTHROTTLED rate for exactly
    one MCS (worth `n_macrophages * chemokine.b_c/2` -- confirmed
    numerically: 4 macrophages * 0.0279 = 0.112, matching the measured
    `ode.C` input at record 1, 0.103, within diffusion/decay rounding) before
    the sig_1-Hill gate (correctly computing scale=0 while sig_1's IC is 0.0)
    ever gets applied on MCS 2+. That one-time accidental burst then DOMINATES
    the entire run's chemokine trajectory (ongoing sig_1-throttled secretion
    is ~1000x smaller by comparison, matching the loop-gain analysis in
    `docs/cc3d-reference/chemokine-recruitment-scale-analysis.md`), driving
    `ref`'s C/T/N/X ODE trajectory and hence its ROS-driven death differently
    from the composite's `ImmuneProcess`, which computes `scale` fresh from
    the actual `sig_1` INPUT value every update (no CPM-engine default-scale
    concept -- it never uses `set_cell_secretion_scale` for macrophages at
    all, secreting via `field_deposit` point sources instead) and so never
    gets this initialization-order artifact. This lives entirely in
    `run_full_model`'s own per-MCS loop ordering (`run.py`, the protected
    source-faithful reference) -- NOT touched, per the task's explicit "do
    NOT touch run_full_model" constraint, and NOT a composite bug: verified
    by isolating pure infection/ifn/death/allee (no ros/killing/chemokine) on
    BOTH drivers (below) and finding an EXACT match, proving the shared
    cascade is correctly ported; the divergence is 100% attributable to the
    ROS/ODE coupling amplifying this reference-side secretion-timing quirk.

    **Fallback used (reviewer-sanctioned):** compare with immune/ODE
    influence on epithelial fate minimized on BOTH sides -- `enable=
    {"infection","ifn","death","allee"}` (drops "ros" so the ODE's X/ROS
    coupling, where the artifact surfaces, cannot affect either series;
    drops "killing"/"chemokine"/"macrophage"/"nk_cd8"/"recruitment"/"ode"
    tokens, which are either not epithelium-relevant or not gate-able
    per-driver without touching source -- `run_full_model` accepts `enable`
    directly; the composite side sets `EpitheliumProcess.enable` directly on
    the constructed instance, since the document's `config["enable"]`
    override is silently ignored by `bigraph_schema`'s `core.fill` for a
    plain "list"-typed config field -- it APPENDS the override onto the
    schema default instead of replacing it (confirmed: `core.fill({"_type":
    "list", "_default": [...]}, ["x"])` returns default+["x"], not ["x"]) --
    a separate, pre-existing `bigraph_schema`/`EpitheliumProcess.config_schema`
    quirk, out of this fix round's scope, flagged in task-3.4-report.md, not
    fixed here since it isn't exercised by `run_full_model_composite`'s own
    default call path (which never sets `enable` in the document at all)).

    ACHIEVED at this minimized-influence regime, ELEMENT-WISE across all 30
    records: `dead` series EXACT MATCH (max abs delta = 0,
    `[0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]` on both
    sides); `uninfected` series EXACT MATCH (max abs delta = 0,
    `[45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,44,44,
    43,43,43,43,43,42,41]` on both sides). Not a degenerate all-zero
    baseline (`dead_ref` is nonzero from record 8 on; `uninfected_ref`
    genuinely declines 45->41), so this is a real, non-trivial, element-wise
    parity check, not a "didn't explode" check."""
    from pbg_cpm_studies.influenza import run
    from pbg_cpm_studies.composites.influenza import full_model_composite_document
    from pbg_cpm_studies.core import build_core
    from pbg_cpm_studies.influenza import types as inf_types

    minimal_enable = {"infection", "ifn", "death", "allee"}
    cells_per_side, steps, seed, init_infection_frac = 8, 30, 2, 0.3
    mcs_per_step = 7

    ref = run.run_full_model(cells_per_side=cells_per_side, steps=steps, seed=seed,
                             init_infection_frac=init_infection_frac, enable=minimal_enable)

    # Drive the composite with the SAME cadence `run_full_model_composite`
    # uses (epithelium+immune every raw MCS, ODE once per mcs_per_step-MCS
    # record), but with `EpitheliumProcess.enable` set directly on the
    # constructed instance to the same minimized set (see docstring: the
    # document's config["enable"] override is silently ignored).
    doc = full_model_composite_document(cells_per_side=cells_per_side, seed=seed,
                                        init_infection_frac=init_infection_frac,
                                        mcs_per_step=mcs_per_step)
    doc["epithelium"]["config"]["mcs_per_step"] = 1
    core = build_core()
    comp = Composite({"state": doc}, core=core)
    epi_proc = comp.state["epithelium"]["instance"]
    epi_proc.enable = set(minimal_enable)
    imm_proc = comp.state["immune"]["instance"]
    ode_proc = comp.state["ode"]["instance"]

    store = {"fates": {}, "field_deposit": [], "ode_state": {}, "immune_kills": [],
             "immune_agents": list(doc["immune_agents"]), "chemo_field": [], "virus_field": [],
             "dims": [], "ode_inputs": {}, "immune_counts": {}, "recruit_drivers": {},
             "sig_1": 0.0, "margin_box": doc["margin_box"], "apply_recruitment": True}

    def _epithelial_counts():
        types_now = epi_proc.world.cell_types()
        H = sum(1 for t in types_now[1:] if t == inf_types.H)
        D = sum(1 for t in types_now[1:] if t == inf_types.D)
        return H, D

    got_dead, got_uninfected = [], []
    H, D = _epithelial_counts()
    got_dead.append(D)
    got_uninfected.append(H)

    for _i in range(1, steps):
        for _mcs in range(mcs_per_step):
            epi_out = epi_proc.update({
                "fates": store["fates"], "field_deposit": store["field_deposit"],
                "ode_state": store["ode_state"], "immune_kills": store["immune_kills"],
            }, 1)
            store["chemo_field"] = epi_out["chemo_field"]
            store["virus_field"] = epi_out["virus_field"]
            store["dims"] = epi_out["dims"]
            store["ode_inputs"] = epi_out["ode_inputs"]
            types_now = epi_out["types"]
            positions_now = epi_out["positions"]
            infected_ids = [cid for cid in range(1, len(types_now)) if types_now[cid] == inf_types.I]
            epithelial_positions = {str(cid): [float(positions_now[cid][0]), float(positions_now[cid][1])]
                                    for cid in infected_ids}
            imm_out = imm_proc.update({
                "chemo_field": store["chemo_field"], "virus_field": store["virus_field"],
                "dims": store["dims"], "immune_agents": store["immune_agents"],
                "epithelial_positions": epithelial_positions, "infected_ids": infected_ids,
                "sig_1": store["sig_1"], "recruit_drivers": store["recruit_drivers"],
                "margin_box": store["margin_box"], "apply_recruitment": store["apply_recruitment"],
            }, 1)
            store["immune_agents"] = imm_out["immune_agents"]
            store["immune_kills"] = imm_out["immune_kills"]
            store["field_deposit"] = imm_out["field_deposit"]
            store["immune_counts"] = imm_out["immune_counts"]

        ode_inputs = store["ode_inputs"]
        ode_out = ode_proc.update({
            "H": ode_inputs.get("H", 0.0), "I": ode_inputs.get("I", 0.0), "DH": ode_inputs.get("DH", 0.0),
            "V": ode_inputs.get("V", 0.0), "F": ode_inputs.get("F", 0.0), "C": ode_inputs.get("C", 0.0),
            "L": ode_inputs.get("L", 0.0), "B_ei": ode_inputs.get("B_ei", 0.0), "G_ki": ode_inputs.get("G_ki", 0.0),
            "M": store["immune_counts"].get("M", 0), "K": store["immune_counts"].get("K", 0),
            "E": store["immune_counts"].get("E", 0),
        }, mcs_per_step)
        store["ode_state"] = ode_out["ode_state"]
        store["recruit_drivers"] = ode_out["recruit_drivers"]
        store["sig_1"] = ode_out["sig_1"]

        H, D = _epithelial_counts()
        got_dead.append(D)
        got_uninfected.append(H)

    ref_dead = ref["counts"]["dead"]
    ref_uninfected = ref["counts"]["uninfected"]
    assert len(ref_dead) == len(got_dead) == steps

    dead_deltas = []
    uninfected_deltas = []
    for i in range(steps):
        d_ref, d_got = ref_dead[i], got_dead[i]
        u_ref, u_got = ref_uninfected[i], got_uninfected[i]
        dead_tol = max(3, int(0.15 * max(d_ref, 1)))
        uninfected_tol = max(3, int(0.15 * max(u_ref, 1)))
        dead_deltas.append(abs(d_got - d_ref))
        uninfected_deltas.append(abs(u_got - u_ref))
        assert abs(d_got - d_ref) <= dead_tol, (
            f"record {i}: epithelial dead-count parity broke: "
            f"ref={d_ref} got={d_got} (tol={dead_tol})")
        assert abs(u_got - u_ref) <= uninfected_tol, (
            f"record {i}: epithelial uninfected-count parity broke: "
            f"ref={u_ref} got={u_got} (tol={uninfected_tol})")

    # Not a degenerate always-0 baseline (the original test's failure mode).
    assert max(ref_dead) > 0
    assert ref_uninfected[-1] < ref_uninfected[0]
    # ACHIEVED (measured): both max deltas are 0 -- exact element-wise match.
    assert max(dead_deltas) == 0
    assert max(uninfected_deltas) == 0


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


def test_run_full_model_composite_init_viral_load_seeds_infection():
    """fig5/fig7 path: init_viral_load lays down a ~uniform virus-field IC with
    NO pre-infected cells; infection then EMERGES from the field (mirrors
    run_full_model's init_viral_load branch, now supported by the composite)."""
    from pbg_cpm_studies.influenza import run
    r = run.run_full_model_composite(cells_per_side=8, steps=18, seed=1,
                                     init_viral_load=1000.0)
    assert r["counts"]["infected"][0] == 0        # no pre-infected cells at t0
    assert max(r["counts"]["infected"]) > 0       # infection emerges from the field
    assert r["params"]["init_viral_load"] == 1000.0


def test_repro_fig5_runs_on_composite_engine():
    """The fig5 viral-load sweep runs end-to-end through the pb-Composite
    (engine='composite') now that init_viral_load is wired."""
    from pbg_cpm_studies.influenza import run
    r = run.repro_fig5(loads=(1, 10000), replicas=1, cells_per_side=8, steps=10,
                       seed0=0, engine="composite")
    assert set(r["by_load"].keys()) == {1, 10000}
    surv = lambda L: r["by_load"][L]["uninfected_final_frac"]
    assert surv(10000) <= surv(1) + 1e-9      # higher load => not more survivors
    assert "band_eval" in r


def test_repro_fig7_runs_on_composite_engine():
    """The fig7 infection-fraction sweep runs end-to-end through the composite."""
    from pbg_cpm_studies.influenza import run
    r = run.repro_fig7(fracs=(0.001, 0.05), replicas=1, cells_per_side=8, steps=10,
                       seed0=0, engine="composite")
    assert set(r["by_frac"].keys()) == {0.001, 0.05}
    assert "band_eval" in r
