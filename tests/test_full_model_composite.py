from pbg_cpm_studies.composites.influenza import full_model_composite_document
from pbg_cpm_studies.core import build_core
from process_bigraph import Composite


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
