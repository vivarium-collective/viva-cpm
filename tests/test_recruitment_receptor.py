import process_bigraph as pb
from pbg_cpm_studies.composites import chemotaxis_receptor as CR


def test_contact_j_invariant_activated_matches_naive():
    """Activated responder (3) shares every contact-J with naive (2): only
    chemotaxis differs, so adhesion is not confounded."""
    spec = CR.build_receptor_spec(cue_rate=10.0, chemo_lambda=14.0, blocked=False, seed=17)
    j = {(min(r["a"], r["b"]), max(r["a"], r["b"])): r["j"] for r in spec["contact"]}
    for other in (0, 1, 2, 3):
        lo2, hi2 = min(2, other), max(2, other)
        lo3, hi3 = min(3, other), max(3, other)
        assert j[(lo2, hi2)] == j[(lo3, hi3)], f"J mismatch for type {other}"


def test_only_activated_type_chemotaxes():
    spec = CR.build_receptor_spec(cue_rate=10.0, chemo_lambda=14.0, blocked=False, seed=17)
    chemo = {c["type"]: c["lambda"] for c in spec["fields"][0]["chemotaxis"]}
    assert chemo.get(CR.NAIVE_TYPE, 0.0) == 0.0
    assert chemo.get(CR.ACTIVATED_TYPE) == 14.0


def test_blocked_zeroes_activated_lambda():
    spec = CR.build_receptor_spec(cue_rate=10.0, chemo_lambda=14.0, blocked=True, seed=17)
    chemo = {c["type"]: c["lambda"] for c in spec["fields"][0]["chemotaxis"]}
    assert chemo.get(CR.ACTIVATED_TYPE, 0.0) == 0.0


def test_document_preinitializes_fates_for_all_responders():
    doc = CR.recruitment_receptor(cue_rate=10.0, chemo_lambda=14.0)
    n_responders = sum(1 for c in doc["cpm"]["config"]["spec"]["cells"]
                       if c["type"] == CR.NAIVE_TYPE)
    assert len(doc["fates"]) == n_responders
    assert all(v == CR.NAIVE_TYPE for v in doc["fates"].values())


def test_smoke_run_baseline_recruits_more_than_blocked():
    from pbg_cpm_studies.chemotaxis import metrics as M
    core = pb.allocate_core()

    def final_index(blocked):
        doc = CR.recruitment_receptor(core=core, blocked=blocked, seed=17)
        comp = pb.Composite({"state": doc}, core=core)
        comp.run(40)
        world = comp.state["cpm"]["instance"].world
        # "recruited" here means responders that have activated (type 3) AND
        # reached the source -- naive (2) responders never chemotax, so a
        # naive-type filter would read ~0 in both conditions and hide the
        # effect this composite is meant to exercise.
        return M.recruitment_index(world, responder_type=CR.ACTIVATED_TYPE)

    assert final_index(blocked=False) > final_index(blocked=True)
