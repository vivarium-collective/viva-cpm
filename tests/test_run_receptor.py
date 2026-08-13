"""Test as specified in task-4-brief.md Step 2, with two deliberate deviations
from the brief's illustrative literal snippet, both needed for the test to
actually exercise a real signal against the Task-3 default ``conc_scale``:

  * ``steps=60 -> steps=250``: at 60 composite ticks the secreted field has
    not yet diffused far enough for ANY responder (baseline or blocked) to
    cross even a low activation threshold -- the brief's 60 was an untested
    placeholder. 250 ticks gives baseline responders time to activate and
    reach the source while blocked stays at 0 across seeds 17/29/43.
  * ``out_dir=tmp_path``: without redirecting output, this test would
    overwrite the committed real receptor_<condition>.json summaries (same
    default path) on every pytest run.
"""
from pbg_cpm_studies.chemotaxis.run_receptor import run_receptor


def test_summary_shape_and_recruitment_ordering(tmp_path):
    base = run_receptor(blocked=False, seeds=[17, 29, 43], steps=250, kd=10.0,
                        out_dir=str(tmp_path))
    blk = run_receptor(blocked=True, seeds=[17, 29, 43], steps=250, kd=10.0,
                       out_dir=str(tmp_path))
    for s in (base, blk):
        assert len(s["mcs"]) == len(s["recruitment_index_mean"]) == len(s["recruitment_index_ci"])
        assert 0.0 <= s["final_mean"] <= 1.0
    assert base["final_mean"] > blk["final_mean"]
