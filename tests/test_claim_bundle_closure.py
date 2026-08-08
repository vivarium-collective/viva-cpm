"""The chemotactic-recruitment claim bundle derives the expected closure verdict.

Guards the semantic-layer target for the reference slice: the bundle is
well-formed and backend-realizable, its studies executed and gated, but the
claim is deliberately NOT complete — an unresolved lambda-calibration hole keeps
the prediction qualitative. Verdict is derived by code, never hand-set.

The deriver is pure YAML logic, so we load it by file path to keep this test
decoupled from the Rust engine (importing the package would pull in composites).
"""

import importlib.util
from pathlib import Path

_MODULE = Path(__file__).resolve().parents[1] / "pbg_cpm_studies" / "claim_bundle_closure.py"


def _derive():
    spec = importlib.util.spec_from_file_location("claim_bundle_closure", _MODULE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.derive_closure()


def test_chemotactic_recruitment_closure():
    c = _derive()

    assert c["well_formed"] is True, c.get("problems")
    assert c["backend_realizable"] is True, c.get("missing_studies")
    assert c["execution_ground"] is True, c.get("ungated_studies")
    assert c["evaluated_by"] == "code"

    # Deliberately incomplete: the lambda-calibration hole is required_for a
    # quantitative prediction, so the claim is demonstrated only qualitatively.
    assert c["claim_complete"] is False
    hole_ids = [h["id"] for h in c["unresolved"]]
    assert "hole:lambda-calibration" in hole_ids
