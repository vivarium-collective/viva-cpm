from pathlib import Path
import yaml

BASE = Path("workspace/studies/recruitment-receptor-baseline/study.yaml")


def _all_cites(node):
    found = set()
    if isinstance(node, dict):
        for k, v in node.items():
            if k == "cites" and isinstance(v, list):
                found |= {str(x) for x in v}
            found |= _all_cites(v)
    elif isinstance(node, list):
        for x in node:
            found |= _all_cites(x)
    return found


def test_baseline_wired_and_cites_kd_source():
    doc = yaml.safe_load(BASE.read_text())
    assert doc["baseline"][0]["composite"].endswith("chemotaxis_receptor.recruitment_receptor")
    assert float(doc["baseline"][0]["params"]["kd"]) == 2.9   # cited receptor Kd
    # the cited Kd source must be linked somewhere in the study (model_settings and/or the test)
    assert "nasser2009" in _all_cites(doc), "study must cite the receptor-Kd source (nasser2009)"
