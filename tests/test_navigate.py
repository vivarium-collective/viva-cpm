from pbg_cpm_studies.model_building import navigate


def test_climbs_static_to_hill_to_adaptive():
    a = [{"id": "receptor_gating", "severity": "hard", "verdict": "mismatch", "margin": -0.5},
         {"id": "recruits_high", "severity": "hard", "verdict": "mismatch", "margin": -0.3}]
    assert navigate.next_mechanism("static_lambda", a) == "hill_occupancy"
    a2 = [{"id": "receptor_gating", "severity": "hard", "verdict": "within_tol", "margin": 0.1},
          {"id": "recruits_high", "severity": "hard", "verdict": "mismatch", "margin": -0.3}]
    assert navigate.next_mechanism("hill_occupancy", a2) == "adaptive_receptor"
    a3 = [{"id": "recruits_high", "severity": "hard", "verdict": "within_tol", "margin": 0.1}]
    assert navigate.next_mechanism("adaptive_receptor", a3) is None


def test_no_stall_when_target_equals_active():
    # already hill but recruits_high still failing -> advance to adaptive, not stall on hill
    a = [{"id": "recruits_high", "severity": "hard", "verdict": "mismatch", "margin": -0.3}]
    assert navigate.next_mechanism("hill_occupancy", a) == "adaptive_receptor"
