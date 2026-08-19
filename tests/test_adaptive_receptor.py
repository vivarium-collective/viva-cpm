import process_bigraph as pb
from cpm.subcellular.adaptive_receptor import AdaptiveReceptorSubcell

def _proc(**over):
    cfg = {"kd": 2.9, "hill": 2.0, "conc_scale": 0.02, "activate_occupancy": 0.5,
           "epsilon": 0.1, "background": 0.0, "naive_type": 2, "activated_type": 3}
    cfg.update(over)
    return AdaptiveReceptorSubcell(cfg, core=pb.allocate_core())

def test_epsilon_zero_is_pure_hill():
    # with epsilon=0, m stays 0 so response == theta (the hill_occupancy rung / knockout)
    p = _proc(epsilon=0.0)
    theta, m_new, response = p.signal(ligand=200.0, m_prev=0.0)
    assert m_new == 0.0
    assert abs(response - theta) < 1e-9

def test_background_adds_to_ligand():
    # a uniform background raises sensed concentration -> higher occupancy
    p = _proc(epsilon=0.0, background=100.0)
    theta_bg, _, _ = p.signal(ligand=100.0, m_prev=0.0)
    theta_nobg, _, _ = _proc(epsilon=0.0, background=0.0).signal(ligand=100.0, m_prev=0.0)
    assert theta_bg > theta_nobg

def test_m_relaxes_toward_theta():
    # repeated updates at constant ligand drive m toward theta (adaptation)
    p = _proc(epsilon=0.3)
    theta, _, _ = p.signal(ligand=300.0, m_prev=0.0)
    m = 0.0
    for _ in range(50):
        _, m, _ = p.signal(ligand=300.0, m_prev=m)
    assert abs(m - theta) < 0.05          # adapted setpoint tracks occupancy
    _, _, response = p.signal(ligand=300.0, m_prev=m)
    assert response < 0.05                 # steady-state response adapts away

def test_update_emits_fate_and_m():
    p = _proc(epsilon=0.1)
    out = p.update({"ligand": 300.0, "m_prev": 0.0}, 1.0)
    assert set(out) == {"fate", "m"}
    assert out["fate"] in (2, 3)
