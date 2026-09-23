"""Task 8.4 tests: ODE-driven immune-cell recruitment (Hill inflow/outflow).

Pure-function tests (dossier Sec.4b(ii), discrepancy #12 -- the three laws are
asymmetric and must NOT be symmetrized) come first; a single fast, seeded
integration assertion on ``run.run_global_coupling`` follows.
"""
import math

import numpy as np

from viva_cpm_studies.influenza import params, price_ode, recruitment

P = params.load_params()
C = price_ode.resolve_constants(P["price_ode"], num_epithelial=900)


# --- pure rate functions --------------------------------------------------------

def test_macrophage_inflow_monotone_in_chemokine():
    lo = recruitment.macrophage_inflow(10.0, C)
    hi = recruitment.macrophage_inflow(1e5, C)
    assert hi > lo                                  # more chemokine -> more inflow
    assert recruitment.macrophage_inflow(0.0, C) == C["mu_m"] * C["b_m"]  # hill(0)=0 -> baseline only


def test_cd8_inflow_apc_driven_no_baseline():
    # CD8 inflow is P-driven and has NO homeostatic baseline (discrepancy #12):
    assert recruitment.cd8_inflow(0.0, C) == 0.0     # hill(0)=0, no +mu*b term
    assert recruitment.cd8_inflow(1e5, C) > 0.0


def test_nk_inflow_has_baseline():
    assert recruitment.nk_inflow(0.0, C) == C["mu_k"] * C["b_k"]


def test_outflow_terms_carry_resistance_load():
    # NK/CD8 outflow rise with the resistance-weighted infected load (G_ki/B_ei)
    assert recruitment.nk_outflow(0.5, C) > recruitment.nk_outflow(0.0, C)
    assert recruitment.cd8_outflow(0.5, C) > recruitment.cd8_outflow(0.0, C)
    assert math.isclose(recruitment.macrophage_outflow(C), C["mu_m"])


def test_hill_zero_at_nonpositive_and_saturates():
    assert recruitment.hill(0.0, 5.0, 3) == 0.0
    assert recruitment.hill(-1.0, 5.0, 3) == 0.0
    assert 0.0 < recruitment.hill(5.0, 5.0, 3) < 1.0
    assert recruitment.hill(1e9, 5.0, 3) > 0.999


def test_cd8_recruitment_field_is_apc_not_chemokine():
    # discrepancy #12 asymmetry made explicit: CD8 inflow ignores C entirely,
    # macrophage/NK inflow ignore P entirely.
    big = 1e5
    assert recruitment.cd8_inflow(big, C) > 0.0          # responds to P
    assert recruitment.cd8_inflow(0.0, C) == 0.0         # does NOT floor on a baseline
    # macrophage/NK have a strictly positive baseline even with no signal:
    assert recruitment.macrophage_inflow(0.0, C) > 0.0
    assert recruitment.nk_inflow(0.0, C) > 0.0


def test_ul_rate_to_prob_and_poisson_sampler_deterministic():
    assert recruitment.ul_rate_to_prob(0.0) == 0.0
    assert 0.0 < recruitment.ul_rate_to_prob(1.0) < 1.0
    # sampler is 0 for a tiny rate and deterministic under a seeded generator
    rng = np.random.default_rng(0)
    assert recruitment.poisson_inflow_count(0.0, rng) == 0
    a = [recruitment.poisson_inflow_count(2.0, np.random.default_rng(k)) for k in range(5)]
    b = [recruitment.poisson_inflow_count(2.0, np.random.default_rng(k)) for k in range(5)]
    assert a == b                                        # fully seeded/deterministic
    assert sum(a) > 0                                    # a non-trivial rate does seed sometimes


# --- fast seeded integration assertion (run_global_coupling) --------------------

def test_recruitment_persists_and_grows_populations():
    """Recruitment wired into run_global_coupling: (1) at the default tiny patch
    scale, macrophage/NK populations persist (do NOT collapse to zero) and the
    driving inflow rates are recorded and positive; (2) at a boosted ODE
    population scale (a documented *scenario* knob -- num_epithelial -- NOT a
    rate-constant tune), the (homeostatic-baseline-dominated) inflow reliably
    seeds new cells, so macrophage & NK counts end strictly higher WITH
    recruitment than with recruitment disabled.

    HONESTY (see task-8.4-report.md): the chemokine Hill term's contribution is
    eta-dependent because a_mc/a_kc scale with eta. At the tiny DEFAULT eta the
    dissociation constants (a_mc~0.08, a_kc~0.20) fall BELOW the chemokine field
    (C ~ 0.2-0.5), so the Hill term IS active and macro/NK inflow with signal
    exceeds the zero-signal baseline -- but the absolute rates are < 0.01/MCS,
    far below one cell per short run, so populations merely persist (no integer
    seeding). At the BOOSTED eta used below (eta=1: a_mc~575, a_kc~1356 >> C) the
    Hill term collapses to ~0 and inflow is baseline-dominated (mu*b ~ 3.6/MCS
    for macrophages) -- large enough to actually seed cells, but the growth is
    then HOMEOSTATIC, not chemokine-signal-driven. Neither regime yields
    chemokine-driven integer growth in a short run: a genuine field-vs-ODE scale
    calibration gap (same category as the Increment-7 chemotaxis-lambda finding),
    left for Increment 9. Constants are NOT tuned here.
    """
    from viva_cpm_studies.influenza import run

    # (1) default tiny-eta run: populations persist; rates recorded.
    base = run.run_global_coupling(side=30, steps=15, seed=5, with_immune=True,
                                   with_recruitment=True)
    assert set(base["recruit"]) >= {"M", "K", "E", "M_nb", "K_nb", "E_nb",
                                    "macro_inflow_rate", "nk_inflow_rate", "cd8_inflow_rate"}
    assert len(base["recruit"]["M"]) == len(base["mcs"]) == 15
    assert min(base["recruit"]["M"]) > 0            # macrophages never collapse
    assert min(base["recruit"]["K"]) > 0            # NK never collapse
    assert base["recruit"]["macro_inflow_rate"][0] > 0    # baseline inflow present
    assert base["recruit"]["nk_inflow_rate"][0] > 0

    # (2) boosted-eta: recruitment ON grows macrophage & NK vs recruitment OFF.
    kw = dict(side=30, steps=15, seed=5, with_immune=True, num_epithelial=250000,
              recruit_pool_per_type=50)
    on = run.run_global_coupling(with_recruitment=True, **kw)
    off = run.run_global_coupling(with_recruitment=False, **kw)
    m_on, k_on = on["recruit"]["M"][-1], on["recruit"]["K"][-1]
    m_off = off["spatial"]["M"][-1]
    k_off = off["spatial"]["K"][-1]
    assert m_on > m_off        # recruitment grows the macrophage population
    assert k_on > k_off        # ... and the NK population
    assert min(on["recruit"]["M"]) > 0 and min(on["recruit"]["K"]) > 0


def test_cd8_recruitment_is_apc_driven_asymmetry_preserved():
    """Discrepancy #12, at the driver level: CD8 recruitment is APC(P)-driven
    with NO homeostatic baseline. With APC present its inflow RATE is strictly
    positive every MCS; forcing the recruitment signal to zero (P->0) makes the
    CD8 inflow rate exactly 0 (no baseline to fall back on), while macrophage/NK
    keep their positive homeostatic baseline. This is the asymmetry a symmetric
    implementation would erase.
    """
    from viva_cpm_studies.influenza import run

    kw = dict(side=30, steps=12, seed=5, with_immune=True, with_recruitment=True)
    signal = run.run_global_coupling(**kw)
    control = run.run_global_coupling(recruit_zero_signal=True, **kw)

    # CD8 inflow responds to APC (P>0) and has no baseline:
    assert all(r > 0 for r in signal["recruit"]["cd8_inflow_rate"])
    assert all(r == 0.0 for r in control["recruit"]["cd8_inflow_rate"])
    # macrophage/NK retain a positive baseline even with the signal zeroed:
    assert all(r > 0 for r in control["recruit"]["macro_inflow_rate"])
    assert all(r > 0 for r in control["recruit"]["nk_inflow_rate"])
