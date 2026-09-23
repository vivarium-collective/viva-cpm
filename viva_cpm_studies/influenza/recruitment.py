"""Task 8.4: ODE-driven immune-cell RECRUITMENT (Hill inflow / outflow rates).

Resolves the recruitment stub carried since Increment 5. These are the exact
per-type inflow/outflow RATE functions of the Sego-2022 spatial-coupling model
(dossier ``docs/cc3d-reference/sego2022-global-ode.md`` Sec.4b(ii), verbatim from
``ViralInfectionVTMSteppables.ImmuneModelSteppable.inflow_rate_by_type`` /
``outflow_rate_by_type``, lines 1183-1203), with
``hill(v,a,h) = 1/(1+(a/v)^h)`` (``nCoVToolkit/nCoVUtils.hill_equation``; ``0`` at
``v<=0``).

**DISCREPANCY #12 -- the three laws are ASYMMETRIC; DO NOT symmetrize:**
  ============  =========================================  ==============  ========
  type          inflow rate                                outflow rate    field
  ============  =========================================  ==============  ========
  Macrophage    ``b_mc*hill(C,a_mc,h_m) + mu_m*b_m``        ``mu_m``        C (chemo)
  NK            ``b_kc*hill(C,a_kc,h_k) + mu_k*b_k``        ``G_ki + mu_k`` C (chemo)
  CD8+          ``b_ep*hill(P,a_ep,h_e)``                   ``B_ei + mu_e`` P (APC)
  ============  =========================================  ==============  ========

Two source-faithful asymmetries (dossier discrepancy #12, promoted from
``sego2022-parameters.md`` Sec.11):
  (a) CD8+ recruitment is driven by APC (``P``), NOT the chemokine field ``C``
      -- a DIFFERENT field from CD8+ chemotaxis (which is chemokine-driven).
  (b) CD8+ inflow has **NO** ``+ mu_e*b_e`` homeostatic baseline term (no
      ``hs_cd8`` constant exists), unlike macrophage/NK -- so ``cd8_inflow`` is
      exactly ``0`` when ``P <= 0``.

Constants come resolved (concrete floats, patch-cellularized) from
``price_ode.resolve_constants(...)``: ``b_mc,a_mc,mu_m,b_m,h_m`` (macrophage),
``b_kc,a_kc,mu_k,b_k,h_k`` (NK), ``b_ep,a_ep,h_e,mu_e`` (CD8+). ``G_ki``/``B_ei``
are the live resistance-weighted infected-load outflow terms pushed in per MCS
(dossier Sec.4a), passed by the caller.

All functions here are PURE (no RNG, no world mutation). The stochastic
rate->count conversion (``ul_rate_to_prob``, the Poisson-style inflow draw) and
the CPM seed/remove wiring live in ``run.run_global_coupling``; the seeded-RNG
sampler :func:`poisson_inflow_count` is provided here as a pure function of an
explicit ``Generator`` so it stays deterministic and testable.
"""
from __future__ import annotations

import math


def hill(v: float, a: float, h: float) -> float:
    """Hill function (``nCoVUtils.hill_equation``): ``0`` at ``v<=0``, else
    ``1/(1+(a/v)^h)``. Monotone increasing in ``v`` on ``v>0``, saturating to 1."""
    if v <= 0:
        return 0.0
    return 1.0 / (1.0 + (a / v) ** h)


# --- inflow rates (dossier Sec.4b(ii), :1183-1194) ------------------------------

def macrophage_inflow(C: float, consts: dict) -> float:
    """Macrophage inflow rate: ``b_mc*hill(C,a_mc,h_m) + mu_m*b_m``.

    Chemokine(``C``)-driven Hill term plus a homeostatic baseline
    ``mu_m*b_m`` (replenishment toward the tissue-resident population ``b_m``).
    ``hill(0)=0`` -> at ``C<=0`` the rate is exactly the baseline ``mu_m*b_m``.
    """
    return consts["b_mc"] * hill(C, consts["a_mc"], consts["h_m"]) + consts["mu_m"] * consts["b_m"]


def nk_inflow(C: float, consts: dict) -> float:
    """NK inflow rate: ``b_kc*hill(C,a_kc,h_k) + mu_k*b_k``.

    Chemokine(``C``)-driven Hill term plus homeostatic baseline ``mu_k*b_k``
    (same structural form as macrophage). ``hill(0)=0`` -> baseline only at
    ``C<=0``.
    """
    return consts["b_kc"] * hill(C, consts["a_kc"], consts["h_k"]) + consts["mu_k"] * consts["b_k"]


def cd8_inflow(P: float, consts: dict) -> float:
    """CD8+ inflow rate: ``b_ep*hill(P,a_ep,h_e)`` -- discrepancy #12.

    Driven by APC (``P``), NOT chemokine, and with **NO** homeostatic baseline
    term (no ``hs_cd8`` constant): the rate is exactly ``0`` when ``P<=0``. Do
    NOT add a ``+ mu_e*b_e`` term to symmetrize with macrophage/NK.
    """
    return consts["b_ep"] * hill(P, consts["a_ep"], consts["h_e"])


# --- outflow rates (dossier Sec.4b(ii), :1196-1203) -----------------------------

def macrophage_outflow(consts: dict) -> float:
    """Macrophage outflow (attrition) rate: ``mu_m`` (constant)."""
    return consts["mu_m"]


def nk_outflow(G_ki: float, consts: dict) -> float:
    """NK outflow rate: ``G_ki + mu_k``.

    Baseline attrition ``mu_k`` plus the live resistance-weighted infected-load
    term ``G_ki`` (dossier Sec.4a) -- NK are consumed faster where the infected
    load is higher.
    """
    return G_ki + consts["mu_k"]


def cd8_outflow(B_ei: float, consts: dict) -> float:
    """CD8+ outflow rate: ``B_ei + mu_e``.

    Baseline attrition ``mu_e`` plus the resistance-weighted infected-load term
    ``B_ei`` (the CD8+ analogue of NK's ``G_ki``).
    """
    return B_ei + consts["mu_e"]


# --- rate -> probability / count (source-faithful stochastic conversion) ---------

def ul_rate_to_prob(rate: float) -> float:
    """``ImmuneModelLib.ul_rate_to_prob``: ``1 - exp(-rate)`` (per-MCS unit-less
    rate -> Bernoulli probability). Used for outflow removals."""
    return 1.0 - math.exp(-rate)


def poisson_inflow_count(rate: float, rng) -> int:
    """Number of new cells to seed for a per-MCS inflow ``rate`` -- the source's
    ``inflow_by_type`` draw (``ViralInfectionVTMSteppables`` :1205-1212),
    reproduced verbatim over an explicit seeded ``numpy`` ``Generator`` for
    determinism.

    Each iteration draws a fresh uniform and accepts (increments the count) with
    probability ``1 - exp(-rate)*sum_{j=0}^{k} rate^j/j!`` = ``P(X >= k+1)`` for
    ``X ~ Poisson(rate)`` -- the source's sequential Poisson-tail sampler. For
    the tiny per-MCS rates at small patch scale this is ~always ``0``.
    """
    if rate <= 0.0:
        return 0
    num_add = 0
    exp_term = math.exp(-rate)
    sum_term = 1.0
    while rng.random() < 1.0 - exp_term * sum_term:
        num_add += 1
        sum_term += rate ** num_add / math.factorial(num_add)
    return num_add
