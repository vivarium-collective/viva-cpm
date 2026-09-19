"""Task 6.1: per-cell secretion-scale factors for the macrophage-released
chemokine + IL-10 fields (Increment 6 crux).

Two pure functions, matching the CC3D-source secretion-regulation forms
recorded in `params.yaml`'s `il10:` block (`ChemokineSecretionSteppable` /
`IL10SecretionSteppable`, Simulation/ViralInfectionVTMSteppables.py):

  - `macrophage_secretion_scale`: the shared IL-10-Hill self-regulation
    factor that gates BOTH macrophage chemokine secretion (`sec = b_c *
    Hill(L_loc)`) and macrophage IL-10 secretion (`sec = b_l * Hill(L_loc)`)
    -- `Hill(L) = sig_1 / (sig_1 + (g_1*L + g_2) / (L + d_2))`, L = the local
    mean IL-10 concentration seen by the macrophage cell. See `fields.py`'s
    `add_chemokine_field`/`add_il10_field` for how the base per-pixel rate
    (`_per_pixel_secretion_rate`) combines with this per-cell scale via
    `world.set_cell_secretion_scale` (the Increment-3 per-cell-scale
    primitive).
  - `uninfected_il10_scale`: the uninfected (H) cell IL-10 secretion gate,
    `1 - resist` (source: `sec_amount = mu_l * b_lh * (1 - resist)`).

**sig_1 STUB (ruling, see params.yaml's `il10.sig_1_stub` comment):** the
source's `sig_1 := a_11*T + a_12*D` is a live Antimony assignment-rule read
from the immune-model ODE solver (T = TNF-like population, D = tot_cell-H-I),
which no increment through 8 implements. Callers pass `params.il10.
sig_1_stub` (a documented constant, NOT the true dynamic value) as this
function's `sig_1` argument. This fixes the TNF-dependence pending
Increment 8; the L-dependent Hill self-regulation shape captured here, and
the resulting chemokine gradient shape (set by diffusion from macrophage
locations, not by sig_1's exact value), are preserved regardless.
"""
from __future__ import annotations


def macrophage_secretion_scale(il10_local: float, sig_1: float, g_1: float,
                                g_2: float, d_2: float) -> float:
    """The shared IL-10-Hill secretion-scale factor gating macrophage
    chemokine + IL-10 secretion: `sig_1 / (sig_1 + (g_1*L + g_2) / (L + d_2))`,
    `L` = `il10_local` (clamped to >= 0 so a negative, e.g. numerically-noisy,
    local IL-10 reading never flips the ratio term negative).

    Monotone DECREASING in `il10_local` for `g_1, d_2 > 0`: more local IL-10
    -> a larger `(g_1*L+g_2)/(L+d_2)` ratio term -> a smaller Hill factor
    (self-limiting secretion as IL-10 builds up). Always >= 0 for `sig_1,
    g_1, g_2, d_2 > 0` and `il10_local >= 0` (a ratio of two non-negative
    quantities); returns 0.0 in the degenerate case where the denominator is
    not positive (e.g. `sig_1 <= 0`), rather than raising or returning a
    negative/undefined value.
    """
    L = max(0.0, il10_local)
    denom = sig_1 + (g_1 * L + g_2) / (L + d_2)
    if denom <= 0.0:
        return 0.0
    return max(0.0, sig_1 / denom)


def uninfected_il10_scale(resist: float) -> float:
    """Uninfected (H) cell IL-10 secretion scale: `1 - resist`, clamped to
    [0, 1] (source: `sec_amount = mu_l * b_lh * (1 - resist)`, applied per
    UNINFECTED cell -- see `params.yaml`'s `il10.mu_l`/`il10.b_lh`).
    `resistance.cell_resistance` already returns values in [0, 1), but this
    function clamps independently so it stays well-defined for any caller."""
    return min(1.0, max(0.0, 1.0 - resist))
