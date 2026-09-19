"""Task 4.3: the cellularized Allee effect -- uninfected cells poorly
surrounded by healthy tissue DIE (H -> D), dead cells well surrounded by
healthy tissue RECOVER (D -> H).

Pure functions matching the exact source expressions, `RecoverySteppable.step`
(see task-4.0-report.md "EXACT RecoverySteppable death + recovery
expressions" and params.yaml's `allee:` block for the full derivation and
provenance):

    srf_thresh = theta * srf_total
    death_rate   = b_h * (1 - resist) * srf_D          * (srf_thresh - srf_uninfected) / srf_total**2   (cell.type == H, only when srf_uninfected < srf_thresh)
    revive_rate  = b_h * (1 - resist) * srf_uninfected * (srf_uninfected - srf_thresh) / srf_total**2   (cell.type == D, only when srf_uninfected > srf_thresh)

`srf_total`/`srf_uninfected`/`srf_D` are per-cell contact-surface areas
restricted to EPITHELIAL neighbor types only (H, I, D -- the source's
`ec_list = [UNINFECTED, INFECTED, INFECTEDRELEASING, DYING]`; MEDIUM and any
immune-type contact is excluded from `srf_total`, matching the source's own
neighbor-type filter in `RecoverySteppable.step`). `b_h` is shared by both
branches (see params.yaml's `allee.a_H_equals_b_h` note -- no separate
recovery-side coefficient exists in the source). `Pr = 1 - exp(-rate)`
(`ImmuneModelLib.ul_rate_to_prob`) is applied by the caller, not here.

Both functions return 0.0 for an isolated cell (`srf_total == 0`, guarding
the `/srf_total**2` division) and 0.0 when the branch's firing condition
isn't met (rather than a negative rate), matching the source's `if`/`elif`
gating -- so callers can call both unconditionally and simply gate the
STOCHASTIC transition on `rate > 0` (or just draw against `Pr = 1 -
exp(-rate)`, which is 0 exactly when `rate == 0`).
"""
from __future__ import annotations


def allee_death_rate(srf_D: float, srf_uninfected: float, srf_total: float,
                      b_h: float, resist: float, theta: float) -> float:
    """Allee death rate for an UNINFECTED (H) cell.

    0 when `srf_total == 0` (isolated cell) or `srf_uninfected >=
    srf_thresh` (well enough surrounded by healthy tissue not to be at
    risk); otherwise the source-literal rate, positive and monotonically
    increasing as the surface deficit `(srf_thresh - srf_uninfected)`
    grows (i.e. the cell is more poorly surrounded by uninfected
    neighbors), and increasing in `srf_D` (more dead-cell contact).
    """
    if srf_total == 0:
        return 0.0
    srf_thresh = theta * srf_total
    if srf_uninfected >= srf_thresh:
        return 0.0
    return b_h * (1.0 - resist) * srf_D * (srf_thresh - srf_uninfected) / srf_total ** 2


def allee_recovery_rate(srf_uninfected: float, srf_total: float,
                         b_h: float, resist: float, theta: float) -> float:
    """Allee recovery rate for a DYING (D) cell.

    0 when `srf_total == 0` (isolated cell) or `srf_uninfected <=
    srf_thresh` (not well enough surrounded by healthy tissue to revive);
    otherwise the source-literal rate, positive and increasing as
    `srf_uninfected` rises above `srf_thresh`. Note: for a DYING cell the
    source's `srf_area_oi` (interface area with UNINFECTED neighbors) is
    the SAME quantity as `srf_uninfected`, so `srf_uninfected` appears
    twice in this formula (once as the `srf_area_oi` factor, once inside
    the `(srf_uninfected - srf_thresh)` deficit) -- see task-4.0-report.md.
    """
    if srf_total == 0:
        return 0.0
    srf_thresh = theta * srf_total
    if srf_uninfected <= srf_thresh:
        return 0.0
    return b_h * (1.0 - resist) * srf_uninfected * (srf_uninfected - srf_thresh) / srf_total ** 2
