"""Per-cell viral resistance from local IFN (Increment 3, Task 3.3).

Pure function matching the CC3D-source convention recorded in params.yaml's
`resistance:` block: `resist = f_bar/(a_rf+f_bar)`, a plain Hill/saturation
function of the local mean IFN concentration seen by a cell (source
`Type1InterferonModelSteppable.update_resistance`). This is NOT the paper's
printed `rho = 1 - f_bar'/(theta*a_f+f_bar')`; every consumer in the source
applies `(1 - resist)` where the paper's rho would appear in a protective
role, which is exactly how `run.run_virus_infection_with_ifn` uses this
value: gating virus release by `(1 - resist)`.
"""
from __future__ import annotations


def cell_resistance(ifn_at_cell: float, a_rf: float) -> float:
    """Cellular viral resistance from the local mean IFN field concentration.

    `resist = ifn_at_cell / (a_rf + ifn_at_cell)`, clamped so a negative
    (e.g. numerically-noisy) `ifn_at_cell` never yields a negative resist:
    0 at `ifn_at_cell == 0`, 0.5 at `ifn_at_cell == a_rf`, -> 1 as
    `ifn_at_cell -> inf`, monotone increasing throughout.
    """
    f_bar = max(0.0, ifn_at_cell)
    return f_bar / (a_rf + f_bar)
