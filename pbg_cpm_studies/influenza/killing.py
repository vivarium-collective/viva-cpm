"""Task 7.2: NK/CD8 contact-mediated killing of infected cells -- the LOCAL
(surface-contact) term of `ContactKillingSteppable` (Simulation/
ViralInfectionVTMSteppables.py, ~lines 212-319; see params.yaml's `nk.killing`/
`cd8.killing` blocks for the full derivation/provenance).

Pure function matching the source-literal expression:

    kill_rate = g_i * tot_ec_ODE * srf_immune * cell_resist / cell.volume

(`g_i` = `nk.g_ik` or `cd8.g_ie`; `tot_ec_ODE` = `scaling.ode_epithelial_
population`, 250000; `srf_immune` = the infected cell's contact-surface area
with NK (`srf_nk`, `types.K`) or CD8+ (`srf_cd8`, `types.E`) neighbors, from
`World.cell_contact_area_by_type` -- the same Increment-4 primitive `allee.py`'s
callers use for `srf_uninfected`/`srf_D`/`srf_total`). The caller (`run.
run_cytotoxic_response`'s killing step) applies `Pr = 1 - exp(-kill_rate)`
(`ImmuneModelLib.ul_rate_to_prob`) and, separately, draws NK and CD8 kill rolls
independently before deciding the cell's fate -- not this module's concern.

**DISCREPANCY #7 (source-faithful, NOT "fixed" here -- see params.yaml's
`nk.killing.kill_rate_resist_direct`/`cd8.killing.kill_rate_resist_direct` and
docs/cc3d-reference/sego2022-parameters.md §7.7):** `ContactKillingSteppable.
step()` multiplies `kill_rate` by `cell_resist` DIRECTLY (`cell.dict[
ImmuneModelLib.im_resist_key]`, read as-is) -- NOT `(1 - cell_resist)` -- unlike
EVERY other resistance consumer in this codebase (`resistance.cell_resistance`'s
module docstring; `run.run_virus_infection_with_ifn`/`run_epithelial_fate`'s
`(1.0 - resist)` secretion-scale gate; `transitions.infected_death_step`'s
`mu_i * (1 - resist)`; `allee.allee_death_rate`/`allee_recovery_rate`'s
`b_h * (1 - resist) * ...`). Under the (1-resist) convention, `resist` is
PROTECTIVE (higher local IFN -> lower rate); here, taken literally, a MORE
"resistant" cell is killed FASTER by contact immune cells -- either an
intentional antigen-presentation-style effect in the source (high local
interferon marking a cell as a stronger target for NK/CD8 recognition) or a
source inconsistency. This function implements the literal source form
(`cell_resist` direct) per the brief -- it is NOT silently "fixed" to
`(1 - cell_resist)`. Flagged again in `run.run_cytotoxic_response`'s docstring
and in task-7.2-report.md.
"""
from __future__ import annotations


def contact_kill_rate(srf_immune: float, cell_resist: float, g_i: float,
                       tot_ec: float, cell_volume: float) -> float:
    """Local contact-kill rate for one infected cell against one immune type
    (NK or CD8+): `g_i * tot_ec * srf_immune * cell_resist / cell_volume`.

    `srf_immune`: contact-surface area (lattice-site units) between the
    infected cell and same-type immune-cell (NK or CD8+) neighbors, e.g. from
    `world.cell_contact_area_by_type(cid).get(types.K, 0)` /
    `.get(types.E, 0)`.
    `cell_resist`: this cell's local resistance (`resistance.cell_resistance`),
    applied DIRECTLY -- see this module's docstring, DISCREPANCY #7.
    `g_i`: the per-MCS killing-rate coefficient (`params.yaml`'s `nk.g_ik` or
    `cd8.g_ie`).
    `tot_ec`: the ODE-scale epithelial population constant (`params.yaml`'s
    `scaling.ode_epithelial_population`, 250000).
    `cell_volume`: the infected cell's current volume (lattice sites, e.g.
    `world.cell_volumes()[cid]`).

    Always >= 0 (every factor is non-negative for valid inputs). 0 when
    `srf_immune == 0` (no contact with that immune type) or `cell_volume == 0`
    (degenerate/absent cell -- division guarded rather than raising, matching
    `allee.py`'s `srf_total == 0` convention). Monotonically non-decreasing in
    both `srf_immune` and `cell_resist` (all other factors held fixed and
    non-negative) -- more contact area or higher (direct) resist means a
    higher kill rate.
    """
    if srf_immune <= 0 or cell_volume <= 0:
        return 0.0
    return g_i * tot_ec * srf_immune * cell_resist / cell_volume
