# Known divergences from Sego 2022 (CC3D ViralInfectionVTM) — tracking register

Living register of where native viva-cpm diverges from the source, each triaged
by whether it is **fixable source-faithfully** (a real bug/mis-mapping to
correct), **irreducible** (source-faithful physics — closing it would require
tuning, which is excluded), or **investigate** (root not yet pinned). Source
authority: `docs/cc3d-reference/sego2022-*`. No constant is tuned to fit figures.

Status key: 🔴 fixable (real bug) · 🟡 investigate · 🟢 source-faithful/irreducible · ⚪ documented approximation (accepted).

| # | Divergence | Source expectation | Current viva-cpm | Root cause | Status | Notes |
|---|---|---|---|---|---|---|
| D1 | **Macrophage count ~2–4× below fig3b band** | fig3b M → 700 (lb 380) | peak ~119 (composite), ~176 (capped bootstrap) | chemokine field never saturates the recruitment Hill (a_mc≈2.82); M sits near ~105 homeostatic | 🟢 irreducible | Verified: not the pool cap (uncapping gives 119<176), IL-10 feedback z-invariant, a_11/a_12 correctly unscaled, z=2 won't move it (recruitment = z-invariant chemokine INTEGRAL). See `chemokine-recruitment-scale-analysis.md`. |
| D2 | **fig5 dose-response absent** | survival varies with viral load | uninfected→0 at EVERY load (1…10⁴) by 3.5 d | **virus field is missing the source's clearance terms** `−g_v·V/(1+a_v·V)` (≈0.137/MCS at low V) and `−g_vh·H·V`; our field has only the tiny linear `mu_v` (0.000286/MCS), ~350× too little → virus over-persists → over-infects | 🔴 **FIXABLE (confirmed)** | Probe: multiplying virus decay ×50–500 restores the dose-response exactly (load=1 un 6→137→144; load=100 →12 at ×500). `g_v`/`g_vh` are SOURCE constants in params.yaml, NEVER wired into the field. Fix = implement the source's virus-clearance reaction (Rust field-reaction term); source-faithful, NOT tuning. reference==composite, so it's a model-level fix. |
| D3 | **~0.5-day temporal lag (fig3b)** | infected peaks ~day 0.5 | infected peaks ~day 1.5 | early infection spread from scattered lesions too slow | 🟡 investigate | g_hv is a source value (speeding = fudging). Possible real cause: virus-field secretion/diffusion coupling vs source. Contrast with D2 (uniform-IC infection is *too* fast) — different ICs. |
| D4 | **Epithelial survivors when immunity present** | un@3.5d ~2 (near-total loss) | un@3.5d ~83–157 | stronger immunity clears infection → lowers ODE X → less ROS death → more survivors | 🟢 coupled to D1 | The paper's near-total loss corresponds to the (source-faithful, lower) immune level; not independently fixable. |
| D5 | **CD8 reaches band; NK/macro below** | M700/NK600/CD8550 (lb 380/280/280) | CD8 330 (in band), NK 186, M 119 | CD8 is APC(P)-driven (not chemokine-limited) → uncapping lifts it; NK/macro chemokine-limited (D1) | 🟢 partial success | The uncapping WORKED for the non-chemokine-limited type. |
| D6 | ROS death driven by global ODE scalar X | X is a global scalar (dossier §4b) | global ODE X | — | 🟢 source-faithful | Confirmed by ablation + dossier; NOT a divergence. |
| D7 | Intra-MCS parallel ordering (composite) | strict sequential per-MCS | epithelium ordered internally; immune kills merged at fixed priority; cross-process reads use MCS-start grid | pb has no within-MCS sequencing | ⚪ accepted | Second-order for low-prob per-MCS events; parity gate Δ=0 with immune off. |
| D8 | `init_viral_load` IC realization | source seeds Virus field uniformly to v0 | H cells as a transient virus source for one field advance | engine has no field-write primitive | ⚪ accepted | Documented in `_seed_uniform_virus`; may contribute to D2 (profile/magnitude). |
| D9 | Hybrid (not full 20-var) ODE | full coupled ODE+spatial | 10 systemic species integrated; spatial species are ODE inputs | discrepancy #9 | ⚪ accepted | Deliberate hybrid design. |
| D10 | resist multiplies kill rate DIRECTLY | discrepancy #7 | as source (direct, not complement) | source string | 🟢 source-faithful | Matches source literally. |
| D11 | spatial-string scalings (s_v/s_l/s_t) | discrepancy #10 | applied per dossier | — | 🟢 source-faithful | Verified for the chemokine/IL-10 chain (D1). |

## Triage & plan

- **🟢 irreducible (D1, D4, D5, D6, D10, D11):** these are source-faithful; closing them would require tuning (excluded) or geometry that's been ruled out (z=2 for D1). Track, don't "fix." D5 is a partial *success*.
- **⚪ accepted approximations (D7, D8, D9):** documented, second-order or by-design. Revisit D8 only if it proves to drive D2.
- **🟡 investigate (D2, D3):** the two with a plausible *real, fixable* root — both about infection dynamics. These are where "fix" effort should go, source-faithfully:
  - **D2 (fig5 dose-response):** check whether the source's viral-load IC is a *concentration* (as we seed) or a *count/multiplier* interpreted differently, and whether virus decay/clearance in the source prevents low-load self-amplification. If our IC magnitude or the virus-field decay differs from source, that's a genuine fix.
  - **D3 (fig3b lag):** audit the virus-field secretion (`secretion_g_vi=0.386` is 2× the source `g_vi=0.193` — flagged earlier; 2× high would make infection FASTER, so not obviously the lag cause, but worth a clean units re-audit) and the diffusion/decay vs source.

**Next fix target:** D2/D3 share the virus-field dynamics, so a focused audit of the virus field (IC magnitude, secretion units, decay/clearance) against the dossier is the highest-leverage source-faithful fix. Neither is confirmed a bug yet — the audit decides.
