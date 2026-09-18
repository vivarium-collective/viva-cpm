# Influenza (Sego et al. 2022) reproduction — design

**Status:** design (approved in brainstorming 2026-09-18)
**Workspace:** viva-cpm (`cpm-studies`)
**Author:** Eran Agmon (with Claude)
**Reference:** `sego2022influenza` — Sego, Mochan, Ermentrout, Glazier (2022),
*J. Theor. Biol.* 532:110918. PDF at `workspace/references/papers/sego-2022-influenza.pdf`.

## 1. Goal

Reproduce, natively in viva-cpm, the multiscale multicellular spatiotemporal
model of local influenza A infection and immune response of Sego et al. (2022).
That paper *cellularizes* the Price et al. (2015) non-spatial in-host influenza
ODE model into a CompuCell3D (CC3D) hybrid CPM + reaction–diffusion + global-ODE
model. We do **not** bridge CC3D; we build viva-cpm's own engine up until it can
run the same model and match the paper's published results.

**Fidelity target (decided in brainstorming):** *quantitative figure match at
paper scale* — reproduce Figs 3B, 5, and 7 numerically at/near paper scale
(1 mm × 1 mm sheets, ~250k epithelial cells, 50 replicas per condition) within
acceptance bands.

**Oracle (decided):** the paper's supplementary **CompuCell3D source** is the
authoritative parameter/spec; the **published figures** are the numeric
acceptance targets. (No live CC3D co-simulation required.)

**Breadth (decided):** the *full* model this round — all four diffusive fields
(virus, type-I IFN, chemokines, IL-10), all local immune types (macrophage, NK,
CD8⁺ T), the cellularized Allee recovery, and the full Price-2015 global ODE
compartment; all initial-condition scenarios (viral-load and infection-fraction
sweeps).

## 2. What the paper's model is (condensed)

A quasi-2D planar epithelial sheet (5×5-site cells on a 2 µm lattice) evolved by
the CPM (Glazier–Graner–Hogeweg): effective energy `H` = volume constraint +
adhesion (contact energies `J`) + chemotaxis; immune cells additionally
haptotax. Coupled to:

- **Four local diffusive fields** (Table 1 PDEs): extracellular virus `v`,
  type-I IFN `f`, chemokines `c`, IL-10 `l`. Each: diffusion + decay + structured
  source/sink terms regulated by cell state, other fields, and global variables.
- **Epithelial type transitions** (Table 2), each a stochastic per-cell event:
  infection H→I (rate ∝ local virus), uninfected death (ROS + Allee `a_D`),
  infected death (ROS + contact-killing by NK/CD8⁺ + Allee), recovery D→H
  (cellularized Allee `a_H`, contact-area dependent).
- **Local immune cells** recruited from outside at ODE-driven inflow rates:
  macrophages (chemotax up virus, phagocytose, release chemokines/IL-10), NK and
  CD8⁺ T (chemotax up chemokines, haptotax, differential adhesion, contact-kill
  infected cells).
- **A global (non-spatial) ODE compartment** (Price 2015) for species that are
  well-mixed or act systemically: antibodies, TNF, IL-12, type-II IFN, CD4⁺ T,
  B cells, blood neutrophils, ROS, APCs. Coupled to the spatial model by
  cellularization scaling coefficients η (global) and θ (local).
- **Cellular viral resistance** `ρ = 1 − f̄/(θ·a_rf + f̄)` per cell from its local
  type-I IFN, modulating virus release and death. The paper notes ρ is the
  dominant computational cost.

Key parameters live in Tables 1–4 (diffusion lengths, `J` matrix, chemotaxis
λ's, volume constraint, Δt = 1 min, seeding fractions). Scenarios: *initial
viral load* (uniform virus, sweep 1–10000) and *initial infection fraction*
(random infected cells, sweep 0.001–0.05); small (0.3 mm) and large (1.0 mm)
patches; Neumann + periodic BCs.

## 3. Capability gap (engine today vs. paper)

| Capability | viva-cpm engine today | Gap |
|---|---|---|
| CPM: volume/surface, adhesion `J`, Boltzmann temp, parallel step | ✅ Rust (`crates/cpm-*`) | — |
| Diffusive fields: diffusion+decay, per-type secretion, per-type chemotaxis λ | ✅ `add_field`/`set_secretion`/`set_chemotaxis` | linear chemotaxis only |
| Log-chemotaxis with center-of-mass saturation (`λc·c/(1+c_CM)`) | ❌ | new energy term (Rust) |
| Haptotaxis (immune cells up adhesion gradients) | ❌ | new energy term (Rust) |
| Field sources/sinks regulated by cell state + other fields + globals | flat per-type rate | generalize source model (Rust + Python params) |
| Cellular viral resistance ρ (per-cell, from local IFN) | ❌ | per-cell field readout + ρ solve (Rust hot path) |
| Stochastic Table-2 type transitions (infection/death/recovery) | `set_cell_type`/`remove_cells` primitives only | transition Step (Python) + rate machinery |
| Contact-mediated killing (NK/CD8⁺ ↔ infected) | contact detection partial | contact-area query (Rust) + killing Step (Python) |
| Immune-cell recruitment/inflow at ODE rates | `add_cell` primitive only | recruitment Step (Python) |
| Global ODE compartment (Price 2015) + η/θ cellularization scaling | ❌ | ODE process (Python/scipy) + coupling |

## 4. Architecture (Approach A)

**Rust engine (`crates/cpm-core`, `crates/cpm-py`)** — the per-MCS hot path at
250k-cell scale:
- New energy terms: saturating/log chemotaxis with COM normalization;
  haptotaxis. Additive, behind per-type/per-field switches so existing linear
  chemotaxis is preserved.
- Generalized field source/sink model: secretion/uptake as a function of local
  concentration and a per-cell scalar (e.g. resistance ρ, cell state), not just a
  flat per-type constant.
- Per-cell field readout + ρ computation each MCS (`field_mean_at_cell` exists;
  extend to a batched ρ solve to contain the dominant cost).
- Contact-area-by-type query (for Allee `a_H` and contact killing).

**Python processes (`cpm/processes`, `pbg_cpm_studies/composites`)** — the
lower-frequency biology, orchestrated by `CPMProcess`-style Steps over the Rust
world:
- `TransitionProcess` — evaluates Table-2 stochastic transitions per cell per
  step (Poisson/exponential rates from local fields, contact counts, globals).
- `KillingProcess` — contact-mediated killing of infected cells by NK/CD8⁺.
- `RecruitmentProcess` — seeds new immune cells at ODE-driven inflow rates.
- `GlobalODEProcess` — integrates the Price-2015 global compartment (scipy),
  fed by spatial aggregates, feeding back transition/recruitment/field terms via
  η/θ scaling.
- `InfluenzaComposite` — wires the CPM world, the four fields, the transition/
  killing/recruitment steps, and the global ODE into one process-bigraph
  composite with the two IC scenarios as configuration.

**Rejected:** (B) pure-Python-on-current-engine — can't hit paper scale and
can't express the field-source/chemotaxis forms without per-MCS Python
callbacks; useful only as a low-scale cross-check. (C) ingest CC3DML directly —
couples us to CC3D serialization and still needs the same physics; instead the
CC3D source is transcribed once into a cited parameter authority.

## 5. Increment ladder

Each increment = one capability + one (or few) study + one paper target, and
lands as **its own PR** (capability code + tests + study + a visualization) so
reproducibility accrues cumulatively and every step is independently reviewable.
Increments 1–8 are validated at reduced scale (0.3 mm patches); Increment 9 is
the full 1 mm² capstone.

**Increment 0 — Spec authority & targets (no engine code).**
- Acquire the paper's supplementary CC3D source. Transcribe Tables 1–4 into a
  single cited `params.yaml` and extend `docs/cc3d-reference/` with a
  `sego2022-parameters.md` (formulas + constants + sources), following the
  existing cc3d-reference convention.
- Digitize the target curves from Figs 3B, 5, 7 into acceptance-band data.
- Scaffold the `influenza-sego2022` investigation shell.
- Study: `parameter-provenance` (documentation study, no sim).
- Provenance via `/viva-cite-bands`.

**Increment 1 — Epithelial sheet at scale.** 5×5-site cells, 2 µm lattice, vol
100 µm², adhesion `J`, Neumann+periodic BCs. Study `epithelial-sheet-baseline`:
cell-size/packing distribution + **MCS throughput at 250k cells** (perf budget is
a first-class acceptance item). Target: Tables 3–4 geometry.

**Increment 2 — Virus field + infection.** Virus field (D, decay → 5-cell-diam
length), infected-cell release gated by ρ (ρ stubbed = 1 until Incr. 3),
stochastic infection H→I ∝ local virus. Study `virus-field-infection`: infection
front speed, diffusion length. Target: Fig 6 virus; Table 3.

**Increment 3 — Type-I IFN field + resistance ρ.** IFN field (2-cell-diam
length), IFN production by infected cells, per-cell ρ modulating virus release &
death. Study `ifn-resistance`. Target: §2.2; Fig 6 IFN.

**Increment 4 — Epithelial death & Allee recovery.** Uninfected death
(ROS/Allee `a_D`), infected death, dead→uninfected recovery `a_H` (contact-area
Allee). Study `epithelial-fate`: death/recovery balance, critical contact area.
Target: Table 2; Fig 4A lesion recovery.

**Increment 5 — Chemokine + IL-10 fields.** Chemokine field (10-cell-diam,
macrophage-released, TNF-gated), IL-10 field (2-cell-diam). Study
`signaling-fields`: field profiles. Target: Fig 6 chemokines/IL-10.

**Increment 6 — Macrophage response.** Macrophage inflow, log-chemotaxis up
virus (λ=5000) with COM saturation, phagocytosis of virus/dead cells, chemokine/
IL-10 release. Study `macrophage-response`: localization to lesions. Target:
Fig 2B; Fig 3A immune layer.

**Increment 7 — NK & CD8⁺ T cells.** Haptotaxis (Rust), chemotaxis up chemokines
(λ=5000 / 10000), differential adhesion (homotypic 25, heterotypic 10,
preferential attach to infected), contact-mediated killing of infected cells.
Study `cytotoxic-killing`: killing rate vs contact, aggregation. Target: Fig 2;
Tables 2–3.

**Increment 8 — Global ODE compartment + cellularization scaling.** Full
Price-2015 global ODE (antibodies, TNF, IL-12, type-II IFN, CD4⁺, B, neutrophils,
ROS, APCs), η/θ scaling, bidirectional coupling. Study `global-coupling`:
ODE-vs-spatial consistency on global species. Target: Fig 3B global panels.

**Increment 9 — Capstone quantitative reproduction.** Full composite, paper
scale, 50 replicas:
- `repro-fig3b` — 5% initial infection, 0.3 mm patch, time series vs ODE.
- `repro-fig5-viral-load` — initial viral-load sweep (1/10/100/1000/10000),
  lethal threshold.
- `repro-fig7-infection-fraction` — infection-fraction sweep (0.001/0.005/0.01/0.05).
Pass iff observables land within the Increment-0 acceptance bands across the
replica ensemble. This is the "reproducibility achieved" deliverable.

## 6. Reproducibility acceptance

Two Increment-0 artifacts gate the capstone: (a) the cited `params.yaml`
transcribed from the CC3D source, and (b) the digitized Fig 3B/5/7 target curves
+ tolerance bands. Every increment's study asserts its own observable against its
paper target; the capstone asserts the full-model observables against the bands
across 50 replicas. Bands are set in Increment 0 (per-figure, per-observable);
where a per-version CC3D constant can't be retrieved verbatim, the *validated
behavior* (not the exact constant) is the fidelity criterion — the established
`docs/cc3d-reference/` convention.

## 7. Compute plan

The capstone (250k cells × ~20k 1-min steps × 50 replicas × multiple scenarios,
with the cost-dominant ρ solve) is heavy. Increments 1–8 are validated at 0.3 mm
scale on this laptop; the Increment-9 full-scale runs are budgeted to the **Mac
mini** (headless task agents, one worktree each). Increment 1 establishes the MCS
throughput budget that makes the capstone feasible; if the budget is missed, the
`step_parallel` path and the ρ batch solve are the optimization targets before
proceeding.

## 8. Sequencing & deliverables

- Branch: `investigation/influenza-sego2022` (this worktree). One PR per
  increment onto `main`.
- Increment 0's reference-PDF add is already committed (`c72821c`).
- After this spec is approved, `writing-plans` produces the detailed
  implementation plan for **Increment 0 and Increment 1** only; later increments
  get their own plan when reached (each is a sub-project: spec-section → plan →
  implement).
- TDD throughout (engine terms get Rust unit tests; processes get pytest;
  studies assert observables against targets).

## 9. Risks & open questions

- **Perf at paper scale.** The dominant risk. Mitigation: Increment-1 throughput
  gate; Rust hot path; mini for capstone. If 1 mm²/50-replica is infeasible in
  budget, fall back to fewer replicas or 0.5 mm patches and say so explicitly in
  the capstone study (no silent scale reduction).
- **Supplementary source availability.** If the exact CC3D source can't be
  retrieved, transcribe from the paper's tables/equations and mark constants
  CC3D-typical per the cc3d-reference convention.
- **Chemotaxis functional form.** The paper's saturating/COM form differs from
  the engine's linear ΔH_chem = −λ(c_dest − c_source); Increment 6/7 must match
  the paper form, not the current linear one.
- **Cellularization scaling (η, θ) calibration.** The ODE↔spatial mapping is the
  subtlest part; Increment 8 may need a calibration sub-study.

## 10. References & provenance

- `sego2022influenza` (this paper) — reproduction target.
- Price et al. 2015 — the underlying ODE model (to be added as a reference in
  Increment 0/8 when the global compartment is implemented).
- `glazier1993` — CPM foundation (already in refs).
- `docs/cc3d-reference/demo-parameters.md` — existing CC3D transcription
  convention this effort extends.
