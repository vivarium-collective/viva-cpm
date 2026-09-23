# Influenza immune-layer process + full_model Composite — design

**Status:** design (brainstorming output). Next: implementation plan via
`writing-plans` → `subagent-driven-development`.

**Motivating findings (this investigation):**
- `docs/cc3d-reference/reproduction-vs-cc3d-analysis.md` — the ~0.5-day lag / spatial comparison.
- `docs/cc3d-reference/chemokine-recruitment-scale-analysis.md` — after the
  Incr-13 cadence fix (PR #60), spatial immune counts are still ~2–4× below the
  fig3b bands; the gap is **source-faithful loop gain**, not a tunable constant.
  Both the chemokine-loop and ROS-ablation analyses reduce to one root cause: in
  the 2D CPM domain immune cells share the epithelial plane, so their density is
  space-capped and the local IL-10/chemokine environment that sets the
  recruitment loop gain is diluted; cytotoxic cells never reach the density to
  matter.

## Goal

Promote the influenza capstone science out of `run_full_model`'s monolithic
hand-written per-MCS Python loop into a real **process-bigraph `Composite`**, with
the immune layer as a **swappable, agent-based `ImmuneProcess`** whose immune
cells live *off the CPM lattice*. Decoupling immune density from epithelial
lattice space removes the root-cause space cap and lets the recruitment loop reach
its ODE-predicted equilibrium, while keeping the epithelium as a faithful CPM
process. The immune layer becomes swappable at the process-bigraph level (drop in
a different `ImmuneProcess` without touching the rest).

## Non-goals (explicit — no fudging)

- **No tuning of source constants** to hit figure bands. Recruitment/kill/
  secretion/ODE constants stay at their dossier values. The agent realization
  changes *geometry*, not rates.
- Not a native z=2 two-layer CPM (Rust Potts across layers) — rejected in favor
  of off-lattice agents (cheaper, uncapped, genuinely swappable).
- Not changing epithelial CPM dynamics, the diffusion engine, or the ODE math.
- Not deleting the current `run_full_model` path until the Composite reaches
  documented parity on Fig-3B.

## Success criteria

1. `full_model` Composite executes the full capstone science through
   process-bigraph (not a bypass), producing the same observable dict shape as
   `run_full_model`.
2. **Parity gate:** on Fig-3B at reduced scale, the Composite reproduces the
   current `run_full_model` epithelial + systemic trajectories within a documented
   tolerance (the only intended behavioral change is higher immune counts).
3. **Immune validation:** at paper scale (cells_per_side=35, 3.5 d), spatial
   immune counts rise materially above the current ~176-macrophage ceiling toward
   the ODE equilibrium (~380 chemokine-saturated), with no constant tuning.
4. Swapping `ImmuneProcess` for a stub/alternative requires editing only the
   Composite document's immune node.
5. Fast test suite stays green; new per-process unit tests added.

## Architecture

Four Processes wired through shared stores, executed by
`process_bigraph.Composite`. All four fields stay inside the Rust `World` (single
diffusion engine, single source of truth).

```
                 ┌───────────────────────────────────────────────┐
   fates ◄───────┤ EpitheliumProcess (extends CPMProcess)         │
   secretion ◄───┤  lattice + 4 fields + epithelial cell-fates    │
   field_deposit─►  (infection, resistance, secretion, apoptosis, │
                 │   Allee, ROS via ODE-X); Potts+diffusion/MCS    │
                 └───┬───────────────┬───────────────┬────────────┘
             types/positions/     field_at_point   field_means(all 4)
             volumes                 │                │
                 ▼                   ▼                ▼
   ┌───────────────────────┐   ┌───────────────────────────────┐
   │ ImmuneProcess ⭐       │   │ SystemicODEProcess            │
   │ off-lattice agents     │   │ GlobalODE: X, A, P, drivers   │
   │ {id,type,x,y}          │◄──┤ (per-record coupling)         │
   │ chemotax / kill / secre│──►│ reads H,I,M,K,E,field integrals│
   │ recruit spawn/remove   │   └───────────────────────────────┘
   └───────────────────────┘
```

**Stores:** `fates`, `secretion_scales`, `field_deposit`, `immune_agents`,
`ode_state`, `recruit_drivers`, plus the CPM readout stores (`types`, `positions`,
`volumes`, `field_at_cell`).

### 1. `EpitheliumProcess` (extends `CPMProcess`)

The spatial substrate. Owns the lattice, the four diffusive fields, and the
**epithelial** cell-fate transitions (moved in from the run driver so they run as
one ordered `update`, preserving intra-cell-scale source order): infection H→I,
per-cell resistance, secretion scales, infected apoptosis I→D, Allee death/
recovery, ROS death via ODE-X.

- **Ports out:** `types`, `positions` (COMs), `volumes`, `field_at_cell_all`
  (per-cell mean for *all four* fields, not just virus), `field_at_point`
  (sampler: nearest-site field value at continuous (x,y)).
- **Ports in:** `fates` (type overrides, incl. immune kills merged at a defined
  priority), `secretion_scales`, `field_deposit` (point sources: `[(field_idx,
  x, y, amount)]` injected into the field each MCS — how off-lattice agents
  secrete chemokine/IL-10).
- **Config:** the scene `spec`, `mcs_per_update`, field set.
- **Rust surface additions (small):** `field_value_at(idx, x, y)` (sample) and
  `field_add_source_at(idx, x, y, amount)` (point injection into a field site).
  `field_mean_at_cell(idx, cid)` already exists for all idx — the process just
  exposes every field, not only field 0 (today's `cpm_process.py:87` limit).

### 2. `ImmuneProcess` ⭐ (new — the swappable seam)

Owns off-lattice immune agents `{id, type ∈ {M,K,E}, x, y}` in a plain list store
(`immune_agents`). No lattice occupancy → **density uncapped**.

- **Per MCS (`update` at the CPM cadence):**
  1. **Chemotax:** sample chemokine (K/E) or virus (M) via `field_at_point` at
     each agent and a small neighborhood; move the agent a bounded step up the
     gradient (biased random walk; step/λ from the source chemotaxis constants,
     incl. the documented Increment-7 engine-scale note).
  2. **Kill:** for each cytotoxic agent (K/E), find infected epithelial cells
     whose COM (`positions`) is within `kill_radius`; contribute I→D to a
     dedicated immune-fate channel at the source `contact_kill_rate`/
     `nearby_kill_rate` magnitude (proximity replaces contact-area — a documented
     spatial approximation; constants unchanged).
  3. **Secrete:** macrophage agents emit `field_deposit` for chemokine + IL-10 at
     their (x,y), scaled by the same `macrophage_secretion_scale(sig_1, L)` used
     today (L sampled at the agent via `field_at_point`).
- **Per record (recruitment):** read `recruit_drivers` from the ODE (inflow/
  outflow rates) → spawn new agents (at tissue margin) / remove agents, per the
  source recruitment Poisson/outflow logic (unchanged from `recruitment.py`), but
  **without a reserve-pool cap** (agents are minted directly).
- **Ports in:** `field_at_point`, `positions`, `types`, `recruit_drivers`.
  **Ports out:** `fates` (immune-kill channel), `field_deposit`, `immune_agents`.
- **Config:** `kill_radius`, chemotaxis params, secretion rates, recruitment
  constants — all resolved from `params.yaml`/`resolve_constants`; none tuned.
- **Swap point:** replace this node's `address` in the Composite document.

### 3. `SystemicODEProcess` (new)

Wraps `price_ode.GlobalODE`. Per record: read spatial aggregates (H, I, M, K, E =
agent counts by type, and the four field integrals) → integrate the 10 systemic
species → write `ode_state` (incl. X for ROS, A for antibody clearance, P for APC)
and `recruit_drivers` (the six inflow/outflow rates). Cadence: once per
`mcs_per_step` MCS (the corrected recruitment-every-MCS behavior is preserved
because `ImmuneProcess` applies the drivers every MCS between ODE updates).

### 4. Composite document

`full_model_composite_document` (`viva_cpm_studies/composites/influenza.py:452`) is
rewritten from the current CPM-only scene to wire the four processes + stores
above. `run_full_model` becomes a thin driver that builds and steps a
`process_bigraph.Composite` from this document and reads the observable series
from the stores — replacing the hand loop.

## Execution & ordering

- **Cadence:** `EpitheliumProcess` and `ImmuneProcess` step every MCS;
  `SystemicODEProcess` couples once per record (`mcs_per_step` MCS).
- **Intra-MCS ordering approximation (accepted, documented):** cross-process
  reads within one MCS use the MCS-start grid/fields; the strict sequential reads
  of the source pipeline (infection→…→killing→apoptosis) are preserved *inside*
  `EpitheliumProcess`, while `ImmuneProcess` kills via a separate fate channel that
  `EpitheliumProcess` merges at a fixed priority. This is second-order for the
  model's low-probability per-MCS transitions and the model is stochastic; a
  **parity test bounds the divergence** from the current sequential pipeline.

## Data flow (one record)

1. For each of `mcs_per_step` MCS: `EpitheliumProcess.update` (Potts+diffusion +
   epithelial fates, consuming `fates`/`secretion_scales`/`field_deposit` from the
   previous MCS) ‖ `ImmuneProcess.update` (chemotax/kill/secrete, producing
   `fates`/`field_deposit`/`immune_agents`).
2. `SystemicODEProcess.update` reads aggregates + field integrals, integrates,
   writes `ode_state` + `recruit_drivers`.
3. `ImmuneProcess` applies `recruit_drivers` (spawn/remove) over the next record's
   MCS; `EpitheliumProcess` reads `ode_state["X"]` for ROS death.

## Migration / phasing

- **Phase 1 — EpitheliumProcess.** Extend `CPMProcess`: all-field per-cell
  readout, `field_at_point`, `field_add_source_at` (Rust), point-deposit input,
  and move the epithelial fate transitions in. Parity test: an EpitheliumProcess-
  only run (immune off) matches `run_full_model(enable=epithelial-only)`.
- **Phase 2 — ImmuneProcess.** Agent chemotaxis/kill/secretion/recruitment. Unit
  tests per behavior; integration test: immune counts exceed the reserve-pool
  ceiling and rise toward the ODE equilibrium.
- **Phase 3 — SystemicODEProcess + Composite driver.** Wire the Composite; make
  `run_full_model` drive it (behind a flag first). End-to-end Fig-3B parity on
  epithelial+systemic series; then the paper-scale immune-count check.
- The plain-function path (`killing.py`, `recruitment.py`, the hand loop) stays
  until the Composite passes the parity gate; then the driver switches over.

## Testing

- Per-process unit tests (ports, one-step behavior) using the `InfectionProcess`
  test as a template.
- **Parity tests** (Phase 1 & 3) vs the current `run_full_model` within a stated
  tolerance — the guardrail against silent behavior change.
- Immune-count integration test (Phase 2/3) at reduced + paper scale.
- Fast suite stays green; heavy paper-scale runs on the mini.

## Global constraints

- **No constant tuning to fit bands.** Only spatial *realization* changes
  (off-lattice agents, proximity killing) — documented as a cellularization
  approximation alongside the existing 2D approximation, with the source dossier
  (`docs/cc3d-reference/sego2022-*`) as authority.
- The agent-immune layer is a **documented departure** from CC3D's z=2 CPM immune
  cells; record it as such (it is the deliberate design choice, not an accident).
- Worktree discipline (dedicated worktree off current `origin/main`); no AI
  attribution in commits/PRs; manual merges.

## Open risks

- **Agent-kill vs contact-kill calibration:** proximity radius replaces contact
  area. Choose `kill_radius` from the CPM cell geometry (a cell's contact
  footprint), not fitted to outcomes; document the mapping. Verify the per-agent
  kill probability matches the source `contact_kill_rate` in the relevant limit.
- **Ordering divergence:** bounded by the Phase-1/3 parity tests; if a transition
  proves order-sensitive, fall back to sub-step scheduling for that pair.
- **Performance:** many agents × per-MCS field samples; keep the field sampler
  O(1) (nearest-site) and vectorize agent updates.
- **Composite/driver bridge:** first end-to-end where `run_full_model` drives a
  `process_bigraph.Composite`; Phase 1 de-risks the Rust/port surface early.
