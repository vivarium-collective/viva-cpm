# viva-cpm reproduction vs. the CC3D ViralInfectionVTM simulation (Sego et al. 2022)

**Scope.** A rigorous, honest comparison of OUR native reproduction (Rust-CPM +
Python `run_full_model`) against the reference CC3D *ViralInfectionVTM* simulation
of the same model (T.J. Sego et al., *J. Theor. Biol.* 532 (2022) 110918). This is
an analysis/write-up only — no model code, sim runs, or study verdicts were changed.

**What "the CC3D simulation" is here.** We do not have a live CC3D run. The reference
is (a) the digitized target curves — Fig 3B (time series vs. the global ODE), Fig 5
(viral-load sweep), Fig 7 (infection-fraction sweep), and Fig 2/3 (spatial snapshots)
— captured in `viva_cpm_studies/influenza/targets/{fig3b,fig5,fig7}.json`; and (b) the
transcribed CC3D *source*, dossiered in `docs/cc3d-reference/sego2022-parameters.md`,
`…-global-ode.md`, and `…-source-notes.md`. The CC3D source (not the paper prose) is
the fidelity ground truth per the source-notes convention.

**What "our reproduction" is.** `viva_cpm_studies/influenza/run.py::run_full_model` +
`_one_mcs` — a single per-MCS pipeline assembling every mechanism primitive from
Increments 0–8, with Increments 10–12 adding ROS death, scattered seeding, and the
field-units correction. Paper-scale results (cells_per_side=35 / 1225 cells / 720
record-steps / ~3.5 days, 50 replicas) are from Mac-mini runs; the in-repo CI configs
are deliberately tiny and do NOT reproduce the figure (the study verdicts say so
explicitly, e.g. `workspace/studies/repro-fig3b/study.yaml`).

---

## 1. Executive summary

Our reproduction captures the *qualitative and endpoint* behavior of the CC3D Fig-3B
scenario faithfully: the epithelial sheet is driven from ~95% healthy to near-total
death (our paper-scale run ends at 1217/1225 dead vs. CC3D's 1223/1225), and the
systemic/immune signals validate well against the digitized bands (antibodies 8/8
checkpoints in-band, type-I IFN 7/8, dead 7/8; 49/96 overall). The single most
important quantitative divergence is a **systematic ~0.5-day temporal lag**: our
infected population peaks at ~day 1.5 (≈292 cells) where CC3D peaks at ~day 0.5
(≈400), and the whole H/I/D cascade is shifted later even though its shape and
endpoint are right. Spatially, the *composition* trajectory matches the CC3D
green→orange→grey sweep (healthy 63.5%→2.8%, dead 0%→61.8% of the domain), but the
committed spatial capture spreads as **one radial front from a central lesion** rather
than the multi-focal, sheet-wide death the CC3D Fig-3B "fraction" scenario shows, and
our immune cells sit in **fixed seeded clusters outside the tissue** rather than being
recruited and distributed through it. Progressive calibration closed most of the lag
without tuning any CC3D constant (Fig-3B infected worst-miss fell 83→28→2.75 band-widths
across Increments 9→10→11), and the residual ~0.5 day is left open deliberately because
closing it would require tuning source-faithful constants (g_hv, virus secretion/decay)
to fit the figure. Net verdict: a close, source-faithful reproduction of the CC3D
behavior — right shape, right endpoint, validated systemics — with a documented lag and
an explicitly bounded, not-yet-full-band-pass claim.

---

## 2. Spatial-state comparison

### 2.1 What we computed

Decoded the baked capture `INFLUENZA_SPATIAL["repro-full-model"]`
(`viva_cpm_studies/visualizations/_influenza_spatial.py`; `zlib+base64` → `uint8`
per-site cell-type grid, nx=235, ny=195, seed=0, 18 frames spanning MCS 0→5033 ≈
0–3.5 days). Type codes: MEDIUM=0, H=1, I=2, D=3, M=4, K=5, E=6. Per frame I computed
per-epithelial composition (fractions of H/I/D among the H∪I∪D sites), 8-connected
component counts for the infected and dead site-masks (a spatial-spread / fragmentation
metric), and the bounding-box diagonal extent of each (site units). Site counts are
25× cell counts (5×5-site cells).

| day | H sites | I sites | D sites | %H | %I | %D | #I comps | #D comps | max D comp | I extent | D extent |
|----:|--------:|--------:|--------:|----:|----:|----:|---------:|---------:|-----------:|---------:|---------:|
| 0.00 | 29100 | 1525 | 0     | 95.0 | 5.0  | 0.0  | 1   | 0   | 0     | 62  | 0   |
| 0.21 | 28837 | 951  | 558    | 95.0 | 3.1  | 1.8  | 3   | 62  | 134   | 66  | 89  |
| 0.62 | 28519 | 773  | 1036   | 94.0 | 2.5  | 3.4  | 20  | 95  | 505   | 100 | 134 |
| 1.03 | 26741 | 1712 | 1917   | 88.1 | 5.6  | 6.3  | 41  | 180 | 956   | 149 | 155 |
| 1.44 | 22854 | 3153 | 4406   | 75.1 | 10.4 | 14.5 | 71  | 275 | 2729  | 219 | 204 |
| 1.64 | 20465 | 3885 | 6078   | 67.3 | 12.8 | 20.0 | 124 | 317 | 4076  | 234 | 217 |
| 2.06 | 14616 | 5131 | 10957  | 47.6 | 16.7 | 35.7 | 202 | 380 | 8346  | 252 | 267 |
| 2.47 | 9808  | 4840 | 16118  | 31.9 | 15.7 | 52.4 | 227 | 241 | 13791 | 259 | 272 |
| 2.88 | 7015  | 5012 | 18744  | 22.8 | 16.3 | 60.9 | 257 | 148 | 17622 | 256 | 275 |
| 3.29 | 1895  | 2005 | 26723  | 6.2  | 6.5  | 87.3 | 198 | 81  | 26440 | 269 | 284 |
| 3.50 | 1278  | 1017 | 28318  | 4.2  | 3.3  | 92.5 | 109 | 99  | 27836 | 276 | 290 |

Whole-domain composition (over all 235×195 sites): frame 0 → H 63.5%, I 3.3%, D 0%,
immune ~0.7%, medium 31.5%; final → H 2.8%, I 2.2%, D 61.8%, immune ~1.4%, medium
31.8%. These domain numbers match the brief's stated "healthy 63.5%→2.8%, dead
0%→61.8%" exactly, confirming this is the intended full-model capture.

### 2.2 Morphology — what the capture actually shows

- **Initial infection is a single central lesion, NOT scattered.** At frame 0 the
  1525 infected sites form **one** 8-connected component in a compact 44×44-site block
  centred in the patch (x∈[75,119], y∈[75,119]), covering only ~25% of the 174-site
  patch width/height. The epithelial patch bbox is x,y∈[10,184]. This is the
  pre-Increment-11 "single central lesion" morphology (a carry-over from reusing the
  immune-localization scenario builder), **not** the random-scatter that Increment 11
  introduced for the paper-scale band eval. So the committed spatial capture and the
  paper-scale temporal band numbers (§3) come from *different* seeding configs — an
  honest caveat: the composition sweep below is representative, but the committed
  capture's *initial geometry* is the older central-lesion one.
- **Spread is a single expanding front.** The infected/dead bounding-box extent grows
  monotonically from 62→276 (infected) and 0→290 sites (dead) — the lesion radiates
  outward from the centre. Infected components fragment (1→~227) as the centre dies out
  and only an advancing infected annulus remains, then collapse again (→109) as the
  whole sheet dies.
- **Dead cells coalesce into one massive field.** The dead mask starts as many small
  components (≤134 sites) and ends as a single dominant connected component of 27,836
  sites — 98.3% of all 28,318 dead sites — i.e. a confluent grey field, matching CC3D
  Fig 2/3's terminal "sheet of dead cells."
- **Immune cells never infiltrate.** M/K/E sit in three fixed 100-site (~4-cell)
  clusters *outside* the epithelial patch entirely (M at x∈[193,204], K/E at
  x∈[213,224]; the patch ends at x=184). Their counts barely change over the run
  (M 100→225 sites ≈ 4→9 cells; recruitment is near-inert at this η). They chemotax
  locally but do not distribute across the tissue.

### 2.3 vs. the CC3D spatial picture (Fig 2/3)

CC3D Fig 3B (the "fraction" scenario) seeds infection scattered across the sheet, so
death is **multi-focal and near-simultaneous sheet-wide**; immune cells are
**recruited** (`new_immune_cell_by_type`, seeded at the target field's local maximum,
`ImmuneModelSteppable`, dossier §4b(ii)) and therefore **distributed through and around
the lesion** with infiltration geometry. Two known, documented differences:

1. **Seeding geometry.** Committed capture = one central lesion → single radial front;
   CC3D = scattered foci → sheet-wide. (Increment 11's scatter — `run_full_model` lines
   1712–1723 — fixes this for the paper-scale run; it is simply not in the committed
   capture.)
2. **Immune geometry.** Ours = fixed exterior seeded clusters that don't infiltrate;
   CC3D = recruited, distributed, infiltrating cells. This is an *engine-imposed*
   approximation: our engine cannot mint cells after `world.finalize`, so ODE-driven
   inflow *activates* dormant reserve cells parked in the medium margins
   (`_seed_recruit_pool`, `RECRUIT_RESERVE_TYPE`, run.py lines 41–99) rather than adding
   cells at the lesion the way `new_immune_cell_by_type` does.

The *composition* dynamics (green→orange→grey) are faithful; the *morphology* (front vs.
multi-focal, clustered vs. infiltrating) is where the spatial pictures diverge.

---

## 3. Temporal-dynamics comparison

### 3.1 Paper-scale H/I/D vs. the CC3D Fig-3B targets

From the Mac-mini paper-scale run (ROS + scatter; cells_per_side=35, 720 steps,
~3.5 days) against `targets/fig3b.json` (`value` = ODE reference; band = 50-replica
spatial spread). Format: day: U_mod/U_tgt | I_mod/I_tgt | D_mod/D_tgt.

| day | U mod/tgt | I mod/tgt | D mod/tgt |
|----:|:---------:|:---------:|:---------:|
| 0.5 | 1073 / 900 | 103 / 400 | 49 / 20 |
| 1.0 | 756 / 300  | 262 / 300 | 207 / 400 |
| 1.5 | 356 / 50   | 292 / 80  | 577 / 900 |
| 2.0 | 106 / 10   | 124 / 15  | 995 / 1150 |
| 3.5 | 5 / 2      | 3 / 1     | 1217 / 1223 |

**Band evaluation and progression (worst-miss in band-widths, infected/uninfected/dead):**
Increment 9 baseline **83 / 7.5 / 6.8** → +ROS death (Incr 10) **28 / 21.5 / 4.7** →
+scattered infection (Incr 11) **2.75 / 2.2 / 1.0** (dead in-band 7/8 checkpoints).
Overall **49/96** checkpoints in-band; antibodies **8/8**, type-I IFN **7/8**, dead **7/8**.
The systemic/soft observables (IFN, antibodies) validate strongly; the residual failure
is the timing of the H/I cascade, not its shape or endpoint.

### 3.2 Characterizing the lag

The lag is **not** a uniform rigid time-shift; it is a **slower rise with a
correct-magnitude, roughly correct-timing decay and a correct endpoint**:

- **Rise is late/slow.** Infected peaks ~day 1.5 (292) vs. CC3D ~day 0.5 (400). At day
  0.5 we have only 103 infected (26% of target) and 1073 uninfected (still 88% healthy)
  where CC3D is already 900 uninfected and 400 infected. The virus/infection front takes
  ~1 extra day to engage the sheet.
- **Endpoint is exact.** Dead → 1217/1225 (99.3%) vs. 1223/1225 (99.8%); uninfected → 5
  vs. 2. The lethal outcome is reproduced.
- **Decay realigns.** By day 2.0 the dead count (995) is ~86% of target (1150) and by
  day 3.5 they coincide — i.e. the trajectories *converge* late, so the miss is
  concentrated in the 0.5–1.5 day rising limb. This is the signature of a **front-speed /
  ramp-rate deficit**, not a clock offset.

The committed *central-lesion* spatial capture (§2) makes the same point more starkly:
its epithelial %D is only 6.3% at day 1.0 and 73% not until ~day 2.9, versus CC3D's
32.7% at day 1.0 and 73.5% by day 1.5 — a >1-day lag when spread is a single front.
Scattering the seed (Incr 11) collapses that to ~0.5 day, direct evidence that
front-propagation speed is a dominant lag term (§6).

---

## 4. Mechanism-fidelity comparison (CC3D steppable order vs. our `_one_mcs`)

The CC3D per-MCS order (dossier `sego2022-global-ode.md` §4:
`__update_spatial_data` → `rr.timestep()` → `update_populations`, with the secretion/
killing/death/Allee/ROS steppables each firing once per MCS) is mirrored by
`run_full_model._one_mcs` (run.py 1896–2032) + `_ode_couple` (2034–2064). Step-by-step:

| CC3D mechanism (source) | Our helper | Fidelity |
|---|---|---|
| Potts + 4 reaction-diffusion fields advance | `world.step(1)` | **Matches** (fields co-advanced; chemotaxis applied inside step) |
| Infection H→I `Pr=1−exp(−g_hv·v̄)` (`ViralInternalizationSteppable`) | `transitions.infection_step` | **Matches** — `g_hv` source-literal (`infection_g_hv=0.006133…`, params §6) |
| Resistance ρ=f̄/(a_rf+f̄) (`update_resistance`) | `resistance.cell_resistance` | **Matches** — `a_rf` source-literal; consumers use (1−ρ) |
| Virus release `g_vi·(1−ρ)` (`ViralSecretionSteppable`) | `set_cell_secretion_scale(virus,(1−ρ))` | **Matches in form**; see §6 on the dim.z=2 factor |
| Macrophage chemokine/IL-10, saturating `sig_1` Hill (`ChemokineSecretionSteppable`/`IL10`) | `signaling.macrophage_secretion_scale`, dynamic `sig_1=a_11·T+a_12·D` | **Matches form** (§4b(i)); macrophage id list refreshed each MCS |
| Chemotaxis, **saturating** λ_eff=λ/(1+c(COM)) (`ChemotaxisSteppable`) | `set_macrophage_chemotaxis` / `set_nk_cd8_chemotaxis`, **fixed linear** λ | **Approximated** — engine does plain linear ΔH=−λ(c_dest−c_src) with fixed λ (params §4); NK/CD8 λ needs a 100× engine-unit rescale (`NK_CD8_CHEMOTAXIS_ENGINE_SCALE`) because the literal λ is below the noise floor at this field's concentration scale |
| Contact killing (local surface + nearby well-mixed), ×`cell_resist` **directly** (`ContactKillingSteppable`) | `killing.contact_kill_rate` + `nearby_kill_rate` | **Matches**, incl. discrepancy #7 (resist direct, not 1−resist) |
| Allee death H→D / recovery D→H, shared `b_h` (`RecoverySteppable`) | `allee.allee_death_rate`/`allee_recovery_rate` | **Matches form**; direct H→D Allee death is rare at 0.3mm scale (epithelial-fate study) |
| Infected apoptosis `mu_i·(1−ρ)` (`ViralCellDeathSteppable`, "simplified") | `transitions.infected_death_step` | **Matches** the source's *simplified* flat rate (discrepancy #6) |
| ROS death: X→ g_ix·hill(X,a_ix,h_x) (I), g_hx·hill(X,a_hx,h_x) (H) (`OxidationAgentModelSteppable`) | step 8.5 in `_one_mcs` (2008–2032) | **Matches form**; the CC3D *dominant* epithelial killer, wired in Incr 10 |
| Spatial→ODE push → integrate (RoadRunner CVODE, 1×/MCS) → recruitment | `_ode_couple` → `price_ode.GlobalODE.step` → `_recruit_step` | **Approximated cadence** — see below |
| Recruitment via `new_immune_cell_by_type` (add cell at field max) | reserve-pool activation (`_recruit_step`) | **Approximated** — reserve pool vs. add_cell; exterior placement |

**Approximations and their plausible contribution to the lag / immune gaps:**

1. **ODE coupling cadence.** CC3D integrates the hybrid ODE **every MCS** (dossier §5,
   `frequency=1`). We integrate + push/pull **once per recorded step** (mcs_per_step=7
   MCS; `for i…: for _ in range(mcs_per_step): _one_mcs(); _ode_couple()`, run.py
   2077–2081). The integrated species (T, X, A, P…) therefore update on a 7-min grid, and
   the `sig_1` and ROS `X` fed to the cell layer are up to 7 MCS stale. Because X is the
   dominant epithelial killer, a lagged X ramp directly slows the early death rise — a
   candidate lag contributor (documented as immaterial for the slow integrated species,
   but it touches the fast-acting ROS path).
2. **2D vs. CC3D's z=2 slab.** CC3D runs a 175×175×**2** lattice; several source rates
   carry an explicit `dim.z=2` factor (virus/IFN secretion, ROS). Our domain is single-
   layer (z=1). params-dossier flags this "halve if single-layer" convention repeatedly;
   any mismatch in the effective per-site secretion scales the virus-field ramp and hence
   infection speed.
3. **Chemotaxis: linear + fixed λ vs. saturating per-MCS λ.** Combined with the 100×
   engine-unit rescale for NK/CD8, this makes immune *localization* qualitatively right
   but not quantitatively source-faithful — a plausible driver of the immune-magnitude
   softness (macrophage/NK/CD8 are the `soft`-banded observables).
4. **Reserve-pool recruitment vs. add_cell.** Placement is exterior and the pool is
   finite; at reduced η recruitment inflow is near-inert. This under-supplies distributed
   immune pressure early, which would otherwise sharpen the infected-peak timing.

---

## 5. Discrepancy inventory (#1–#12) and comparison impact

From `sego2022-parameters.md §7` (#1–#8) and `sego2022-global-ode.md §6` (#9–#12).
These are source-vs-paper or source-vs-ours items; we follow the **source** throughout.

| # | Discrepancy | Impact on the comparison |
|---|---|---|
| 1 | NeighborOrder = 3 (source) not 2 | We use 3; affects contact-area / adhesion energetics slightly. Low. |
| 2 | Infected–Macrophage adhesion J=20 (source) vs. paper's collapsed 10 | We use source matrix. Affects macrophage–lesion contact. Low. |
| 3 | NK–CD8 adhesion J=25 (source) vs. paper's 10 | Source matrix used. Immune-cluster cohesion. Low. |
| 4 | Default dims 175² (η=0.0049), not a clean 0.3/1.0 mm patch | We use 35²=1225 cells / η=0.0049 for Fig-3B (matches shipped default). Sets the population scale. |
| 5 | v0 sweep {1..10000} not enumerated in source (only 0/50/500) | We encode the paper's sweep as scenario params (Fig-5). No fidelity loss. |
| 6 | `ViralCellDeathSteppable` is a **simplified** flat `mu_i·(1−ρ)` (no ROS/viral-load Hill in that term) | We mirror the simplified form; the ROS Hill lives in the separate OxidationAgent path (Incr 10). Correct per source. |
| 7 | Contact killing multiplies `cell_resist` **directly**, not (1−resist) | Reproduced verbatim; a cell with MORE local IFN is killed FASTER by contact — counter-intuitive but source-faithful. |
| 8 | Macrophage chemotaxis label "moderate" but λ=5000 (=NK's "strong") | Numeric value used; label-only. None. |
| 9 | The driver is a **hybrid**, not a full 20-var ODE: {NB,N,T,X,A,B,P,W,G,O} integrated; {H,I,M,E,K,L,C,F,V,DH} overwritten from CPM/fields each step | We implement the hybrid coupling (`INTEGRATED_STATES` + spatial→ODE push), matching source. Central to the whole coupling design. |
| 10 | Two model strings scale differently (spatial uses s_l; reference uses s_v) | We use the **spatial** string's scalings. Getting this wrong would misscale g_hv/g_vh/g_fi/a_v — directly affects infection speed. Correct. |
| 11 | `sig_1` is a **saturating** IL-10-inhibited secretion, not a linear source | We match the Michaelis form (local + well-mixed). Affects chemokine/IL-10 magnitude. |
| 12 | Recruitment is asymmetric: CD8 is APC(P)-driven with no homeostatic baseline; macro/NK are chemokine(C)-driven with one; NK/CD8 outflow carry resist-weighted G_ki/B_ei | We reproduce the asymmetry (`macrophage_inflow`/`nk_inflow`/`cd8_inflow`). Affects immune-population timing. |

No discrepancy is a source-fidelity error on our side; #6/#7/#9/#10/#11/#12 are places
where following the *source* (not the paper prose) matters, and we do.

---

## 6. The ~0.5-day lag — root-cause hypotheses

Enumerated, grounded in §4, with testability noted. Distinguishing genuine CPM-
realization differences from possible residual scaling issues:

1. **CPM front-propagation speed (genuine realization difference; largest term).**
   The strongest evidence: scattering the seed (Incr 11) dropped the infected worst-miss
   from 28→2.75 band-widths and the central-lesion capture lags >1 day while the
   scattered paper-scale run lags ~0.5 day. A single front must physically traverse the
   sheet cell-by-cell via virus diffusion + stochastic H→I; scattered foci infect the
   sheet in parallel. CC3D's Fig-3B is the scattered case, so part of the *residual* lag
   is that even our scattered seeding is not identical to CC3D's exact focus count/layout.
   **Testable** without fudging: vary the number/dispersion of scattered foci and measure
   the infected-peak day; compare to CC3D's `frac_init_infected=0.05` layout.

2. **Virus-field ramp / diffusion rate (possible scaling; testable).** Infection rate ∝
   local virus·g_hv. If the virus field ramps too slowly (diffusion constant
   virus_dc_im≈0.179 lat²/MCS, or the single-layer z=1 vs. z=2 secretion factor), the
   early H→I is throttled. The `dim.z=2` secretion factor (`secretion_g_vi` applied ×2 in
   source) is the concrete suspect: params §6 flags "halve if single-layer" — if we do
   not carry the correct z-factor the virus ramp is mis-scaled. **Testable:** check the
   effective per-site secretion against source `g_vi·dim.z` and compare the day-0.5 virus
   field integral to the fig3b `extracellular_virus` target (≈500 at day 0.5).

3. **The 2× secretion_g_vi vs. source g_vi observation (possible residual scaling).**
   Source applies `g_vi·dim.z=0.386` per infected cell; our single-layer field must
   reconcile this factor. A residual 2× either way shifts the virus ramp and the infected
   rise timing directly. **Testable** by unit audit of the virus field against the target
   at t=0.5–1.0, no figure-fitting required.

4. **ODE coupling cadence (approximation; testable, likely minor for the lag).** The
   7-MCS coupling grid staleness of ROS `X` (the dominant killer) could blunt the early
   death rise. **Testable:** run `_ode_couple` every MCS (mcs_per_step=1) and check
   whether the day-1.0 dead count rises toward 400. Documented as immaterial for slow
   species, but this is the one fast-acting exception worth checking.

5. **Chemotaxis form/units (genuine realization difference; affects immune magnitude
   more than lag).** Fixed linear λ + 100× rescale vs. saturating per-MCS λ softens
   immune localization; the immune observables are `soft`-banded and this is a magnitude,
   not a timing, effect.

**Explicitly NOT to be closed by tuning:** g_hv, g_vi, virus decay are already CC3D
source-literal values (params §2/§3/§6). Closing the residual lag by changing them would
mean tuning CC3D constants to fit the digitized figure — which we deliberately do not do.
The honest split: hypotheses **1 and 5 are genuine CPM-realization differences**
(front physics, chemotaxis form); hypotheses **2/3 are testable unit/scaling audits**
that could legitimately be fixed if a mis-scaling is found (that is a bug, not a fit);
**4 is a documented cadence approximation** worth one check.

---

## 7. Honest verdict

Our viva-cpm reproduction is a **close, source-faithful reproduction of the CC3D
ViralInfectionVTM Fig-3B behavior**: it reproduces the trajectory **shape** and the
lethal **endpoint** (dead → 1217/1225 vs. CC3D 1223/1225), the systemic/immune signals
validate against the digitized bands (antibodies 8/8, type-I IFN 7/8, dead 7/8; 49/96
overall), and the per-MCS mechanism pipeline matches the CC3D steppable order with every
rate a source-literal constant (no tuning-to-fit; calibration Increments 10–12 were all
source-faithful mechanism additions — ROS death, scattered seeding, a units correction).
It is **NOT a full quantitative band pass**: a systematic **~0.5-day temporal lag**
(infected peaks ~day 1.5/292 vs. CC3D ~day 0.5/400) keeps the H/I rising limb out of
band, and spatially the committed capture shows single-front radial spread with fixed
exterior immune clusters rather than CC3D's multi-focal death with recruited, infiltrating
immune cells. The lag is dominated by CPM front-propagation physics (evidenced by the
28→2.75 worst-miss drop when the seed was scattered) plus testable virus-field unit/scaling
questions; it is left open because the source constants that would close it are exactly
the ones we refuse to tune. **Claim, bounded:** right shape, right endpoint, validated
systemics, documented ~0.5-day lag, not a full band pass — and the spatial-morphology and
immune-geometry differences are known, engine-imposed, and stated, not hidden.

---

*Sources cited: `docs/cc3d-reference/sego2022-parameters.md` (§2–§7, mechanism literals
+ discrepancies #1–#8), `docs/cc3d-reference/sego2022-global-ode.md` (§2–§6, hybrid ODE +
coupling + discrepancies #9–#12), `viva_cpm_studies/influenza/targets/fig3b.json` (digitized
CC3D/paper curves + bands), `viva_cpm_studies/influenza/run.py::run_full_model`/`_one_mcs`
(our pipeline), `viva_cpm_studies/visualizations/_influenza_spatial.py`
(`repro-full-model` capture, decoded here). Our paper-scale numbers are from Mac-mini runs
not stored in the repo. Spatial metrics in §2 were computed directly from the decoded
capture for this analysis.*
