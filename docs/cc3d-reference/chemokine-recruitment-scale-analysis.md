# Where the immune response diverges: the chemokine→recruitment loop gain

**Context.** After the Increment-13 fixes (reserve-pool uncapping + recruitment
applied every MCS, PR #60), paper-scale fig3b macrophages rose ~20× (9 → 176 at
3.5 d) but still land ~2× below the fig3b lower band (M 380, target 700). This
note traces *why*, and concludes the remaining gap is **not** a mis-scaled
constant — so it must not be closed by tuning.

## The chain, measured (cells_per_side=35, bootstrap, seed 0)

Recruitment inflow is `b_mc·hill(C, a_mc, h_m) + μ_m·b_m` — a chemokine(`C`)-driven
Hill on a homeostatic baseline. At fig3b scale `a_mc = 2.82`, so recruitment
half-saturates only when `C ≳ 2.82`.

Measured realized chemokine field `C` (= `chemo_field_integral / dim_z`, the exact
quantity fed to recruitment):

| t (d) | M | C | C/a_mc | hill(C) | inflow | M_eq = inflow/μ_m |
|------|----|------|-------|---------|--------|------|
| 0.2 | 16 | 0.22 | 0.08 | 0.001 | 0.0176 | 105 |
| 0.5 | 23 | 0.34 | 0.12 | 0.002 | 0.0177 | 105 |
| 0.8 | 29 | 0.73 | 0.26 | 0.017 | 0.0184 | 109 |

`C` never approaches `a_mc`, so `hill(C) ≈ 0` — **recruitment runs on the
homeostatic baseline only** (M_eq ≈ 105), which is exactly where the macrophage
count settles. The chemokine term is inert because `C` is ~5–10× too small.

## Why C is small: the macrophage secretion scale

Chemokine is secreted by macrophages at `b_c·M·scale`, where (source
`ImmuneModelLib`, dossier §4)

```
scale = sig_1 / (sig_1 + (g_1·L + g_2)/(L + d_2)),   sig_1 = a_11·T + a_12·D
```

Instrumenting the actual per-macrophage calls (14 644 samples over a 0.58-day
run):

- per-macrophage local IL-10 `L ≈ 0` (median 0.0000, max 1e-4) → the ratio term
  sits at its floor `g_2/d_2 = 0.062`
- `sig_1 ≈ 0.0006` (median), max ~0.024 even at total epithelial death
  (`a_12·D`, D≤1225, a_12 = 1.93e-5; the `a_11·T` term is negligible, T~0.03)
- ⇒ **applied secretion scale ≈ 0.0025** (median), i.e. macrophage chemokine
  output is throttled to ~0.3% because `sig_1 ≪ g_2/d_2`.

So the loop is: little dead tissue → tiny `sig_1` → throttled chemokine → `C ≪
a_mc` → baseline-only recruitment → few macrophages → little chemokine. It only
opens up as `D` grows, and by then (in our ROS-dominated model) the infection has
already cleared.

## This is source-faithful scaling, not a bug

The tempting "fix" is to rescale `a_11`/`a_12` (recorded `scale: null`) or the
IL-10 constants. **That would be fudging.** The cellularization is self-consistent:

- `sig_1 = a_11·T + a_12·D`. The ODE integrates `T, D` in cellularized units
  (`resolve_constants(num_epithelial = tot_cell = 1225)`), so `T_cell ≈ s_v·T_src`
  and `D_cell = D_src·s_v` (1225 = 250000·s_v). Hence `sig_1_cell ≈ s_v·sig_1_src`.
- The ratio floor `g_2/d_2 = il10.g_2·s_v / il10.d_2` also carries one `s_v`
  (source floor 12.6 → 0.062, ×s_v).

Both numerator and denominator scale by the same `s_v`, so
`scale = sig_1/(sig_1+ratio)` is **scale-invariant** — the source sees the same
throttled secretion (source-equiv sig_1 ≈ 0.12 vs ratio floor 12.6 → scale ≈
0.01). `a_11`/`a_12` are correctly `scale: null` (dossier lines 152–153), and
`g_1`/`a_mc`/`g_2` carry their documented `s_v` factors. There is no isolated
mis-scaled constant.

## Conclusion — where viva-cpm actually diverges

The ~2× immune shortfall is **emergent loop gain + timing**, not a broken constant:

1. **Geometry (2D vs CC3D z=2).** Local IL-10 `L ≈ 0` at macrophages and the
   chemokine field are 2D-diluted. CC3D's immune cells live in a z=2 layer with a
   different local-field environment; the same secretion scale acting in that
   geometry, with macrophages packed off the epithelial plane, plausibly yields a
   higher realized `C` per macrophage. **This is the lever the planned z=2
   engine change addresses** — not because it removes a *space cap* (that cap is
   not binding; M=176 ≪ pool 438), but because it changes the field/local-IL-10
   geometry that sets the loop gain.
2. **ROS-vs-immune death balance.** Chemokine builds only as `D` grows, but the
   ROS(X) pathway kills the tissue and clears infection before the recruitment
   feedback fully engages, capping the chemokine driver. Stronger immunity then
   *reduces* ODE oxidant X → less ROS death → more survivors, which is why
   uninfected@3.5d = 157 vs the paper's ~2. The two divergences are coupled.

**Recommendation:** do not tune the chemokine/secretion constants. Pursue the z=2
geometry and the ROS-vs-immune death balance together, since they set the same
loop.

## Follow-up: ROS-vs-immune death balance (ablation)

To test whether ROS is over-weighted, the epithelial death was attributed by
mechanism ablation (fig3b bootstrap, cells_per_side=35, seed 0, to t≈0.87 d; at
this point M/K/E ≈ 29/4/2):

| variant | un | inf | dead |
|---|----|-----|------|
| FULL | 889 | 168 | 168 |
| NO_ROS | 518 | **564** | 143 |
| NO_KILLING (NK/CD8 contact+nearby) | 886 | 180 | 159 |
| NO_ROS_NO_KILLING | 515 | 579 | 131 |
| ONLY_APOPTOSIS (no ros/kill/allee) | 473 | 535 | 217 |

- **Removing ROS lets the infection run away** (inf 168 → 564): the ODE oxidant
  `X` (neutrophil-`N`-driven, `dX/dt = b_xn·N/(N+a_xn) − g_xi·I·X − g_xh·H·X −
  mu_x·X`) is the dominant infected-cell clearer — source-faithful (dossier §4b,
  `OxidationAgentModelSteppable`; X is a **global scalar** in the source too, not
  an immune-cell-local field).
- **Removing spatial NK/CD8 killing changes almost nothing** (un 889 → 886):
  at ~4 NK / ~2 CD8, cytotoxic killing is a bystander to epithelial fate. It
  becomes material only at the hundreds-of-cells counts the paper reaches —
  which we under-produce (§ above).

**There is no ROS over-weighting to correct.** Both the chemokine-loop analysis
and this ablation reduce to one root cause: the spatial immune count is ~2–4×
low, so cytotoxic immunity can't influence epithelial fate and ROS necessarily
dominates. The count is set by the chemokine loop gain, whose only source-faithful
lever is the **z=2 geometry** (local-IL-10/chemokine environment + immune-cell
density). Constant tuning is ruled out on both fronts.

## Task 3.4: the Composite driver removes the reserve-pool cap (off-lattice agents)

The chain above explains why `run_full_model`'s macrophage count under-shoots
the fig3b *target*. Separately, `run_full_model` also imposes a *mechanical*
ceiling of its own, independent of loop gain: immune cells are on-lattice CPM
cells activated one-at-a-time from a fixed-size `RECRUIT_RESERVE_TYPE` reserve
pool (`run.py`'s Task-8.4 module note) -- once the pool of dormant reserve
cells is exhausted, further ODE-predicted inflow silently fails to seed, no
matter how strong the recruitment signal gets. Default `run_full_model`
args cap the macrophage population at `n_macrophages(4) + recruit_pool_per_type(6)
= 10`, *regardless of tissue scale*, unless the caller hand-tunes a bigger
pool (as `bootstrap_immune_config` does for the paper-scale fig3b runs).

The `influenza-immune-process-composite` epic (Increments this task closes)
replaces on-lattice immune CPM cells with OFF-LATTICE `ImmuneProcess` agents
(chemotaxing, killing, secreting, and recruiting via the same ODE-driven
rate functions in `recruitment.py`, but with no CPM lattice slot to run out
of -- see `immune_process.py`'s "Recruitment" docstring section). Task 3.4
verifies this end to end:

- **Epithelial parity** (`test_composite_epithelial_parity_reduced_scale`):
  at `cells_per_side=8, steps=15, seed=2, init_infection_frac=0.1`, the
  composite's terminal DEAD count (1) matches `run_full_model`'s reference
  (0) within tolerance (delta 1 <= max(3, int(0.15*1))=3) -- no debugging or
  suppressed immune influence was needed; at this reduced scale/short window
  both drivers stay almost entirely in the infection/no-death regime and the
  one-cell gap is consistent with the composite's immune agents landing one
  proximity kill the reference run's locked-at-seed-count NK/CD8 did not.
- **Uncapped growth** (`test_composite_immune_not_pool_capped`): at
  `cells_per_side=48, steps=35, mcs_per_step=7 (238 raw MCS), seed=1,
  init_infection_frac=0.05`, the composite's macrophage count rises from its
  seeded 6 to a measured max of **14** --
  `[6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8,
  10, 10, 10, 12, 12, 12, 13, 14, 14, 14]` -- comfortably clearing
  `run_full_model`'s default 10-cell ceiling at ANY scale (the pool cap is a
  fixed constant, not a fraction of tissue size). This is the direct proof
  the mechanical reserve-pool cap this epic exists to remove is gone.

  Two smaller scales were tried first and rejected as too noisy to
  demonstrate uncapping cleanly (not because they failed, but because
  reduced-scale recruitment inflow is itself "near-inert", per
  `run_full_model`'s own docstring): `cells_per_side=12/steps=40` (the
  brief's original example) never left the seeded 5-6 macrophage band over
  273 raw MCS (expected inflow < 1 agent at that scale); `cells_per_side=24/
  steps=175` (1044 raw MCS) only reached max=10, inside Poisson-draw
  outflow/inflow noise at that rate. The ODE's macrophage inflow baseline
  term `mu_m*b_m` scales with `b_m ~ cells_per_side**2` but per-record
  simulation cost only grows sub-quadratically with `cells_per_side` in this
  composite (immune agents are off-lattice, so the CPM lattice only needs to
  hold the epithelial patch + a margin box, not on-lattice immune clusters),
  so a bigger `cells_per_side` with fewer raw MCS reaches a decisive signal
  faster than a smaller `cells_per_side` run for longer -- `cells_per_side=48`
  was the scale where this trade-off paid off within a practical test
  runtime (~140s).
