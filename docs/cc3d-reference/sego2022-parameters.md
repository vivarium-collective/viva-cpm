# Sego et al. 2022 — influenza model parameters (for faithful reimplementation)

Source paper: T.J. Sego, Ericka D. Mochan, G. Bard Ermentrout, James A. Glazier, "A multiscale
multicellular spatiotemporal model of local influenza infection and immune response,"
*J. Theor. Biol.* 532 (2022) 110918. Tables 1-4, Fig. 2, Sec. 2.2-2.4. Local copy:
`workspace/references/papers/sego-2022-influenza.pdf`.

Source code: the CC3D 4.2.2 **ViralInfectionVTM** package, supplementary "Source Code" (`mmc3.zip`
→ `Source Code.7z`), retrieved per `docs/cc3d-reference/sego2022-source-notes.md`. Re-fetched to a
scratch dir for this task (not vendored in-repo — re-run the retrieval commands in that file to
reproduce). Per that note's fidelity convention and this task's brief, **literal constants from the
source are preferred over the paper's Tables 1-4 wherever they overlap**; where source and paper
disagree, both are recorded explicitly rather than silently picking one (see "Source vs. paper
discrepancies" at the end).

Following the `demo-parameters.md` convention: every value below carries its source (paper
table/eq, or source file + Python symbol, or "CC3D-typical, marked"). Where the exact constant
can't be pinned, the *validated behavior* is the fidelity criterion, not the literal number.

---

## 1. CPM constants (Tables 3-4 + source)

| Parameter | Value | Source |
|---|---|---|
| Lattice width (`um_to_lat_width`) | 2 µm/lattice site | Paper Table 4; source `ViralInfectionVTMModelInputs.py:um_to_lat_width = 2.0` |
| Experimental cell diameter (`exp_cell_diameter`) | 10 µm | Source `ViralInfectionVTMModelInputs.py:exp_cell_diameter = 10.0` |
| Unitless cell diameter (`cell_diameter`) | 5 lattice sites (`10/2`) | Source, derived: `exp_cell_diameter / um_to_lat_width` |
| Cell footprint / volume constraint v_c | 25 sites = 100 µm² (5×5) | Paper Table 3 (`v_c = 100 µm²`, "chosen for an average epithelial cell diameter of 10 µm"); source `cell_volume = cell_diameter**2 = 25` |
| Volume multiplier λ_v (`volume_lm`) | **9** | Paper Table 3; source `ViralInfectionVTMModelInputs.py:volume_lm = 9`, applied per cell as `cell.lambdaVolume = volume_lm` in `ViralInfectionVTMSteppables.py:1291` |
| Intrinsic random motility / CPM temperature H* | **10.0** | Paper Table 4 ("Graner and Glazier, 1992"); source XML `ViralInfectionVTM.xml: <Potts><Temperature>10.0</Temperature>` |
| Time step Δt | **1 min/step** (`s_to_mcs = 60`) | Paper Table 4; source `ViralInfectionVTMModelInputs.py:s_to_mcs = 1 * 60  # s/mcs` |
| Seeding fraction | **1%** | Paper Table 4 ("chosen for efficient local immune response") — see §6, no literal `0.01` constant found in the base `ModelInputs` files; used in immune-cell recruitment logic in the steppables, not a single named constant. Marked CC3D-typical/paper-sourced, to confirm exact call site if a later increment needs the literal recruitment code path. |
| Neighbor order | **3** (both `Potts` and `Contact` plugins) | Source XML: `<Potts><NeighborOrder>3</NeighborOrder>`, `<Plugin Name="Contact">...<NeighborOrder>3</NeighborOrder>` — **DISCREPANCY vs. this task's brief**, see §7. Paper Tables 3-4 do not state a neighbor order explicitly. |
| Boundary conditions | **Periodic in-plane (X, Y); Neumann (zero-flux) out-of-plane (Z)** | Paper text (Sec. 2.4): "Neumann and periodic conditions were applied to boundaries parallel and orthogonal, respectively, to the epithelial sheet." Source XML confirms for both `Potts` (`Boundary_x/y = Periodic`, no Z entry) and every diffusion field (`Plane Axis="X"/"Y"` → `<Periodic/>`; `Plane Axis="Z"` → `<ConstantDerivaive .../>` value 0.0, i.e. zero-flux/Neumann). The domain is a thin slab (`z=2` lattice sites = 4 µm) representing the epithelial sheet in the XY plane, so "orthogonal to the sheet" = the Z faces = Neumann, "parallel to the sheet" boundaries = the X/Y side walls = periodic. |
| Cell-type freeze | Epithelial types (`Uninfected`, `Infected`, `InfectedReleasing`, `Dying`) are `Freeze=""` (static lattice footprint, do not move under CPM copy attempts); immune types (`Macrophage`, `NKcell`, `CD8Tcell`) are not frozen (motile) | Source XML `<Plugin Name="CellType">`. Not explicitly called out in the paper's Tables but consistent with the paper's description of a fixed epithelial sheet with motile immune cells (Fig. 2). |

## 2. Adhesion (Contact) J matrix (Table 3 + source)

Paper Table 3 gives only immune-cell-facing values, collapsing epithelial subtypes into
"uninfected / infected / dead" and immune subtypes into one "immune" category:

| Pair (paper's categories) | J (paper) |
|---|---|
| Uninfected – immune | 20 |
| Infected – immune | 10 |
| Dead – immune | 20 |
| Homotypic immune | 25 |
| Heterotypic immune | 10 |

The brief flagged epithelial-epithelial and epithelial-medium J as **a gap in Table 3**. The
source XML's `Contact` plugin gives the full, literal 8×8 matrix (Medium, Uninfected, Infected,
InfectedReleasing, Dying, Macrophage, NKcell, CD8Tcell), which resolves that gap directly:

| Type1 \ Type2 | Medium | Uninfected | Infected | InfectedReleasing | Dying | Macrophage | NKcell | CD8Tcell |
|---|---|---|---|---|---|---|---|---|
| **Medium** | 10.0 | 25.0 | 25.0 | 25.0 | 25.0 | 10.0 | 10.0 | 10.0 |
| **Uninfected** | — | 5.0 | 5.0 | 5.0 | 5.0 | 20.0 | 20.0 | 20.0 |
| **Infected** | — | — | 5.0 | 5.0 | 5.0 | **20.0** | 10.0 | 10.0 |
| **InfectedReleasing** | — | — | — | 5.0 | 5.0 | 10.0 | 10.0 | 10.0 |
| **Dying** | — | — | — | — | 5.0 | 20.0 | 20.0 | 20.0 |
| **Macrophage** | — | — | — | — | — | 25.0 | 10.0 | 10.0 |
| **NKcell** | — | — | — | — | — | — | 25.0 | **25.0** |
| **CD8Tcell** | — | — | — | — | — | — | — | 25.0 |

Source: `Source Code/Simulation/ViralInfectionVTM.xml`, `<Plugin Name="Contact">`.

**Resolved gap:** epithelial-epithelial J = **5.0** uniformly (all Uninfected/Infected/
InfectedReleasing/Dying pairs); epithelial-Medium J = **25.0** uniformly (except Medium-Medium =
10.0). Both values are literal source constants, not CC3D-typical placeholders.

**Two cells bolded above disagree with the paper's collapsed Table 3** — see §7 for the
discrepancy writeup (Infected–Macrophage = 20 vs. paper's "infected–immune = 10"; NKcell–CD8Tcell
= 25 vs. paper's "heterotypic immune = 10").

## 3. Diffusion coefficients & lengths (Table 3, from de Jong et al. 2006 — cross-checked against source)

Paper's stated diffusion length convention: **diffusion length L = √(D/decay)**, in units of cell
diameters (10 µm each). The source computes it the other way — decay rates are literal ODE-model
constants (calibrated, from Price et al. 2015 via `ImmuneModelLib`), and diffusion coefficients are
*derived* as `D = decay · L²` for an assumed diffusion length in cell diameters:

| Field | ODE decay rate (source, `ImmuneModelInputs.py`) | decay (s⁻¹) | Diffusion length (cell diam.) | D = decay·L² (µm²/s) | Paper Table 3 D (µm²/s) | Match |
|---|---|---|---|---|---|---|
| Extracellular virus v | `virus_decay_ODE = 0.412015488642712` /day | 4.7687e-6 | 5 (`exp_virus_dl`) | **0.01192** | 0.0119 | ✓ |
| Type-I IFN f | `t1ifn_decay_ODE = 112.230629642229` /day | 1.2987e-3 | 2 (`exp_t1ifn_dl`) | **0.5195** | 0.520 | ✓ |
| Chemokines c | `chemo_decay_ODE = 8.99003894588318` /day | 1.0401e-4 | 10 (`exp_chemo_dl`) | **1.0401** | 1.04 | ✓ |
| IL-10 l | `il10_decay_ODE = 2.82296559789435` /day | 3.2673e-5 | 10 (`exp_il10_dl`, "assumed same as chemokines" per paper text) | **0.3267** | 0.327 | ✓ |

All four literal source-derived values reproduce the paper's Table 3 diffusion coefficients to the
stated 3 significant figures — strong cross-validation between source and paper for this section.

Unitless (lattice/MCS) values actually fed to the CC3DML `GlobalDiffusionConstant`/
`GlobalDecayConstant` fields at runtime (overriding the XML skeleton's placeholder defaults of
`0.1`/`1E-6` via `id`-tagged substitution, e.g. `virus_dc`, `virus_decay`): unitless decay = decay(s⁻¹)·60
(s/MCS); unitless D = D(µm²/s)·60/`um_to_lat_width²` = D·60/4. Computed: virus decay ≈ 2.861e-4/MCS,
D ≈ 0.1788 lattice²/MCS; type-I IFN decay ≈ 0.07792/MCS, D ≈ 7.792; chemokines decay ≈ 6.240e-3/MCS,
D ≈ 15.601; IL-10 decay ≈ 1.960e-3/MCS, D ≈ 4.901. Source: `ImmuneModelInputs.py`
(`virus_decay_im`, `virus_dc_im`, `t1ifn_decay`, `t1ifn_dc`, `chemo_decay`, `chemo_dc`, `il10_decay`,
`il10_dc`), derived from `BaseInputs.s_to_mcs` / `BaseInputs.um_to_lat_width`.

Explicit-scheme stability note (per `demo-parameters.md`'s DiffusionSolverFE convention):
Δt·D·(2·ndim) ≤ 1 must still hold for the unitless per-MCS values above; not separately checked
in the source (CC3D's `ReactionDiffusionSolverFE`/`DiffusionSolverFE` sub-step internally if
needed).

**Increment 2 — exact virus-field literals confirmed against source** (re-fetched per
`sego2022-source-notes.md`, `ImmuneModel/ImmuneModelInputs.py`, `um_to_lat_width=2.0`,
`exp_cell_diameter=10.0`, `s_to_mcs=60`):

```
exp_virus_decay_im = virus_decay_ODE(0.412015488642712 /day) / 86400 = 4.768697785216574e-06 /s
exp_virus_dl        = exp_cell_diameter * 5.0 = 50.0 um   (source `exp_virus_dl`; = 5 cell diameters)
exp_virus_dc_im      = exp_virus_decay_im * exp_virus_dl**2 = 0.011921744463041435 um^2/s
virus_decay_im       = exp_virus_decay_im * s_to_mcs = 0.00028612186711299444  /MCS  (unitless decay)
virus_dc_im          = exp_virus_dc_im * s_to_mcs / um_to_lat_width**2 = 0.17882616694562153 lattice^2/MCS
```

CONFIRMED against `ImmuneModelInputs.py:virus_dc_im` and `virus_decay_im` — matches the ≈0.1788 /
≈2.861e-4 approximations stated above to full precision (those were already correct; this records
the exact literals). `diffusion_length_cell_diam = exp_virus_dl / exp_cell_diameter = 5` (source
`exp_virus_dl`). Recorded verbatim in `params.yaml`'s `virus:` section.

**Increment 3 — exact type-I IFN field literals confirmed against source** (re-fetched per
`sego2022-source-notes.md`, `ImmuneModel/ImmuneModelInputs.py`, same conversion constants as above):

```
exp_t1ifn_decay = t1ifn_decay_ODE(112.230629642229 /day) / 86400 = 0.001298965620859132 /s
exp_t1ifn_dl     = exp_cell_diameter * 2.0 = 20.0 um   (source `exp_t1ifn_dl`; = 2 cell diameters)
exp_t1ifn_dc     = exp_t1ifn_decay * exp_t1ifn_dl**2 = 0.5195862483436527 um^2/s
t1ifn_decay      = exp_t1ifn_decay * s_to_mcs = 0.07793793725154792  /MCS  (unitless decay)
t1ifn_dc         = exp_t1ifn_dc * s_to_mcs / um_to_lat_width**2 = 7.793793725154791 lattice^2/MCS
```

CONFIRMED against `ImmuneModelInputs.py:t1ifn_dc` and `t1ifn_decay` — matches the ≈7.792 / ≈0.07792
approximations in §3's table above to full precision (those were already correct; this records the
exact literals). `diffusion_length_cell_diam = exp_t1ifn_dl / exp_cell_diameter = 2` (source
`exp_t1ifn_dl`). Recorded verbatim in `params.yaml`'s `ifn:` section.

## 4. Chemotaxis λ and functional form (Table 3, source confirms + extends)

| Cell type | Field | λ_c (paper Table 3 & source) | Source |
|---|---|---|---|
| Macrophage | Virus | **5,000** | Paper Table 3 ("chosen for moderate chemotaxis..." — NB paper text says "moderate" for macrophage-virus, "strong" for macrophage-virus is NOT what the table caption says; see raw table image, macrophage row says "Chosen for moderate chemotaxis according to typical field values"); source `ImmuneModelInputs.py: chemotaxis_v_macro = 5E3` |
| NK cell | Chemokines | **5,000** | Paper Table 3 ("chosen for strong chemotaxis..."); source `chemotaxis_v_nk = 5E3` |
| CD8+ T cell | Chemokines | **10,000** | Paper Table 3 ("chosen for strong chemotaxis..."), CD8+ sensitivity "twice that of NK cells" per paper text (Sec. 2.3); source `chemotaxis_v_cd8 = 10E3` |

All three λ values match exactly between paper Table 3 and source literals — no discrepancy here.

**Paper's stated functional form** (Sec. 2.3 text, saturating, evaluated at the moving cell's
center of mass y'): `Δ_chem = λ_c · c(y') / (1 + c_CM(σ))` — a *saturating* chemotactic
sensitivity, not a plain linear gradient-following term.

**Source's literal implementation** (`ViralInfectionVTMSteppables.py`, `ChemotaxisSteppable`,
lines ~340-351):
```python
concentration = field[cell.xCOM, cell.yCOM, cell.zCOM]   # field value AT the cell's own COM
cd.setLambda(chemotax_val / (1.0 + concentration))        # dynamically-scaled lambda
```
i.e. the source does **not** implement a custom saturating ΔH formula directly; instead it
re-sets the *coefficient* fed to CC3D's standard `Chemotaxis` plugin every MCS, as
`λ_effective(t) = λ_c / (1 + c(cell's own COM))`. CC3D's built-in `Chemotaxis` plugin then applies
its standard linear pixel-copy energy `ΔH_chem = -λ_effective · (c(x_dest) - c(x_source))` (see
`demo-parameters.md`'s Chemotaxis section) using this saturating, per-cell, per-MCS-updated λ. This
is functionally close to — but not textually identical to — the paper's stated
`λ_c·c(y')/(1+c_CM(σ))` form (the source's saturation denominator uses the concentration at the
cell's own current COM rather than a distinct `c_CM(σ)` domain-center-of-mass-normalized term).

**Known gap for viva-cpm (per `demo-parameters.md`):** viva-cpm's engine currently implements only
the *plain linear* CC3D form `ΔH_chem = -λ·(c_dest - c_source)` with a **fixed** λ per cell type —
it does not yet support a per-MCS-recomputed, concentration-dependent λ (the saturating mechanism
above). Closing this gap (dynamic λ = λ_c/(1+c) recomputed each step) is deferred to a later
increment; Increment 0 should use the fixed-λ CC3D form as a documented approximation and treat the
saturating behavior as a validated-behavior target once implemented.

## 5. Field PDE source/decay terms (Table 1, general form)

General reaction-diffusion form for each field z̄ ∈ {v̄ (virus), f̄ (type-I IFN), c̄ (chemokines),
l̄ (IL-10)}, per paper Table 1 (cellularized from the ODE model via a global scaling coefficient η
and local scaling coefficient θ — see §6):

    ∂_t z̄ = (D_z ∂_i² − q) z̄ + r

where **q** is a field-specific decay rate and **r** is a field-specific source rate, both built
from ODE-model rate parameters (subscripted `b_*`, `g_*`, `μ_*`) evaluated at each lattice site's
cell type τ(σ,t), scaled by the binary indicator function 𝓑(x,y) (1 when x=y, else 0), cell volume
|𝒱(σ,t)|, and — for virus decay/uptake — the per-cell resistance ρ(s,t) (§ below). All symbols
with subscripts are ODE-model parameters (Price et al. 2015 calibration, carried through
`ImmuneModel/ImmuneModelLib.py`'s Antimony/SBML model string) — this task does not re-derive every
individual ODE rate constant (e.g. `g_hv`, `g_vi`, `mu_i`, `b_h`, `a_rf`, `g_ik`, `g_ie`, ...); they
are available verbatim in `ImmuneModelLib.py` for a later task that needs the full ODE parameter
table.

Field-by-field summary (decay q / source r structure, per paper Table 1 and cross-checked against
the steppable that implements each in the source):

| Field | Decay q (structure) | Source r (structure) | Implementing steppable |
|---|---|---|---|
| Virus v̄ | Natural decay `μ_v` + mucociliary/antibody clearance `+ b_a·a/(a_v+A)` + uninfected-cell uptake `+ 𝓑(τ,Ĥ)·g_hv/|𝒱|` | Release by infected/releasing cells, resistance-scaled: `𝓑(τ,Î)·g_vi·(1−ρ)/|𝒱|` | `ViralSecretionSteppable` (secretion `g_vi*(1-resist)`), `ViralInternalizationSteppable` (uptake `g_hv`) |
| Type-I IFN f̄ | Natural decay `μ_f` + uptake by infected cells `+ 𝓑(τ,Î)·g_fi/(θ|𝒱|)` | Release by infected/releasing cells `𝓑(τ,Î)·g_fp/|𝒱|` (basal), plus APC-driven amplification | `Type1InterferonSecretionSteppable` |

**Increment 3 — exact `g_fp`/`g_fi`/resistance `a_rf` literals confirmed against source**, same
cellularized-instance convention as Increment 2's `g_hv`/`g_vi` (`ImmuneModel/ImmuneModelLib.py`'s
`immune_model_string()`, instantiated with real `scale_time=s_t`, `scale_loc=s_l`,
`scale_vol=s_v=η` from `ViralInfectionVTMSteppables.py:1134-1144` — **not** the unscaled
`immune_model_string_old()` copy of the same raw `/day` symbols):

```
s_t (day/MCS) = s_to_mcs / 86400 = 6.944444e-4          (as in Increment 2)
s_l            = 1 / tot_ec_ODE / cell_volume = 1.6e-7   (as in Increment 2)
s_v            = η = num_epithelial / tot_ec_ODE          (get_pop_scale_factor; SCENARIO-dependent,
                                                             0.0049 for 0.3mm patch, 0.04 for 1.0mm)
```

**Naming clarification (do not conflate the two "basal IFN" terms):** no symbol literally spelled
`g_fp` exists anywhere in the source. `Type1InterferonSecretionSteppable.step` has a code *comment*
"Add global secretion by APCs: `g_fp*P`", but the Antimony parameter it actually reads is `b_fp`
(`sec_amount_global = num_APCs * get_model_val('b_fp')`), applied as a spatially-**uniform**
boundary secretion (via 8 XML elements `t1ifn_secr0`..`t1ifn_secr7`), **population(η)-scaled**
(`b_fp = 0.221363566856883 * s_t * s_l / s_v`) — i.e. this is the table row's "plus APC-driven
amplification" clause above, not the basal per-infected-cell term, and is **not** recorded as a
`params.yaml` key (out of scope; flagged for whoever wires immune recruitment/amplification later).
The literal per-infected-cell, resistance-gated term that structurally matches the table's basal
`𝓑(τ,Î)·g_fp/|𝒱|` clause is source symbol **`b_fi`** (`Type1InterferonSecretionSteppable.step`:
`sec_amount = (1-resist)*b_fi`, for `cell in cell_list_by_type(INFECTED, INFECTEDRELEASING)`,
`secreteInsideCell(cell, sec_amount/cell.volume)`) — this is what `params.yaml`'s `ifn.secretion_g_fp`
records, under the brief's `g_fp` key name but the source's `b_fi` symbol:

```
b_fi (raw, ODE-calibrated) = 0.196756617697923 /day   (ImmuneModelLib.py immune_model_string(): "b_fi = 0.196756617697923 * s_t")
b_fi_percell_per_mcs = b_fi_raw * s_t = 0.00013663654006800208
  -> * self.dim.z(=2) [Type1InterferonSecretionSteppable.step: `b_fi = get_model_val('b_fi') * self.dim.z`]
  -> secretion_g_fp (as literally applied) = 0.00013663654006800208 * 2 = 0.00027327308013600416
```

`g_fi` (uptake by infected cells, matches the brief's `uptake_g_fi` key and the table's decay-side
`g_fi/(θ|𝒱|)` clause directly — population-independent, no `s_v`/η factor):

```
g_fi (raw, ODE-calibrated) = 0.00181375452827859 /day   (ImmuneModelLib.py immune_model_string(): "g_fi = 0.00181375452827859 * s_t / s_l")
g_fi_percell_per_mcs = 0.00181375452827859 * s_t / s_l = 7.87219847343138
  -> * self.dim.z(=2) [Type1InterferonSecretionSteppable.step: `up_amount = get_model_val('g_fi') * self.dim.z`]
  -> uptake_g_fi (as literally applied) = 7.87219847343138 * 2 = 15.74439694686276
```

Both `secretion_g_fp` (`b_fi`) and `uptake_g_fi` (`g_fi`) bake in the same `dim.z=2` slab-thickness
factor as Increment 2's `secretion_g_vi` — the same z=1-vs-z=2 convention note applies (halve if the
viva-cpm IFN field is single-layer). Recorded verbatim in `params.yaml`'s `ifn:` section; wiring
uptake into the field-consumer logic is deferred to a later task per the brief.
| Chemokines c̄ | Natural decay `μ_c` | Release regulated by TNF and dead-cell presence `𝓑(τ,D̂)·b_p·ρ'` | `ChemokineSecretionSteppable` |
| IL-10 l̄ | Natural decay `μ_l` | Release by macrophages regulated by TNF and dead cells, resistance-scaled `(1−ρ)` | `IL10SecretionSteppable` (`sec_amount = mu_l*b_lh*(1-resist)`) |

Paper's auxiliary form: `ζ = 1 + η(g₁l + g₂c)/((a₁1²+a₁2)(1+a₀2))` (mean-field IL-10/chemokine
cross-term; Table 1 footnote) and `z̄̄ = (1/|𝒱(σ,t)|)∫_{𝒱(σ,t)} z̄ dV` (mean cellular field
measurement). Transcribed as given in the paper; not independently re-derived from source (the
exact Antimony rate-law strings live in `ImmuneModelLib.py`'s `immune_model_string`/
`immune_model_string_ode` functions for a later task that needs bit-exact ODE rate laws).

**Rate → probability conversion used throughout the source** (applies to every stochastic event
below): `Pr(event in one MCS) = 1 − exp(−rate)` — source `ImmuneModelLib.py:ul_rate_to_prob`.

## 6. Type-transition rates (Table 2)

Paper Table 2 gives dense multi-term Hill-function expressions per transition (virus/ROS Hill
terms with exponents n_H, n_I, thresholds, Allee terms, contact-killing γ terms). Rather than
risk mis-transcribing the exact printed subscripts from the PDF, the **structure actually run** is
taken from the source steppables (ground truth per the source-notes fidelity convention), with the
paper's transition names/symbols preserved:

| Transition | Symbol | Rate structure (source) | Implementing steppable |
|---|---|---|---|
| Infection Ĥ→Î | — | `rate = g_hv · v̄(s)` (mean local extracellular virus concentration at the cell, no saturation term in source) → `Pr = 1−exp(−rate)`; on success cell type → `InfectedReleasing` directly (source has no separate "Infected-not-yet-releasing" transition state distinct from `InfectedReleasing` for this event — see note below) | `ViralInternalizationSteppable.do_cell_internalization` |
| Uninfected death Ĥ→D̂ (Allee, `a_D`) | a_D(s) | Cellularized Allee-effect death: `rate = b_h·(1−ρ)·srf_area_oi·(srf_thresh−srf_uninfected)/srf_area_total²`, only when local uninfected surface fraction < threshold `θ_local = θ_ODE/num_epithelial` scaled by neighbor surface area | `RecoverySteppable` (docstring: "Implements cellularized Allee effect death and recovery") |
| Infected death Î→D̂ (viral-load term) | — | `rate = μ_i·(1−ρ)` (docstring: "Virally-induced apoptosis, **simplified version**" — see §7 discrepancy note: this omits the paper's printed ROS/viral-load Hill saturation term) | `ViralCellDeathSteppable` |
| Infected death Î→D̂ (NK contact-killing, γ term) | γ(s; g_ik, K̂, H₀A_s) | Two-part: population-scaled "nearby" killing `rate = (g_ik·tot_ec_ODE/pop_scale)·num_nk_nearby·ρ`, plus local surface-contact killing `rate = g_ik·tot_ec_ODE·srf_nk·ρ/|𝒱|` | `ContactKillingSteppable` |
| Infected death Î→D̂ (CD8+ contact-killing, γ term) | γ(s; g_ie, Ê, H₀A_s) | Same structure as NK, with `g_ie` and CD8+ population/surface | `ContactKillingSteppable` |
| Recovery D̂→Ĥ (Allee, `a_H`) | a_H | Cellularized Allee-effect recovery, structurally symmetric to `a_D` above (swap uninfected↔dying roles in the surface-fraction calculation) | `RecoverySteppable` |

**Exact numeric literals for the Allee death/recovery and infected-death rows above** (`μ_i`,
`b_h`, the surface threshold `θ_local`, and the confirmed recovery-branch form/coefficient) are
in the "Increment 4" subsection below, after the `g_hv`/`g_vi`/`g_fp`/`g_fi`/`a_rf` derivations
for the other rows in this table.

**Increment 2 — exact `g_hv`/`g_vi` literals confirmed against source.** Both are ODE-calibrated
rate constants defined in `ImmuneModel/ImmuneModelLib.py`'s `immune_model_string()` (the
**cellularized**, per-cell Antimony model generator instantiated with real `scale_time`/
`scale_vol`/`scale_loc` from `ViralInfectionVTMSteppables.py:1134-1144` — not the separate
organism-scale `immune_model_string_ode()` copy of the same symbols, which uses a different
`s_l`/`s_v` scaling and is not what the CPM steppables read via `get_model_val`):

```
g_hv (raw, ODE-calibrated) = 1.41324585239137E-06 /day   (per unit local-scale virus)
g_vi (raw, ODE-calibrated) = 278.068781202644 /day        (per InfectedReleasing cell)

s_t (day/MCS)  = s_to_mcs / 86400 = 60/86400 = 6.944444e-4          (ViralInfectionVTMSteppables.py:1140/1544)
s_l            = 1 / tot_ec_ODE / cell_volume = 1/(250000*25) = 1.6e-7  (…:1142, cell_volume=25 sites^2)

infection_g_hv (per-MCS) = g_hv * s_t / s_l = 0.006133879567670876
  -> used directly, no further factor: ViralInternalizationSteppable.do_cell_internalization,
     `rate = g_hv * viral_amount_com`, `viral_amount_com = secretor.amountSeenByCell(cell)/cell.volume`,
     `Pr = 1 - exp(-rate)` (nCoVUtils/ImmuneModelLib.ul_rate_to_prob). On success: cell type ->
     InfectedReleasing directly.

secretion_g_vi (per-MCS-per-cell, pre z-factor) = g_vi * s_t = 0.19310332027961388
secretion_g_vi (as literally applied)           = 0.19310332027961388 * dim.z(=2) = 0.38620664055922777
  -> ViralSecretionSteppable.step: `g_vi = self.im_steppable.get_model_val('g_vi') * self.dim.z`,
     then `sec_amount = g_vi * (1 - resist)`, `secretor.secreteInsideCell(cell, sec_amount/cell.volume)`.
     `dim.z = 2` is the lattice's fixed z-thickness (`ViralInfectionVTM.xml <Dimensions .../>`,
     z=2 in both the active 175x175x2 config and the commented-out 500x500x2 config — i.e. fixed
     across both 0.3mm/1.0mm patch scenarios, not itself a scenario parameter).
```

Both `infection_g_hv` and `secretion_g_vi` are **independent of domain size / epithelial population**
(no `s_v`/η scaling factor in either formula) — a single canonical value applies across the 0.3mm
and 1.0mm patch scenarios. Recorded verbatim (with full derivation) in `params.yaml`'s `virus:`
section; Task 2.1/2.2 must pick and document consistently whether to use the pre- or post-`dim.z`
`secretion_g_vi` value depending on whether the viva-cpm virus field models a z=2 slab like the
source or a single-layer (z=1) field — see the CONVENTION note in `params.yaml`.

**Initial conditions** (source `ImmuneModel/ImmuneModelInputs.py`): `v0_ODE` is one of three
hand-picked, mutually-commented options — `0` (default active; means "seed by
`frac_init_infected` instead"), `50` ("non-lethal according to calibrated ODE model"), `500`
("lethal according to calibrated ODE model"); `frac_init_infected = 0.05` is the fraction of
epithelial cells randomly set to `InfectedReleasing` at t=0 when `v0_ODE=0`
(`CellsInitializerSteppable`, `num_to_infect = int(num_epithelial * frac_init_infected)`). When
`v0_ODE>0`, `v0 = v0_ODE/tot_ec_ODE` (virus per epithelial cell at ODE scale) and the virus field
is seeded **uniformly** at every lattice site to `v0_sites = v0 * num_epithelial / (dim.x*dim.y)`
(`FieldInitializerSteppable`, `Simulation/ViralInfectionVTMSteppables.py`) — i.e. "initial viral
load" is a total-virus quantity converted to a per-site field concentration by the actual seeded
cell count and domain area, not a literal field value itself.

**Increment 4 — exact `mu_i`/`b_h`/Allee-surface-threshold/recovery literals confirmed against
source** (re-fetched per `sego2022-source-notes.md`; `ImmuneModel/ImmuneModelLib.py`
`immune_model_string()`, `Simulation/ViralInfectionVTMSteppables.py`
`ViralCellDeathSteppable`/`RecoverySteppable`):

```
s_t (day/MCS) = s_to_mcs / 86400 = 6.944444e-4          (as in Increments 2-3)
tot_ec_ODE     = 250,000                                  (ImmuneModelInputs.py, as in §6/Cellularization)

mu_i (raw, ODE-calibrated)  = 4.05042485998488 /day
  (ImmuneModelLib.py immune_model_string(): "mu_i = 4.05042485998488 * s_t")
mu_i_per_mcs = mu_i_raw * s_t = 0.002812795041656167     (population-independent; no s_l/s_v factor)
  -> ViralCellDeathSteppable.step: `mu_i = get_model_val('mu_i')`,
     `pr_death = ul_rate_to_prob(mu_i * (1 - cell_resist))` for cell in
     cell_list_by_type(INFECTED, INFECTEDRELEASING). No dim.z factor.

b_h (raw, ODE-calibrated) = 0.0000851174198534486 /day
  (ImmuneModelLib.py immune_model_string(): "b_h = 0.0000851174198534486 * s_t")
b_h_model_val (cellularized) = b_h_raw * s_t = 5.910931934267264e-08
b_h (as literally applied)   = b_h_model_val * tot_ec_ODE = 0.01477732983566816
  -> RecoverySteppable.step: `b_h = self.im_steppable.get_model_val('b_h') * tot_ec_ODE`
     (population-independent: the *only* population-dependent factor cancels, see below).

theta_ODE_raw = 13217.8105366859   (ImmuneModelLib.py immune_model_string(): "theta = 13217.8105366859 * s_v")
theta_local (srf_threshold) = get_model_val('theta') / num_epithelial
                             = (theta_ODE_raw * s_v) / num_epithelial,  s_v = eta = num_epithelial/tot_ec_ODE
                             = theta_ODE_raw / tot_ec_ODE               (num_epithelial cancels algebraically)
                             = 13217.8105366859 / 250000 = 0.0528712421467436
  -> SCENARIO-INDEPENDENT: identical for the 0.3mm patch (eta=0.0049, num_epithelial=1225) and the
     1.0mm patch (eta=0.04, num_epithelial=10000) — verified numerically, both reduce to
     0.0528712421467436. `srf_thresh` (the per-cell absolute cutoff each step) = theta_local *
     srf_area_total (the cell's total epithelial-neighbor interface area that step).
```

**RecoverySteppable — EXACT death and recovery branches** (literal, from
`Simulation/ViralInfectionVTMSteppables.py`, class `RecoverySteppable.step`, ~lines 1594-1636):

```python
b_h = self.im_steppable.get_model_val('b_h') * tot_ec_ODE
theta = self.im_steppable.get_model_val('theta') / self.num_epithelial
for cell in self.cell_list_by_type(self.UNINFECTED, self.DYING):
    # ... accumulate srf_area_total, srf_area_oi, srf_uninfected over epithelial neighbors ...
    # (type_oi = DYING if cell.type == UNINFECTED else UNINFECTED)
    srf_thresh = theta * srf_area_total
    resist = cell.dict[ImmuneModelLib.im_resist_key]

    if cell.type == self.UNINFECTED and srf_uninfected < srf_thresh:
        death_rate = b_h * (1 - resist) * srf_area_oi * (srf_thresh - srf_uninfected) / srf_area_total ** 2.0
        pr_death = ImmuneModelLib.ul_rate_to_prob(death_rate)
        if random.random() < pr_death:
            cell.dict[ImmuneModelLib.im_dead_healthy_key] = True
            cell.type = self.DYING
    elif cell.type == self.DYING and srf_uninfected > srf_thresh:
        revive_rate = b_h * (1 - resist) * srf_area_oi * (srf_uninfected - srf_thresh) / srf_area_total ** 2.0
        pr_revive = ImmuneModelLib.ul_rate_to_prob(revive_rate)
        if random.random() < pr_revive:
            cell.type = self.UNINFECTED
```

**IMPORTANT — no separate `a_H` symbol exists in the source.** `grep -rn "a_H\|a_h\b"` over the
full re-fetched source tree finds nothing besides `b_h`/`a_hx` (an unrelated ODE Hill-threshold
constant used only inside the *non-cellularized* ODE H-dynamics, never referenced by
`RecoverySteppable`). Both the death branch (`UNINFECTED`→`DYING`) and the recovery branch
(`DYING`→`UNINFECTED`) read the **same** local variable `b_h`, computed once at the top of
`step()`. The paper's printed Table 2 uses `a_D`/`a_H` as two distinct symbol names for the
death/recovery Allee terms, but the source implements both with one shared coefficient — this
transcription doc's own §6 table (row: "Recovery D̂→Ĥ (Allee, `a_H`)... structurally symmetric to
`a_D`") anticipated the symmetry; this increment confirms the coefficient itself is literally
identical (`b_h`), not merely structurally analogous. `params.yaml`'s `allee:` section records
`b_h` once and an explicit `a_H_equals_b_h: true` flag rather than inventing a distinct `a_H`
numeric value.

Recorded verbatim in `params.yaml`'s new `cell_death:` and `allee:` sections.

Local immune-type inflow/outflow (macrophage, NK, CD8+ recruitment — Table 2's "local immune type"
rows) use Hill-equation recruitment (`nCoVUtils.hill_equation`) driven by chemokines C / APCs P,
plus Poisson-like stepwise addition (`inflow_by_type`) and removal (`outflow_by_type`) per MCS,
implemented in `ImmuneModelSteppableBasePy`/`ViralInfectionVTMSteppables.py`. Exact per-parameter
values (`b_mc`, `a_mc`, `h_m`, `b_kc`, ...) are ODE-model constants from `ImmuneModelLib.py`, not
re-derived here (out of scope for the CPM-level parameter transcription; flagged for `params.yaml`
authors to pull directly from `ImmuneModelLib.py` if the immune-recruitment sub-model is
implemented in this increment).

## Cellular viral resistance ρ

Paper's printed form (Sec. 2.2): `ρ = ρ(s,t) = 1 − f̄'/(θ·a_f + f̄')`, where `a_f` is an ODE model
parameter and `f̄` is the mean type-I IFN in the cell's domain.

Source's literal implementation (`ViralInfectionVTMSteppables.py`, `Type1InterferonModelSteppable.
update_resistance`):
```python
t1ifn_cell = t1i_secretor.amountSeenByCell(cell) / cell.volume   # local mean IFN, f̄(s,t)
cell.dict[im_resist_key] = t1ifn_cell / (a_rf + t1ifn_cell)       # resist = f̄ / (a_rf + f̄)
```
i.e. source computes `resist = f̄/(a_rf + f̄)` (a plain Hill/saturation function, **no leading
`1 −`, no θ scaling factor**), and every consumer of this value in the source (`ViralCellDeathSteppable`,
`RecoverySteppable`, `ChemokineSecretionSteppable`/`IL10SecretionSteppable`) applies it as
`(1 − resist)` wherever the paper's formulas use `ρ` in a protective role — e.g. infected-death rate
`μ_i·(1−resist)`, Allee death rate `b_h·(1−resist)·...`, IL-10 secretion `mu_l·b_lh·(1−resist)`.
Net effect: `(1 − resist_source) ≈ ρ_paper` when the paper's own `1 −` and source's `(1−resist)`
wrapper are composed — i.e., the source's stored `resist` variable behaves like the paper's `f̄/(θa_f+f̄)`
term (the thing being subtracted from 1), not like `ρ` itself. **One exception:** `ContactKillingSteppable`
multiplies `kill_rate` by `cell_resist` *directly* (not `1−cell_resist`) — see §7, flagged as a
possible naming/sign inconsistency worth double-checking against the paper's contact-killing γ term
before `params.yaml` encodes it.

**Increment 3 — exact `a_rf` literal confirmed against source**, same cellularized-instance
convention as the `g_fp`/`g_fi` derivation above (`ImmuneModelLib.py`'s `immune_model_string()`,
`scale_loc=s_l=1.6e-7`):

```
a_rf (raw, ODE-calibrated) = 53.2223922879035   (ImmuneModelLib.py immune_model_string(): "a_rf = 53.2223922879035 * s_l")
a_rf (cellularized)        = 53.2223922879035 * s_l = 8.51558276606456e-06
```

Applied with **no** `dim.z` factor — `Type1InterferonModelSteppable.update_resistance` reads
`a_rf = self.im_steppable.get_model_val('a_rf')` as-is, then `cell.dict[im_resist_key] =
t1ifn_cell / (a_rf + t1ifn_cell)` where `t1ifn_cell = t1i_secretor.amountSeenByCell(cell) /
cell.volume` — i.e. `a_rf` is in the same `s_l`-scaled internal field-concentration units as that
per-cell mean-IFN measurement, unlike `secretion_g_fp`/`uptake_g_fi` above which both get an
additional `dim.z=2` factor. Recorded verbatim (with `formula: "resist = f_bar/(a_rf + f_bar)"`) in
`params.yaml`'s `resistance:` section.

## Cellularization scaling (Table 4 + Sec. 2.3 text + source)

| Quantity | Value | Source |
|---|---|---|
| ODE-model organism-level epithelial population (`tot_ec_ODE`) | **250,000** cells | Paper text ("Scaling was performed by epithelial cell population... 250 k"); source `ImmuneModelInputs.py: tot_ec_ODE = 2.5e5`. **Used only for scaling** — the spatial CPM domain holds far fewer cells (see η below); it is not simulated directly. |
| Local scaling coefficient θ | **4×10⁻⁸ µm⁻²** | Paper text: "calculated from the total number of epithelial cells according to the ODE model (250 k) and cell volume constraint v_c (Lucas et al., 2020)" — derivable as `θ = 1/(tot_ec_ODE · v_c) = 1/(250,000 · 100 µm²) = 4×10⁻⁸ µm⁻²`. Not found as a single literal Python constant of this name in the base `ModelInputs` files inspected (the ODE-level `theta` in `ImmuneModelLib.py`, e.g. `13217.81...`, is a **different** quantity — the ODE Allee population threshold, not this PDE-cellularization scaling coefficient; do not conflate the two "theta" symbols). Marked paper-sourced; formula shown, not independently re-derived from a source literal. |
| Global scaling coefficient η | **0.0049** (0.3 mm patch) / **0.04** (1.0 mm patch) | Paper text: "the ratio of the number of epithelial cells in the simulation domain to the ODE model parameters" — 1,225/250,000 ≈ 0.0049 (0.3 mm, ~150×150 lattice ÷ 25 sites/cell ≈ 900 cells at minimum, up to 1225 depending on packing) and 10,000/250,000 = 0.04 (1.0 mm, 500×500 lattice ÷ 25 sites/cell = 10,000 cells exactly). Not found as a single named Python literal in the inspected source files (likely computed per-run from actual seeded cell count via `num_epithelial`, e.g. `theta = get_model_val('theta') / self.num_epithelial` in `RecoverySteppable` — this divides the ODE-level `theta` by the *actual* epithelial count each run, which is the dynamic analogue of a fixed η). Marked paper-sourced. |
| Seeding fraction | **1%** | Paper Table 4. See §1 note — not isolated as a single literal constant in the inspected files. |
| Local immune-population fractions | Macrophage **100%**, NK **75%**, CD8+ **75%** | Paper Table 4; source `ImmuneModelInputs.py`: `local_ratio_macro = 1.0`, `local_ratio_nk = 0.75`, `local_ratio_cd8 = 0.75` — **exact match**. |
| Domain population note | ODE model's 250k is an **organism-scale** epithelial population used only to calibrate the scaling coefficients θ/η above; the spatial CPM domain itself holds **~10,000 cells at 1.0 mm²** (500×500 lattice ÷ 25 sites/cell) and **~900-1,225 cells at 0.3 mm²** (≈150×150 lattice ÷ 25 sites/cell, depending on packing efficiency at the domain boundary) | Paper text (Sec. 2.4) + arithmetic from lattice width (§1) and patch sizes (§8) |

## 7. Source vs. paper discrepancies (flagged, not silently resolved)

1. **Neighbor order.** Source XML (`Potts` and `Contact` plugins) uses `NeighborOrder = 3`. This
   task's brief (and `demo-parameters.md`'s general CC3D convention) expected second-order Manhattan
   (`NeighborOrder = 2`). The paper's Tables 3-4 do not state a neighbor order. **Use `3` (source
   ground truth) for `params.yaml`**, not `2`.
2. **Infected–Macrophage adhesion.** Source: `Infected`–`Macrophage` J = **20.0** (matches the
   *Dead*-immune and *Uninfected*-immune value, not the *Infected*-immune value). Paper Table 3
   states a single "Infected–immune = 10" for all three immune types. The source's
   `InfectedReleasing`–`Macrophage` = 10.0 **does** match the paper. Net: the paper's collapsed
   category appears accurate for `InfectedReleasing` but not for the pre-release `Infected` subtype,
   which is 2× the paper's stated value against macrophages specifically (NK/CD8+ contacts with
   `Infected` do match the paper's 10). Use the literal per-subtype source matrix in §2 for
   `params.yaml`.
3. **NKcell–CD8Tcell adhesion.** Source: **25.0** (equal to the homotypic-immune value). Paper
   Table 3 states "heterotypic immune = 10" (and the Macrophage–NK / Macrophage–CD8+ source values
   do match 10.0). So NK↔CD8+ specifically is *not* heterotypic-reduced in the source, unlike the
   other two heterotypic immune pairs. Use the literal source matrix.
4. **Default lattice dimensions.** Source XML's active `<Dimensions x="175" y="175" z="2"/>`
   (350×350 µm ≈ 0.35 mm) matches **neither** paper patch size exactly; a commented-out
   `<Dimensions x="500" y="500" z="2"/>` (1000×1000 µm = 1.0 mm, matching the paper's larger patch
   and giving exactly 10,000 cells at 25 sites/cell) is present but inactive. The 0.3 mm patch
   (150×150 lattice, 300×300 µm) is not present in this file at all — likely set by hand-editing
   this XML per scenario run, outside what this package ships as a labeled config. `params.yaml`
   should treat 0.3 mm / 1.0 mm as explicit scenario overrides of `Dimensions`, not assume the
   shipped XML default encodes either paper scenario.
5. **Initial viral load sweep.** Paper text describes an initial-viral-load sweep of
   {1, 10, 100, 1000, 10000} (per this task's brief) and an initial-infection-fraction sweep of
   {0.001, 0.005, 0.01, 0.05}. Source `ImmuneModelInputs.py` only exposes `v0_ODE` as one of three
   hand-picked, mutually-commented options — `0` (use `frac_init_infected` instead), `50`
   ("non-lethal according to calibrated ODE model"), `500` ("lethal according to calibrated ODE
   model") — with `0` (i.e., infected-fraction-seeded) active by default, and `frac_init_infected =
   0.05` (matching the top of the paper's infection-fraction sweep, not literally one of the four
   listed values' full set — 0.05 is present, 0.001/0.005/0.01 are not encoded as alternatives in
   this file). The full sweep values from the paper are **scenario/run configuration**, not present
   as an enumerated list in the base `ModelInputs` — `params.yaml` should encode the paper's sweep
   values directly as scenario parameters rather than expect them literally in this source file.
6. **`ViralCellDeathSteppable` is an explicitly "simplified version".** Its own docstring says
   "Virally-induced apoptosis, simplified version" and implements only a flat `μ_i·(1−resist)` rate,
   omitting the viral-load/ROS Hill-saturation term that the paper's printed Table 2 "Infected
   death" row appears to include alongside the contact-killing γ terms. Treat the paper's fuller
   printed formula as the *conceptual* target and the source's flat term as the literal, *simplified*
   behavior actually run to produce the paper's results — both should be recorded; `params.yaml`
   should mirror the source's simplified form unless a later increment intentionally reinstates the
   fuller Hill term.
7. **`ContactKillingSteppable` multiplies by `cell_resist` directly**, not `(1 − cell_resist)`,
   unlike every other consumer of the resistance value (`ViralCellDeathSteppable`, `RecoverySteppable`,
   `ChemokineSecretionSteppable`, `IL10SecretionSteppable`, all of which use `(1 − resist)`). This
   may be intentional (contact-killing efficacy could plausibly scale differently than viral-load
   effects) or an inconsistency; flagged rather than resolved — see the "Cellular viral resistance ρ"
   section above.
8. **Macrophage chemotaxis "moderate" vs. paper table value 5,000.** Paper Table 3's
   Source/Justification column literally reads "Chosen for moderate chemotaxis according to typical
   field values" for macrophage–virus (vs. "strong" for the other two rows), yet the numeric λ
   (5,000) is identical to NK's ("strong"). Not a numeric discrepancy (source confirms 5,000 exactly)
   — noted only because the qualitative label and the numeric value don't obviously match each
   other within the paper table itself.

## 8. Scenarios (Table 4 text + source)

| Scenario axis | Values |
|---|---|
| Initial viral load sweep | 1, 10, 100, 1000, 10000 (paper) — source only ships 0/50/500 as example `v0_ODE` picks, see discrepancy §7.5 |
| Initial infection fraction sweep | 0.001, 0.005, 0.01, 0.05 — source default `frac_init_infected = 0.05` |
| Patch sizes | 0.3 mm and 1.0 mm square domains, discretized at 2 µm/site (150×150 and 500×500 lattice sites respectively) |
| Boundary conditions | Neumann (zero-flux, out-of-plane/Z) + periodic (in-plane, X/Y) — see §1 |
| Run length / termination | Paper Sec. 3: "executed for two weeks of simulation time at most," terminated early on a "determined lethal scenario" (all epithelial cells dead) or an "assumed non-lethal scenario" (uninfected and total extracellular virus < 0.001, "several orders of magnitude less than typical values during infection") |
| Replicates | Multiple stochastic replicas per initial condition/parameter set (paper Sec. 3, exact replicate count not stated in the excerpted pages) |

---

## Fidelity convention

Per `demo-parameters.md` and the source-notes fidelity convention: values in §1-§6 and §8 marked
with an explicit source-file citation are **exact literals**, safe for `params.yaml` to encode
verbatim. Values marked "CC3D-typical" or "paper-sourced, formula only" (θ, η, seeding fraction,
individual ODE rate constants beyond those listed) should be treated as **structurally correct but
not independently re-derived from a source literal** — where `params.yaml` needs a specific number
for one of these, pull it directly from `ImmuneModelLib.py`/`ImmuneModelInputs.py` at implementation
time rather than trusting an approximation here. The eight discrepancies in §7 are not resolved in
favor of either the source or the paper; `params.yaml`'s author should pick one explicitly (this doc
recommends the source's literal values, per the source-notes convention, except where noted) and
record which was chosen.
