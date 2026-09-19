# Sego 2022 — Global (non-spatial) Price-2015 ODE and its η/θ cellularization coupling

**Authoritative source dossier for Increment 8** of the influenza-sego2022 native reproduction.
Extraction-only document: every equation, constant, and coupling below is transcribed verbatim from
the CC3D *ViralInfectionVTM* source (the codebase Sego et al. 2022, *J. Theor. Biol.* 532:110918,
used to cellularize the Price, Mochan-Keef et al. 2015 in-host influenza ODE). Nothing here is
fabricated; where the source does not state something, that is called out explicitly.

Citations are `path:symbol` with the line region in the extracted source. Line numbers are from the
re-fetched *Source Code* tree (Feb-2021 mtimes, identical layout to the regions cited in
`docs/cc3d-reference/sego2022-parameters.md`).

---

## 1. Source provenance

**How obtained (this increment).** The Elsevier CDN `mmc3.zip` → `"Source Code.7z"` → py7zr
extraction (the procedure recorded in `docs/cc3d-reference/sego2022-source-notes.md`) had already
been performed in this session's scratchpad and was found **cached and intact**; the CDN was not
re-hit. Extraction root used:

```
/private/tmp/claude-502/-Users-eranagmon-code-viva-cpm/8592ff3e-.../scratchpad/cc3d-src/unpacked/Source Code/
```

Files intact (Feb 16 2021), byte sizes match the parameters-dossier's cited regions. No repo files
were fetched into; extraction/reading happened entirely in scratch.

**Files read for this dossier (all under `.../Source Code/`):**

| File | Role |
|---|---|
| `ImmuneModel/ImmuneModelLib.py` | The Price-2015 Antimony model **strings** (spatial-coupling model, full reference ODE, legacy), scaling/homeostatic helpers. |
| `ImmuneModel/ImmuneModelInputs.py` | Immune-model constants: `tot_ec_ODE`, `v0`, local-population ratios, field decay/diffusion, chemotaxis λ, cell sizes, `im_step_size`. |
| `ImmuneModel/ImmuneModelSteppableBasePy.py` | Base steppable; posts/reads the `im_steppable` shared reference. |
| `Simulation/ViralInfectionVTMSteppables.py` | All coupling logic: `ImmuneModelSteppable` (ODE step + recruitment), secretion/killing/death steppables. |
| `Simulation/ViralInfectionVTMModelInputs.py` | Base conversion factors: `s_to_mcs`, `um_to_lat_width`, `exp_cell_diameter`, `cell_diameter`, `cell_volume`, `volume_lm`. |
| `Simulation/ViralInfectionVTM.xml` | `Potts` domain `<Dimensions x="175" y="175" z="2"/>` (active), `500×500` (commented). |
| `nCoVToolkit/nCoVUtils.py` | `hill_equation(val, diss_cf, hill_cf) = 1/(1+(diss_cf/val)^hill_cf)`, `= 0` at `val==0`. |

**Three model strings exist** in `ImmuneModelLib.py` — this is architecturally load-bearing:

- `immune_model_string(...)` (lines **45–229**), model name `FluODE_20vars` — **the model actually
  integrated by the spatial simulation.** In it, the epithelial / immune-cell / spatial-field ODEs
  are **commented out** (they are handled by CPM + reaction-diffusion steppables); only the
  non-spatialized "global" compartments are integrated. This is the object of Increment 8.
- `immune_model_string_ode(...)` (lines **232–410**), model name `FluODE_20varsODE` — **the FULL,
  self-contained Price-2015 20-variable ODE with every equation active.** Used by `SimDataSteppable`
  as a *side-by-side reference/validation* model, not as the spatial driver (see §6).
- `immune_model_string_old(...)` (lines **413–589**) — legacy; **ignore** (superseded).

> **Discrepancy #9 (new).** The "global ODE for Increment 8" is **not** a single full ODE bolted
> onto the spatial model. The source runs a **hybrid**: the spatial-coupling model integrates only
> `{NB, N, T, X, A, B, P, W, G, O}` (+ nearby surrogates), while `{H, I, M, E, K, L, C, F, V, DH}`
> are overwritten every step from CPM/field state, and `{D, DI, Σ1, Σ2, R}` are algebraic
> assignment rules. The complete 20-var ODE (`immune_model_string_ode`) only runs as a *reference*.
> Increment 8 should reproduce the **hybrid coupling**, not integrate all 20 states as the driver.

---

## 2. The Price-2015 global ODE system

### 2.1 State-variable roster (20 states + `DH`)

Names are the single-letter Antimony symbols. Identity **[C]** = confirmed from source (field map
`ImmuneModelSteppable.__spat_fields` / type map `__type_map` / steppable usage); **[I]** = inferred
from the equation structure + paper/spec (source gives only the letter).

| Sym | Identity | Confirmed? | In spatial model integrated or spatialized? |
|---|---|---|---|
| `NB` | Neutrophils, blood/available | [I] spec | **integrated** |
| `N` | Neutrophils, tissue | [I] | **integrated** |
| `M` | Macrophages | [C] `__type_map[MACROPHAGE]='M'` | **spatialized** (CPM count) |
| `L` | IL-10 | [C] `__spat_fields["L"]="il10"` | **spatialized** (field integral) |
| `T` | TNF | [C] spec; enters Σ1 | **integrated** |
| `X` | ROS / oxidative agent | [C] `OxidationAgentModelSteppable` reads `X` | **integrated** |
| `A` | Antibodies | [C] neutralizes virus `g_va·V·A` | **integrated** |
| `C` | Chemokine | [C] `__spat_fields["C"]="chemo"` | **spatialized** (field integral) |
| `I` | Infected epithelial | [C] `__type_map[INFECTEDRELEASING]='I'` | **spatialized** (CPM count) |
| `H` | Uninfected (healthy) epithelial | [C] `__type_map[UNINFECTED]='H'` | **spatialized** (CPM count) |
| `F` | Type-I interferon | [C] `__spat_fields["F"]="type1interferon"` | **spatialized** (field integral) |
| `V` | Virus | [C] `__spat_fields["V"]="Virus"` | **spatialized** (field integral) |
| `DH` | Dead-from-healthy count | [C] `__update_cell_counts` sets `DH=num_dead_uninfected` | **spatialized** (CPM count) |
| `E` | CD8⁺ effector T cells | [C] `__type_map[CD8TCELL]='E'` | **spatialized** (CPM count) |
| `B` | B cells | [C] source of antibody `b_ab·B` | **integrated** |
| `K` | NK cells | [C] `__type_map[NKCELL]='K'` | **spatialized** (CPM count) |
| `P` | APCs (antigen-presenting) | [C] drives CD8 recruit `b_ep`, IFN-I `b_fp·P` | **integrated** |
| `W` | IL-12 | [I] APC→drives IFN-γ | **integrated** |
| `G` | Type-II IFN (IFN-γ) | [I] spec; amplifies APC via `b_pg·G/(a_pg+G)` | **integrated** |
| `O` | CD4⁺ T-helper (Th1) | [I] spec "CD4⁺ T"; produced by APC Hill, drives W & G | **integrated** |

Nearby well-mixed surrogates (integrated in the spatial model only): `M_nb`, `K_nb`, `E_nb`.

### 2.2 Algebraic (assignment `:=`) rules — `ImmuneModelLib.py:immune_model_string`, lines **60–65**

```
D      := tot_cell - H - I                       # total dead epithelial
Sigma1 := a_11*T + a_12*D                         # <-- the "sig_1" stub (see §4)
Sigma2 := Sigma1 + (a_21*V)/(a_22 + V)
DI     := tot_cell - H - I - DH                   # dead-from-infected
R      := F/(a_rf/s_l*s_v + F)                    # resistance rule (spatial-model form)
```

`tot_cell = {tot_ec}` (= `num_epithelial`, templated at build time, line 47).

### 2.3 The FULL 20-variable d/dt system

Verbatim from `ImmuneModelLib.py:immune_model_string_ode` (lines **367–386**), which is the
self-contained reference with **all equations active** (the spatial model's equations are identical
in form; the subset it comments out is noted in the last column). Antimony `-> S; expr` means
`dS/dt = expr`.

```
dNB/dt = b_nt*T/(a_nt + a_nl*L + T) - NB*C*g_nc/(C + a_nc) - mu_n*NB          # [integrated]
dN/dt  = NB*C*g_nc/(C + a_nc) - mu_n*N                                        # [integrated]
dM/dt  = b_mc*C^h_m/(C^h_m + a_mc^h_m) - mu_m*(M - b_m)                       # [spatialized→ recruitment §4]
dL/dt  = b_l*M*Sigma1/(Sigma1 + (g_1*L+g_2)/(L+d_2)) - mu_l*(L - b_lh*(1-R)*H)# [spatialized→ il10 field §4]
dT/dt  = b_t*M*Sigma2/(Sigma2 + (Sigma2 + (g_1*L+g_2)/(L+d_2))*(k_1*L+k_2)/(L+d_1)) - mu_t*T   # [integrated]
dX/dt  = b_xn*N/(N + a_xn) - g_xi*I*X - g_xh*H*X - mu_x*X                      # [integrated]
dA/dt  = b_a + b_ab*B - g_av*A*V - mu_a*A                                     # [integrated]
dC/dt  = b_c*M*Sigma1/(Sigma1 + (g_1*L+g_2)/(L+d_2)) - mu_c*C                 # [spatialized→ chemo field §4]
dI/dt  = g_hv*V*H - g_ix*I*X^h_x/(X^h_x+a_ix^h_x) - g_ik*R*I*K - g_ie*R*I*E - mu_i*(1-R)*I   # [spatialized §4]
dH/dt  = b_h*(1-R)*H*D*(H-theta)/tot_cell - g_hv*V*H - g_hx*H*X^h_x/(X^h_x+a_hx^h_x)          # [spatialized §4]
dF/dt  = b_fi*(1-R)*I + b_fp*P - g_fi*I*F - mu_f*F                            # [spatialized→ t1ifn field §4]
dV/dt  = g_vi*(1-R)*I - g_vh*H*V - g_va*V*A - g_v*V/(1+a_v*V) - mu_v*V        # [spatialized→ Virus field §4]
dDH/dt = g_hx*H*X^h_x/(X^h_x+a_hx^h_x) - b_h*(1-R)*H*DH*(H-theta)/tot_cell    # [spatialized count §4]
dE/dt  = b_ep*P^h_e/(P^h_e+a_ep^h_e) - b_ei*R*I*E - mu_e*E                    # [spatialized→ recruitment §4]
dB/dt  = b_b + b_bp*W*P*(b_0-B) - mu_b*B                                      # [integrated]
dK/dt  = b_kc*C^h_k/(C^h_k+a_kc^h_k) - g_ki*R*I*K - mu_k*(K - b_k)            # [spatialized→ recruitment §4]
dP/dt  = p_0*(g_pv*V/(a_pv+V) + g_pi*DI)*(g_p + b_pg*G/(a_pg+G)) - mu_p*(P - b_p)  # [integrated]
dW/dt  = b_wo*O/(a_wo+O)*P - mu_w*W                                           # [integrated]
dG/dt  = b_go*W/(a_go+W)*O + b_gk*W/(a_gk+W)*K - mu_g*G                       # [integrated]
dO/dt  = b_op*P^h_o/(P^h_o+a_op^h_o) - mu_o*O                                 # [integrated]
```

**Hill exponents** (both strings, lines 48–52 / 254–258): `h_m=3, h_x=2, h_k=2, h_e=3, h_o=2`.

**Spatial-coupling model** (`immune_model_string`, lines 164–187) — the SAME right-hand sides, but
only these are left active (uncommented): `NB, N, T, X, A, B, P, W, G, O`, plus the three nearby
surrogates. The `T` equation there uses `g_2_T`, `d_2_T` in place of `g_2`, `d_2` (unit correction
because `L` is now a spatial-field integral; lines 160–161, 169):
`g_2_T = g_2/s_l*s_v`, `d_2_T = d_2/s_l*s_v`.

### 2.4 Rate constants — verbatim, spatial-coupling model (`immune_model_string`, lines 68–158)

Values are the **calibrated ODE defaults** (the `* s_t / * s_v / * s_l` factors are the
cellularization scalings applied at build time — see §3). Constants with **no** scaling factor are
pure dimensionless rate ratios.

| Const | Value × scaling | Const | Value × scaling |
|---|---|---|---|
| `a_11` | `0.00061091213762049` | `mu_i` | `4.05042485998488 * s_t` |
| `a_12` | `0.0000192792875130612` | `g_vi` | `278.068781202644 * s_t` |
| `a_21` | `1.45266909729275 * s_v` | `g_vh` | `0.0000121539069969014 * s_t / s_l` |
| `a_22` | `1478.77349187226 * s_v` | `g_va` | `0.00109087356190614 * s_t / s_v` |
| `g_1` | `498.217901666102 * s_v` | `g_v` | `197.445019092656 * s_t` |
| `g_2` | `5462.05818489935 * s_v * s_l` | `a_v` | `9.33445943442858 / s_l` |
| `d_2` | `428.253164311871 * s_l` | `mu_v` | `0.412015488642712 * s_t` |
| `b_mc` | `13620.4859808076 * s_t * s_v` | `b_fi` | `0.196756617697923 * s_t` |
| `a_mc` | `574.707608135459 * s_v` | `b_fp` | `0.221363566856883 * s_t * s_l / s_v` |
| `mu_m` | `0.242202409519884 * s_t` | `g_fi` | `0.00181375452827859 * s_t / s_l` |
| `b_m` | `{hs_macro}` (templated) | `mu_f` | `112.230629642229 * s_t` |
| `b_t` | `61.1099298442215 * s_t` | `a_rf` | `53.2223922879035 * s_l` |
| `k_1` | `0.919560555407162` | `b_k` | `{hs_nk}` (templated) |
| `k_2` | `251.401895880346 * s_v` | `b_kc` | `212683.830671333 * s_t * s_v` |
| `d_1` | `57.0468005728478 * s_v` | `a_kc` | `1356.29443420953 * s_v` |
| `mu_t` | `559.245454719045 * s_t` | `g_ki` | `8.83835904194133E-09 * s_t / s_v` |
| `b_l` | `2.32377371988728 * s_t` | `mu_k` | `2.20226570821557 * s_t` |
| `mu_l` | `2.82296559789435 * s_t` | `b_go` | `0.00267078457036687 * s_t` |
| `b_lh` | `0.000463779438894315` | `a_go` | `0.701154807175855 * s_v` |
| `b_c` | `40.203310453199 * s_t` | `b_gk` | `1.61001492527962 * s_t` |
| `mu_c` | `8.99003894588318 * s_t` | `a_gk` | `11.3762919738437 * s_v` |
| `b_nt` | `1831199.07978509 * s_t * s_v` | `mu_g` | `89.9382159197401 * s_t` |
| `a_nt` | `177.556178464038 * s_v` | `p_0` | `6246.50442297409 * s_t * s_v` |
| `a_nl` | `0.0466177133703488` | `g_pv` | `52.361340187773` |
| `g_nc` | `1565.77531103506 * s_t` | `a_pv` | `953.016291300287 * s_v` |
| `a_nc` | `199.135063775527 * s_v` | `g_pi` | `0.0000216210269485218 / s_v` |
| `mu_n` | `0.776781658388633 * s_t` | `g_p` | `0.00973849483459614` |
| `b_xn` | `609.739684621594 * s_t * s_v` | `b_pg` | `0.126601176697398` |
| `a_xn` | `27277.1658907811 * s_v` | `a_pg` | `902.666448958217 * s_v` |
| `g_xi` | `0.0000390938622570566 * s_t / s_v` | `mu_p` | `0.46572454272428 * s_t` |
| `g_xh` | `4.60945452030246E-07 * s_t / s_v` | `b_p` | `{hs_apc}` (templated) |
| `mu_x` | `235.25392576208 * s_t` | `b_ep` | `101318.358506851 * s_t * s_v` |
| `b_h` | `0.0000851174198534486 * s_t` | `a_ep` | `14522.3188063909 * s_v` |
| `theta` | `13217.8105366859 * s_v` | `b_ei` | `2.87556725383114E-07 * s_t / s_v` |
| `g_hv` | `1.41324585239137E-06 * s_t / s_l` | `mu_e` | `0.441073272611722 * s_t` |
| `g_hx` | `29.0380758578089 * s_t` | `b_op` | `317070.082528766 * s_t * s_v` |
| `a_hx` | `15.6825240521184 * s_v` | `a_op` | `14823.8275674019 * s_v` |
| `a_ix` | `0.0611971734649656 * s_v` | `mu_o` | `0.921060461394949 * s_t` |
| `g_ix` | `4.6524853526203 * s_t` | `b_wo` | `0.0092284656047177 * s_t` |
| `g_ik` | `0.0000308183563969158 * s_t` | `a_wo` | `1987.31111667007 * s_v` |
| `g_ie` | `0.000984984039016579 * s_t` | `mu_w` | `0.936640894395778 * s_t` |
| `b_b` | `9.99860060276184 * s_t * s_v` | `b_bp` | `1.05601614882507E-06 * s_t / s_v / s_v` |
| `b_0` | `56777.3326082613 * s_v` | `mu_b` | `0.405022167434195 * s_t` |
| `b_a` | `0.00487789881403722 * s_t * s_v` | `b_ab` | `0.0685940586330701 * s_t` |
| `g_av` | `0.000159387382062868 * s_t / s_v` | `mu_a` | `5.43215036471264 * s_t` |
| `V0` | `{v_ini}` (templated) | | |

> **Discrepancy #10 (new).** The scaling **differs between the two model strings.** In
> `immune_model_string` (spatial, integrated by the sim): `g_2 = ...*s_v*s_l`, `d_2 = ...*s_l`,
> `a_v = .../s_l`, `g_hv = .../s_l`, `g_vh = .../s_l`, `g_fi = .../s_l`, `b_h = ...*s_t` (no `/s_v`),
> `R := F/(a_rf/s_l*s_v+F)`. In `immune_model_string_ode` (reference): the *same* constants use
> `s_v` instead of `s_l` (`g_2 = ...*s_v*s_v`, `d_2 = ...*s_v`, `a_v = .../s_v`, `g_hv = .../s_v`,
> `b_h = ...*s_t/s_v`), and `R := F/(a_rf+F)`. Increment 8 must pick the **spatial** string's
> scalings for the coupled model and keep the ODE string's scalings only for a reference run.

### 2.5 Initial conditions — `immune_model_string`, lines 190–226

```
NB = {hs_neutro_blood}            # default 0
N  = {hs_neutro}                  # default 0
M  = b_m                          # = hs_macro
L  = b_lh*(1-R)*tot_cell
T  = 0
X  = 0
A  = (b_a + b_ab*B)/mu_a
C  = 0
I  = 0
H  = tot_cell
F  = b_fp*b_p/mu_f
V  = V0
DH = 0
E  = b_ep/mu_e * b_p^h_e/(b_p^h_e+a_ep^h_e)
B  = (b_b + b_bp*b_p*W*b_0)/(b_bp*b_p*W + mu_b)
K  = b_k                          # = hs_nk
P  = b_p                          # = hs_apc
W  = b_p/mu_w * O*b_wo/(O+a_wo)
G  = (b_go*W/(a_go+W)*O + b_gk*W/(a_gk+W)*K)/mu_g
O  = b_op/mu_o * b_p^h_o/(b_p^h_o+a_op^h_o)
# homogenized spatial coupling vars, default 0, overwritten each step:
B_ei = 0 ;  G_ki = 0
# nearby surrogates, set at init by init_fresh_immune_model():
M_ratio = 0.0 ; M_nb = 0
K_ratio = 0.0 ; K_nb = 0.0
E_ratio = 0.0 ; E_nb = 0.0
```

**Templated homeostatic values** (`get_homeostatic_pops`, lines 597–616): the calibrated
ODE-scale homeostatic populations, scaled by `pop_scale_factor = tot_ec/2.5e5`:
`hs_macro_ODE = 21347.2655837073`, `hs_nk_ODE = 225.999183404974`,
`hs_apc_ODE = 2064.80195826911`, `hs_neutro_blood_ODE = 0`, `hs_neutro_ODE = 0`.
So `b_m = 0.0049·21347.3 ≈ 104.6`, `b_k = 0.0049·226.0 ≈ 1.107`, `b_p = 0.0049·2064.8 ≈ 10.12`
at the active 175² domain (η = 0.0049; see §3).

**At-init seeding of infection** (`ImmuneModelSteppable.init_fresh_immune_model`, lines 1099–1105):
if `v0 == 0`, `I ← int(num_epithelial·frac_init_infected)` and `H ← H − I`
(`frac_init_infected = 0.05`, `ImmuneModelInputs.py:34`).

---

## 3. η / θ cellularization scaling

The Antimony strings are built by `ImmuneModelSteppable.generate_ode_model_instance`
(`ViralInfectionVTMSteppables.py`, lines **1129–1156**) with four scaling coefficients:

```python
model_string = ImmuneModelLib.immune_model_string(
    tot_ec       = num_epithelial,
    hs_macro/... = get_homeostatic_pops(num_epithelial)[...],
    scale_time   = s_to_mcs / 24 / 60 / 60,                 # s_t
    scale_vol    = ImmuneModelLib.get_pop_scale_factor(num_epithelial),   # s_v  == η
    scale_loc    = 1 / tot_ec_ODE / cell_volume,            # s_l  == θ
    scale_cell_vox = cell_volume,                           # s_c
    v_ini        = get_field_secretor("Virus").totalFieldIntegral())
```

`get_pop_scale_factor(tot_ec) = tot_ec / 2.5e5` (`ImmuneModelLib.py`, lines 592–594);
`tot_ec_ODE = 2.5e5` (`ImmuneModelInputs.py:27`). Base conversions
(`ViralInfectionVTMModelInputs.py`): `s_to_mcs = 60` (line 31), `um_to_lat_width = 2.0` (33),
`exp_cell_diameter = 10.0` (37) ⇒ `cell_diameter = 5` (43), `cell_volume = 25` (45).

### η = `scale_vol` = `s_v` — the **global→patch** (organism→domain) scale

`η = num_epithelial / 250000`. The organism-scale calibrated model represents **250,000** epithelial
cells; the CC3D patch holds `num_epithelial = (L/cell_diameter)²` cells.

| CC3D domain | epithelial cells | η = `s_v` | patch |
|---|---|---|---|
| **175×175** (active in shipped XML, line 19) | 35² = **1225** | **0.0049** | ~0.35 mm |
| 500×500 (commented, line 20) | 100² = **10000** | **0.04** | 1.0 mm (paper's large patch) |
| 150×150 (paper's 0.3 mm, not in XML) | 30² = 900 | 0.0036 | 0.3 mm |

> **CONFIRMED:** the spec's η = 0.04 (1 mm², 250k organism) and η = 0.0049 numbers are exactly the
> `pop_scale_factor` at the 500² (1 mm) and the shipped 175² domains respectively. Note the shipped
> default is 175² (η = 0.0049), **not** a clean 0.3 mm — see existing discrepancy #4. η is **not** a
> hardcoded literal; it is derived from the chosen `Dimensions`.

**Where η multiplies (global→patch):** it scales all *extensive/population* rate constants in the
ODE — every constant carrying `* s_v` in §2.4 (`b_mc, b_kc, b_ep, b_nt, b_xn, p_0, b_op, a_*` Hill
dissociation constants, `theta`, homeostatic `b_m/b_k/b_p`, etc.), so the ODE's populations live at
patch scale rather than organism scale. `get_pop_scale_factor` also scales the homeostatic
populations (`get_homeostatic_pops`).

### θ = `scale_loc` = `s_l` — the **global→site** (organism→lattice-site) scale

`θ = 1 / tot_ec_ODE / cell_volume = 1/(250000·25) = 1.6e-7`.

**Where θ multiplies (global→local/site):** the *concentration/site-density* constants —
`a_rf = 53.2223922879035 * s_l` (resistance half-saturation), `d_2 = ...*s_l`, `g_2 = ...*s_v*s_l`,
and the per-site virus/IFN interaction rates `g_hv, g_vh, g_fi, a_v` (all carrying `/s_l` or `*s_l`).

**Resistance ρ.** The operative per-cell resistance is computed in
`Type1InterferonModelSteppable.update_resistance` (lines 1445–1450):

```python
a_rf = self.im_steppable.get_model_val('a_rf')          # = 53.2223922879035 * s_l
cell.dict[im_resist_key] = t1ifn_cell / (a_rf + t1ifn_cell)   # ρ = f̄ / (a_rf + f̄)
```

where `t1ifn_cell = t1i_secretor.amountSeenByCell(cell)/cell.volume` (`f̄`, local mean type-I IFN).
So the source form is **ρ = f̄ / (θ-scaled a_rf + f̄)**, and every consumer uses **(1 − ρ)** (except
contact-killing — discrepancy #7). The spec's "ρ = 1 − f̄/(θ·a_rf + f̄)" is the same quantity written
as its complement; the θ factor lives inside `a_rf = a_rf_base·s_l`. (`s_c = cell_volume` is passed
but only used to define `s_l`; no other consumer.)

---

## 4. Bidirectional coupling map

Driver: `ImmuneModelSteppable.step` (lines **1070–1086**) each MCS calls `__update_spatial_data()`
(push SPATIAL→ODE), then `self.__rr.timestep()` (integrate), then `update_populations()`
(pull ODE→SPATIAL recruitment). Get/set go through `set_model_val`/`get_model_val` (1319–1327).

### 4a. SPATIAL → ODE (fed in every step)

| ODE var set | From (spatial aggregate) | Code |
|---|---|---|
| `H` | count of `UNINFECTED` CPM cells | `__update_cell_counts`, 1311 |
| `I` | count of `INFECTEDRELEASING` CPM cells | 1311 |
| `M` | `num_immune_by_type(MACROPHAGE)` (local CPM + nearby surrogate) | 1312 |
| `K` | `num_immune_by_type(NKCELL)` | 1312 |
| `E` | `num_immune_by_type(CD8TCELL)` | 1312 |
| `DH` | count of `DYING` cells flagged `im_dead_healthy` | 1313 |
| `V` | `Virus` field total integral `/ dim.z` | `__update_field_integrals`, 1315–1317 |
| `F` | `type1interferon` field integral `/ dim.z` | 1315–1317 |
| `C` | `chemo` field integral `/ dim.z` | 1315–1317 |
| `L` | `il10` field integral `/ dim.z` | 1315–1317 |
| `B_ei` | `(Σ resist over infected cells)·b_ei` | `__update_spatial_data`, 1300–1304 |
| `G_ki` | `(Σ resist over infected cells)·g_ki` | 1305–1306 |

`D := tot_cell−H−I` and `DI := tot_cell−H−I−DH` are then recomputed by the solver from these; `Σ1`,
`Σ2`, `R` follow as assignment rules. So **`V` → `Σ2` → `T` ODE**, and **dead count `D` → `Σ1`**.

### 4b. ODE → SPATIAL (read back each step) — resolves the three stubs

#### (i) `sig_1` = `Sigma1 = a_11·T + a_12·D`   ✅ RESOLVED

`a_11 = 0.00061091213762049`, `a_12 = 0.0000192792875130612` (dimensionless, lines 68–69).
`T` = the **integrated TNF ODE state**; `D = tot_cell − H − I` = **total dead epithelial from spatial
counts**. Read via `get_model_val('Sigma1')` and used to drive **both** chemokine and IL-10
macrophage secretion:

- **Chemokine** (`ChemokineSecretionSteppable.step`, lines 1465–1501). Local, per macrophage cell:
  ```
  sec = b_c · sig_1 / ( sig_1 + (g_1·ĩl10 + g_2)/(ĩl10 + d_2) )      # ĩl10 = local il10 amountSeenByCell/vol
  ```
  Global (well-mixed) term for the "nearby" macrophage surrogate:
  ```
  sec_global = (b_c/(dim.x·dim.y·dim.z)) · num_macro_nearby · sig_1 / ( sig_1 + (g_1·L + g_2_g)/(L + d_2_g) )
  ```
  with `g_2_g = g_2·loc_to_global`, `d_2_g = d_2·loc_to_global`,
  `loc_to_global = scale_vol/scale_loc = s_v/s_l` (line 1476), `L` = total IL-10 ODE var,
  `num_macro = num_immune_nearby_by_type(MACROPHAGE)`.
- **IL-10** (`IL10SecretionSteppable.step`, lines 1523–1571): identical structure with `b_l` in place
  of `b_c`; plus an uninfected-cell homeostatic source `sec = mu_l·b_lh·(1−resist)` (line 1570).

  Constants: `b_c = 40.203310453199*s_t`, `b_l = 2.32377371988728*s_t`, `g_1 = 498.217901666102*s_v`,
  `g_2 = 5462.05818489935*s_v*s_l`, `d_2 = 428.253164311871*s_l` (all read live via
  `get_model_val`). (`b_c`/`b_l` are further multiplied by `dim.z` in the steppable, line 1483/1545.)

> This **resolves the Increment 5/6 `sig_1` stub**: earlier increments treated `sig_1` as a
> placeholder constant because `T` (TNF) and `D` (dead count) had no producer. `T` is now supplied by
> the integrated ODE and `D` by the spatial dead count. **Discrepancy #11 (new, resolution note):**
> the field secretion is NOT a plain `sig_1`-proportional source — it is the **saturating**
> `b·sig_1/(sig_1 + (g_1·L+g_2)/(L+d_2))` Michaelis-type form with IL-10 (`L`) self-inhibition in the
> denominator, split into a per-cell **local** term (IL-10 seen by that macrophage) and a well-mixed
> **global/nearby** term. Any increment that stubbed `sig_1` as a linear chemokine source is wrong
> in functional form, not just value.

#### (ii) Immune-cell RECRUITMENT (Hill inflow / outflow)   ✅ RESOLVED

`ImmuneModelSteppable.update_populations` (1158–1181) pulls per-type inflow/outflow **rates** from
the ODE and stochastically adds/removes CPM cells (`inflow_by_type` Poisson-style draw 1205–1212;
`outflow_by_type` per-cell removal 1214–1221; `local_ratio` splits local vs nearby).

`inflow_rate_by_type` (1183–1194) and `outflow_rate_by_type` (1196–1203), with
`hill(v,a,h)=1/(1+(a/v)^h)`:

| Type | inflow rate | outflow rate | driving field |
|---|---|---|---|
| **Macrophage** | `b_mc·hill(C, a_mc, h_m) + mu_m·b_m` | `mu_m` | **C** (chemokine) |
| **NK** | `b_kc·hill(C, a_kc, h_k) + mu_k·b_k` | `G_ki + mu_k` | **C** (chemokine) |
| **CD8⁺** | `b_ep·hill(P, a_ep, h_e)` | `B_ei + mu_e` | **P** (APC) |

Constants (§2.4): `b_mc=13620.49·s_t·s_v, a_mc=574.71·s_v, mu_m=0.2422·s_t, b_m=hs_macro, h_m=3`;
`b_kc=212683.83·s_t·s_v, a_kc=1356.29·s_v, mu_k=2.2023·s_t, b_k=hs_nk, h_k=2`;
`b_ep=101318.36·s_t·s_v, a_ep=14522.32·s_v, mu_e=0.4411·s_t, h_e=3`.
`G_ki`, `B_ei` are the live homogenized outflow terms pushed in at §4a (resistance-weighted infected
load). Local vs nearby split by `local_ratio_macro=1.0, local_ratio_nk=0.75, local_ratio_cd8=0.75`
(`ImmuneModelInputs.py` 39–44); nearby fraction accumulates in `M_nb/K_nb/E_nb` surrogates.

> **Findings worth flagging (already noted in `sego2022-parameters.md` §11, promoted here to
> Discrepancy #12):** (a) **CD8⁺ recruitment is APC(`P`)-driven, not chemokine-driven** — different
> field from its chemotaxis (which IS chemokine-driven), and (b) CD8⁺ has **no `+ mu_e·b_e`
> homeostatic baseline** term (no `hs_cd8` constant exists), unlike macrophage/NK. Increment 8 must
> not symmetrize the three recruitment laws.

#### (iii) NEARBY-population contact-killing   ✅ RESOLVED

`ContactKillingSteppable` (212–318). Coefficients precomputed in `start` (232–252):
```
pr_cf_nk_nb  = g_ik / pop_scale_factor            # nearby NK   (only if local_ratio_nk < 1)
pr_cf_cd8_nb = g_ie / pop_scale_factor            # nearby CD8  (only if local_ratio_cd8 < 1)
pr_cf_nk_loc  = g_ik * tot_ec_ODE                 # local NK   (surface-contact term)
pr_cf_cd8_loc = g_ie * tot_ec_ODE                 # local CD8
```
`pop_scale_factor = get_pop_scale_factor(num_epithelial) = η`. Per step (261–288), for each infected
cell:
```
pr_cf_nk  = (g_ik/η) · num_nk_nearby              # num_nk_nearby = num_immune_nearby_by_type(NKCELL)
p_death_nk_nearby = ul_rate_to_prob( pr_cf_nk · cell_resist )       # 1 - exp(-x)
```
identically for CD8 with `g_ie`. This is the spec's `g_i/pop_scale · num_nearby` form **exactly**.
`num_immune_nearby_by_type` (1336–1337) = `K_nb`(or `E_nb`) surrogate + queue count — a **well-mixed
scalar**, distinct from the CPM-neighbor surface areas `srf_nk/srf_cd8` used by the LOCAL term
(295–317: `kill_rate = pr_cf_*_loc · srf_* · cell_resist / cell.volume`).
Constants: `g_ik = 0.0000308183563969158·s_t`, `g_ie = 0.000984984039016579·s_t`, `tot_ec_ODE=2.5e5`.

> **Ties to existing discrepancy #7:** both nearby AND local killing multiply by `cell_resist`
> **directly** (not `1−resist`), unlike every other resistance consumer. Preserve verbatim.

#### Other ODE→SPATIAL back-couplings (needed to close the loop, not previously stubbed)

| Effect | Source ODE var | Spatial application | Code |
|---|---|---|---|
| ROS death of uninfected | `X`, via `g_hx·hill(X,a_hx,h_x)` | `pr_death` for `UNINFECTED` cells (flags `im_dead_healthy`) | `OxidationAgentModelSteppable.step`, 1363–1376 |
| ROS death of infected | `X`, via `g_ix·hill(X,a_ix,h_x)` | `pr_death` for infected cells | 1364–1381 |
| Antibody neutralizes virus | `A`, via `g_va·V·A` | added to `Virus` field decay: `virus_decay = min(1, virus_decay_im + g_va·A·dim.z)` | `ViralSecretionSteppable.step`, 200–202 |
| APC-sourced type-I IFN | `P`, via `b_fp·P` | global `t1ifn_secr` term | `Type1InterferonSecretionSteppable.step`, 1418–1422 |
| Virus internalization | `g_hv` | uptake prob `ul_rate_to_prob(g_hv·V̄)` → cell becomes infected | `ViralInternalizationSteppable`, 154–164 |
| Virus release by infected | `g_vi·(1−resist)` | `Virus` secretion per infected cell | `ViralSecretionSteppable`, 205–209 |
| Infected apoptosis | `mu_i·(1−resist)` | death prob per infected cell | `ViralCellDeathSteppable`, 115–128 |
| Allee death/recovery of epithelium | `b_h`, `theta` | surface-area Allee rule; `b_h·tot_ec_ODE`, `theta/num_epithelial` | `RecoverySteppable.step`, 1601–1636 |
| Type-II IFN (`G`) → APC amplification | `G`, via `b_pg·G/(a_pg+G)` factor in `dP/dt` | stays inside ODE (no direct spatial field); reaches spatial layer only through `P`→CD8 recruitment & IFN-I | `immune_model_string_ode` line 383 |

Note: **type-II IFN `G` has no diffusive spatial field** — its only route to the spatial model is
indirect (amplifying APC `P`). Antibody `A` and ROS `X` are likewise **global scalars** applied
uniformly, not fields.

---

## 5. ODE integrator settings

- **Solver:** CC3D free-floating Antimony/SBML → **libRoadRunner** (RoadRunner's default CVODE
  stiff/BDF integrator; deterministic). Instantiated by `add_free_floating_antimony(model_string,
  model_name, step_size=im_step_size)` (`ViralInfectionVTMSteppables.py` 1147–1149), then retrieved
  from `persistent_globals.free_floating_sbml_simulators`.
- **`im_step_size = 1`** (`ImmuneModelInputs.py:24`) — the integrator step size argument.
- **Stepping cadence:** `ImmuneModelSteppable` is registered at `frequency=1`
  (`Simulation/ViralInfectionVTM.py:28`), so `self.__rr.timestep()` runs **once per MCS**
  (one `timestep()` call per MCS, lines 1078–1079). On any solver exception the sim flushes and
  stops (1080–1083).
- **Δt:** `s_to_mcs = 60` s/MCS = **1 minute per MCS** (`ViralInfectionVTMModelInputs.py:31`). The
  ODE constants carry `s_t = scale_time = s_to_mcs/86400 = 6.9444e-4` day, converting the
  Price-2015 **per-day** calibrated rates to per-MCS, so one `timestep()` advances **1 min** of
  biological time. (The reference `FluODE_20varsODE` in `SimDataSteppable` is stepped the same way,
  once per MCS — line 733.)

---

## 6. Discrepancy log (continuing 0–7's #1–#8; new entries #9+)

- **#1–#8** — see `docs/cc3d-reference/sego2022-parameters.md §7` (neighbor order 3; Infected–Macro
  adhesion 20 vs 10; NK–CD8 adhesion 25 vs 10; default dims 175² vs paper patches; v0 sweep not
  enumerated in source; `ViralCellDeathSteppable` "simplified"; contact-killing uses `resist`
  directly; macrophage chemotaxis "moderate" label vs value 5000). Unchanged.
- **#9 (new) — Hybrid, not full-ODE, driver.** The spatially-coupled model
  (`immune_model_string`) integrates only `{NB, N, T, X, A, B, P, W, G, O}` + nearby surrogates;
  `{H, I, M, E, K, L, C, F, V, DH}` are overwritten from CPM/field state each step and their ODE
  lines are commented out. The full 20-var ODE (`immune_model_string_ode`) runs **only** as a
  reference inside `SimDataSteppable`. Increment 8 must implement the hybrid coupling, not a
  free-standing global ODE.
- **#10 (new) — Divergent scaling between the two model strings.** `immune_model_string` scales
  `g_2, d_2, a_v, g_hv, g_vh, g_fi` with `s_l` and sets `b_h = ...*s_t` (no `/s_v`),
  `R := F/(a_rf/s_l*s_v+F)`; `immune_model_string_ode` uses `s_v` for those and `b_h = ...*s_t/s_v`,
  `R := F/(a_rf+F)`. Use the **spatial** string's scalings for the coupled model.
- **#11 (new) — `sig_1` resolves as a saturating IL-10-inhibited secretion, not a linear source.**
  The Increment 5/6 `sig_1` stub is resolved by `Sigma1 = a_11·T + a_12·D` (TNF from ODE, dead count
  from CPM), but the field secretion driven by it is the Michaelis form
  `b·sig_1/(sig_1 + (g_1·L+g_2)/(L+d_2))` split into local + well-mixed-global halves. Match the
  functional form, not just the input.
- **#12 (new) — Asymmetric recruitment laws.** CD8⁺ recruitment is **APC(`P`)-driven** (not
  chemokine) and lacks a homeostatic baseline term; macrophage/NK are chemokine(`C`)-driven with a
  `mu·b` baseline. NK/CD8 outflow carry extra resistance-weighted terms `G_ki`/`B_ei`. Do not
  symmetrize.
- **Resolution note (no discrepancy):** the **nearby-killing** stub is resolved exactly as the
  brief's `g_i/pop_scale·num_nearby` with `pop_scale = η` and `num_nearby` = the well-mixed
  `K_nb`/`E_nb` surrogate — consistent with what Increment 7 recorded as
  `nk.killing.nearby_term_deferred`. No conflict; Increment 8 simply supplies the surrogate
  bookkeeping (`num_immune_nearby_by_type`) and the recruitment loop that maintains `K_nb`/`E_nb`.

---

## 7. For Increment 8 `params.yaml`

Proposed `price_ode:` / `coupling:` block. **Templated values** (`{...}`) are runtime/domain
dependent — record the *formula*, not a fake literal. η shown at both the shipped-175² and the
1 mm²-500² domains. Every key cites `path:symbol`.

```yaml
price_ode:
  # -- provenance --
  source: "ImmuneModel/ImmuneModelLib.py:immune_model_string (lines 45-229)"     # spatial-coupling model
  reference_full_ode: "ImmuneModel/ImmuneModelLib.py:immune_model_string_ode (232-410)"
  model_name: "FluODE_20vars"                          # ImmuneModelLib.py:14
  integrated_states: [NB, N, T, X, A, B, P, W, G, O]   # uncommented eqns, lines 164-187
  integrated_nearby_surrogates: [M_nb, K_nb, E_nb]     # lines 220,223,226
  spatialized_states: [H, I, M, E, K, L, C, F, V, DH]  # set from CPM/fields each step, §4a
  hill_exponents: {h_m: 3, h_x: 2, h_k: 2, h_e: 3, h_o: 2}   # lines 48-52

  # -- integrator --  (§5)
  solver: roadrunner_cvode          # add_free_floating_antimony, ViralInfectionVTMSteppables.py:1147
  step_size: 1                      # ImmuneModelInputs.py:24  im_step_size
  steps_per_mcs: 1                  # ImmuneModelSteppable freq=1, rr.timestep() once/MCS, :1078
  seconds_per_mcs: 60               # ViralInfectionVTMModelInputs.py:31  s_to_mcs  (Δt = 1 min)

  # -- scaling (§3) --
  scaling:
    ode_epithelial_population: 250000        # ImmuneModelInputs.py:27  tot_ec_ODE
    eta_scale_vol_formula: "num_epithelial / 250000"   # get_pop_scale_factor, ImmuneModelLib.py:592
    eta_at_175_domain: 0.0049                # 1225/250000 (shipped XML dims)
    eta_at_500_domain_1mm: 0.04              # 10000/250000 (paper 1mm patch, commented dims)
    theta_scale_loc: 1.6e-7                  # 1/tot_ec_ODE/cell_volume, generate_ode_model_instance:1142
    scale_time_s_t_formula: "s_to_mcs/86400" # = 6.9444e-4 day/MCS
    scale_cell_vox: 25                       # cell_volume, ViralInfectionVTMModelInputs.py:45
    loc_to_global: "s_v/s_l"                 # ChemokineSecretionSteppable:1476

  # -- rate constants: store the calibrated *base* values + their scaling factor tokens from §2.4 --
  # (full table transcribed in §2.4; a_11/a_12 below because they are load-bearing for sig_1)
  a_11: 0.00061091213762049       # ImmuneModelLib.py:68   (dimensionless)
  a_12: 0.0000192792875130612     # ImmuneModelLib.py:69
  a_21: {base: 1.45266909729275, scale: s_v}       # :70
  a_22: {base: 1478.77349187226, scale: s_v}       # :71
  # ... (remaining ~90 constants per §2.4, same {base, scale} convention)

  homeostatic_pops_ode:           # get_homeostatic_pops, ImmuneModelLib.py:597
    macro: 21347.2655837073
    nk:    225.999183404974
    apc:   2064.80195826911
    neutro_blood: 0
    neutro: 0
    # cellularized: multiply each by eta

coupling:
  spatial_to_ode:               # §4a; set each MCS via set_model_val
    H:  count(UNINFECTED)                          # :1311
    I:  count(INFECTEDRELEASING)                   # :1311
    M:  num_immune_by_type(MACROPHAGE)             # :1312
    K:  num_immune_by_type(NKCELL)                 # :1312
    E:  num_immune_by_type(CD8TCELL)               # :1312
    DH: count(DYING & dead_healthy)                # :1313
    V:  Virus_field_integral / dim.z               # :1317
    F:  type1interferon_field_integral / dim.z     # :1317
    C:  chemo_field_integral / dim.z               # :1317
    L:  il10_field_integral / dim.z                # :1317
    B_ei: "sum_infected(resist) * b_ei"            # :1303
    G_ki: "sum_infected(resist) * g_ki"            # :1305

  sig_1:                        # §4b(i) — RESOLVED
    definition: "Sigma1 = a_11*T + a_12*D"         # ImmuneModelLib.py:62
    D: "tot_cell - H - I"                           # :61
    drives: [chemo_secretion, il10_secretion]
    secretion_form_local:  "b * sig_1 / (sig_1 + (g_1*il10_local + g_2)/(il10_local + d_2))"   # :1489/:1551
    secretion_form_global: "b/(dimx*dimy*dimz) * num_macro_nearby * sig_1 / (sig_1 + (g_1*L + g_2_g)/(L + d_2_g))"  # :1499/:1561
    g_2_g: "g_2 * loc_to_global"                    # :1497
    d_2_g: "d_2 * loc_to_global"                    # :1498
    b_chemo: {base: 40.203310453199, scale: s_t}    # b_c, :87
    b_il10:  {base: 2.32377371988728, scale: s_t}   # b_l, :84

  recruitment:                  # §4b(ii) — RESOLVED; hill(v,a,h)=1/(1+(a/v)^h)
    macrophage:
      inflow:  "b_mc*hill(C, a_mc, h_m) + mu_m*b_m"  # :1185-1187
      outflow: "mu_m"                                # :1198
      field: chemo
    nk:
      inflow:  "b_kc*hill(C, a_kc, h_k) + mu_k*b_k"  # :1189-1191
      outflow: "G_ki + mu_k"                         # :1200
      field: chemo
    cd8:
      inflow:  "b_ep*hill(P, a_ep, h_e)"             # :1193  (NO homeostatic baseline term)
      outflow: "B_ei + mu_e"                         # :1202
      field: apc_P
    local_ratios: {macro: 1.0, nk: 0.75, cd8: 0.75}  # ImmuneModelInputs.py:40-44
    seeding_fraction: 0.01                           # new_immune_cell_by_type:1240

  nearby_killing:               # §4b(iii) — RESOLVED
    nk:  "ul_rate_to_prob( (g_ik/eta) * num_nk_nearby  * cell_resist )"   # ContactKilling:250,265,280
    cd8: "ul_rate_to_prob( (g_ie/eta) * num_cd8_nearby * cell_resist )"   # :252,270,285
    local_nk:  "ul_rate_to_prob( g_ik*tot_ec_ODE * srf_nk  * cell_resist / cell_volume )"   # :238,305
    local_cd8: "ul_rate_to_prob( g_ie*tot_ec_ODE * srf_cd8 * cell_resist / cell_volume )"   # :240,313
    g_ik: {base: 0.0000308183563969158, scale: s_t}  # :107
    g_ie: {base: 0.000984984039016579, scale: s_t}   # :108
    note: "multiplies cell_resist DIRECTLY, not (1-resist) — discrepancy #7"

  resistance:                   # §3
    form: "rho = f_bar / (a_rf + f_bar)"             # update_resistance:1450 ; consumers use (1-rho)
    a_rf: {base: 53.2223922879035, scale: s_l}       # :120  (theta scaling)

  other_ode_to_spatial:         # §4b (loop-closing, not previously stubbed)
    ros_death_uninfected: "ul_rate_to_prob(g_hx*hill(X,a_hx,h_x))"   # OxidationAgent:1369
    ros_death_infected:   "ul_rate_to_prob(g_ix*hill(X,a_ix,h_x))"   # :1370
    antibody_virus_clearance: "virus_decay += g_va * A * dim.z"      # ViralSecretion:202
    apc_type1ifn_source: "b_fp * P"                                  # T1IFNSecretion:1420
    infected_apoptosis: "ul_rate_to_prob(mu_i*(1-resist))"          # ViralCellDeath:125
    allee_death: "b_h*tot_ec_ODE ...  theta/num_epithelial"          # Recovery:1601-1602
    type2_ifn_G_route: "indirect only — amplifies dP/dt via b_pg*G/(a_pg+G); no diffusive field"  # ODE line 383
```

---

### Appendix — key symbol → CC3D field / cell-type map (confirmed)

`__spat_fields` (`ViralInfectionVTMSteppables.py:1024-1027`): `V→Virus, F→type1interferon,
C→chemo, L→il10`. `__type_map` (`1043-1047`): `H→UNINFECTED, I→INFECTEDRELEASING, M→MACROPHAGE,
E→CD8TCELL, K→NKCELL`. `hill_equation` (`nCoVToolkit/nCoVUtils.py:42`):
`1/(1+(diss_cf/val)^hill_cf)`, `0` at `val==0`. `ul_rate_to_prob` (`ImmuneModelLib.py:619`):
`1 - exp(-rate)`.
