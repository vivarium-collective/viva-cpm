"""``SystemicODEProcess`` -- a process-bigraph `Process` wrapping the global
Price-2015 ODE (`price_ode.GlobalODE`) so the systemic species (oxidant X,
antibody A, APC P, ...) and the recruitment driver rates become process-
bigraph stores the other processes (epithelium, immune recruitment) consume.

Mirrors exactly how `run.py`'s `_ode_couple` (~run.py:2098-2128) feeds and
steps the ODE:
  - assembles the input dict `H,I,M,K,E,DH,V,F,C,L,B_ei,G_ki` from `state`
    (no nearby-surrogate correction here -- unlike the CPM driver, this
    Process's M/K/E inputs are the uncapped agent-level counts already, so
    they are used as given, no `+ state["M_nb"]` surrogate addition);
  - steps `self.ode` with `dt_seconds` and persists the returned ODE state
    across updates (`self.state`);
  - derives the six recruitment driver rates via `recruitment.*` at the
    freshly-integrated `C` (from the input) and `P` (from the new ode_state);
  - derives the dynamic `sig_1 = a_11*T + a_12*D` with `D = num_epithelial -
    H - I` (verbatim from `run.py`'s `_ode_couple`:
    ``a_11*state["T"] + a_12*(tot_cell - H - I)``), using the freshly-
    integrated `T` and the *input* `H,I`.
"""
from __future__ import annotations

from process_bigraph import Process

from bigraph_schema.contract import ProcessContract

from . import price_ode, recruitment

_ODE_INPUT_INT_KEYS = ("H", "I", "M", "K", "E", "DH")
_ODE_INPUT_FLOAT_KEYS = ("V", "F", "C", "L", "B_ei", "G_ki")


class SystemicODEProcess(Process):
    """Global (non-spatial) Price-2015 immune ODE as a process-bigraph Process.

    Advances the 10 systemic species (NB, N, T, X, A, B, P, W, G, O) of the
    Mochan-Keef/Price-2015 influenza immune model one step per ``update``, with
    the spatialized species (H, I, M, K, E, L, C, V, DH) and the derived
    algebraic variables held fixed over the step as inputs. Each step it also
    emits the six per-type recruitment driver rates (macrophage/NK/CD8+
    inflow/outflow) and the T-cell activation signal ``sig_1``. It is the
    swappable global-coupling alternative to ``run.py``'s per-MCS CPM driver
    loop (``_ode_couple``), wiring the ODE's outputs into stores the epithelium
    and immune-recruitment processes consume.
    """

    description = (
        "Global Price-2015 influenza immune ODE. Integrates the 10 systemic "
        "species {NB, N, T, X, A, B, P, W, G, O} with scipy LSODA over "
        "dt_seconds (default one 60 s MCS); the spatial species {H,I,M,K,E,L,"
        "C,V,DH} and derived vars {D, Sigma1, Sigma2, DI, R} are frozen as "
        "inputs. Emits the integrated state, the six Hill recruitment driver "
        "rates, and the T-activation signal sig_1 = a_11*T + a_12*D."
    )

    contract = ProcessContract(
        summary=(
            "Advances the ten systemic species of the Mochan-Keef/Price-2015 influenza "
            "immune-response ODE one biological step per update, holding the spatialized "
            "cell/field species fixed as inputs. It integrates the scaled per-MCS system "
            "with scipy LSODA, then derives and emits the six immune-cell recruitment "
            "driver rates and the T-cell activation signal that the spatial epithelium and "
            "recruitment processes consume — the swappable global-coupling alternative to "
            "the per-MCS CPM driver loop."
        ),
        description=(
            "State is persisted across updates (self.state): each update assembles the "
            "fixed inputs, steps price_ode.GlobalODE over dt_seconds / seconds_per_mcs "
            "model-time units, clamps tiny numerical negatives to 0, and carries the three "
            "nearby surrogates (M_nb, K_nb, E_nb, d/dt = 0) through unchanged. Rate "
            "constants arrive pre-resolved and patch-cellularized from "
            "price_ode.resolve_constants (time scale s_t = s_to_mcs/86400, population "
            "scale s_v = num_epithelial/250000, density scale s_l = 1/250000/cell_volume)."
        ),
        inputs={
            "H": "Healthy (uninfected) epithelial-cell count; frozen over the step. Feeds D = num_epithelial - H - I, DI, and the oxidant loss term g_xh*H*X.",
            "I": "Infected epithelial-cell count; frozen. Feeds D, DI, and the oxidant loss term g_xi*I*X.",
            "M": "Macrophage count; frozen. Drives T-cell activation (b_t*M*Sigma2/denom_T in dT/dt).",
            "K": "NK-cell count; frozen. Drives the b_gk*W/(a_gk+W)*K term of dG/dt.",
            "E": "CD8+ effector-cell count; accepted for interface symmetry with the per-MCS driver but not consumed by the integrated RHS in this build.",
            "DH": "Dead (lysed) epithelial-cell count; frozen. Feeds DI = num_epithelial - H - I - DH, the APC viral/debris drive.",
            "V": "Extracellular virus level (field); frozen. Feeds Sigma2, antibody neutralization g_av*A*V, and APC production g_pv*V/(a_pv+V).",
            "F": "Interferon (IFN) level (field); frozen. Passed to derived_inputs to set resistance R = F/(a_rf' + F); R is computed but does not enter the 10 integrated ODEs in this build.",
            "C": "Chemokine level (field); frozen. Drives neutrophil extravasation g_nc*NB*C/(C+a_nc) and the macrophage/NK recruitment-inflow Hill terms.",
            "L": "Lymph/antigen-transport level (field); frozen. Appears in the neutrophil (a_nl*L) and T-cell (g_1, g_2_T, k_1, k_2, d_1, d_2_T) denominators.",
            "B_ei": "Live resistance-weighted infected-load term for CD8+; frozen. Added to the CD8+ outflow rate (B_ei + mu_e).",
            "G_ki": "Live resistance-weighted infected-load term for NK; frozen. Added to the NK outflow rate (G_ki + mu_k).",
        },
        outputs={
            "ode_state": "overwrite[map[float]] of the 10 freshly-integrated systemic species {NB, N, T, X, A, B, P, W, G, O} (price_ode.INTEGRATED_STATES).",
            "recruit_drivers": "overwrite[map[float]] of the six per-type recruitment rates: macro_inflow, nk_inflow, cd8_inflow, macro_outflow, nk_outflow, cd8_outflow (per-MCS).",
            "sig_1": "overwrite[float] T-cell activation signal sig_1 = a_11*T + a_12*D (freshly-integrated T, input-derived D = num_epithelial - H - I).",
        },
        config={
            "consts": "Pre-resolved, patch-cellularized rate-constant dict from price_ode.resolve_constants (concrete floats plus the underscore-keyed helpers _s_t/_s_v/_s_l; typed 'tree' so those keys survive Composite realization).",
            "num_epithelial": "Total epithelial-cell count N_ep; sets the population scale s_v and the D/DI conservation totals.",
            "dt_seconds": "Wall-clock seconds advanced per update (default 60.0, one MCS); divided by seconds_per_mcs to get the per-MCS integration span.",
        },
        math=[
            "dNB/dt = b_nt*T/(a_nt + a_nl*L + T) - g_nc*NB*C/(C + a_nc) - mu_n*NB",
            "dN/dt  = g_nc*NB*C/(C + a_nc) - mu_n*N",
            "dT/dt  = b_t*M*Sigma2 / [ Sigma2 + (Sigma2 + (g_1*L + g_2_T)/(L + d_2_T)) * (k_1*L + k_2)/(L + d_1) ] - mu_t*T",
            "dX/dt  = b_xn*N/(N + a_xn) - g_xi*I*X - g_xh*H*X - mu_x*X",
            "dA/dt  = b_a + b_ab*B - g_av*A*V - mu_a*A",
            "dB/dt  = b_b + b_bp*W*P*(b_0 - B) - mu_b*B",
            "dP/dt  = p_0*(g_pv*V/(a_pv + V) + g_pi*DI)*(g_p + b_pg*G/(a_pg + G)) - mu_p*(P - b_p)",
            "dW/dt  = b_wo*O/(a_wo + O)*P - mu_w*W",
            "dG/dt  = b_go*W/(a_go + W)*O + b_gk*W/(a_gk + W)*K - mu_g*G",
            "dO/dt  = b_op*P^{h_o}/(P^{h_o} + a_op^{h_o}) - mu_o*O",
            # derived algebraic inputs, held fixed over the ODE step:
            "D = N_ep - H - I",
            "Sigma1 = a_11*T + a_12*D",
            "Sigma2 = Sigma1 + a_21*V/(a_22 + V)",
            "DI = N_ep - H - I - DH",
            # recruitment driver rates (hill clamped to 0 at v<=0):
            "hill(v,a,h) = 1/(1 + (a/v)^h)",
            "macro_inflow = b_mc*hill(C,a_mc,h_m) + mu_m*b_m",
            "nk_inflow = b_kc*hill(C,a_kc,h_k) + mu_k*b_k",
            "cd8_inflow = b_ep*hill(P,a_ep,h_e)",
            # recruitment outflow rates:
            "macro_outflow = mu_m",
            "nk_outflow = G_ki + mu_k",
            "cd8_outflow = B_ei + mu_e",
            "sig_1 = a_11*T + a_12*D",
        ],
        symbols={
            "NB": "blood neutrophil pool (cells, patch-cellularized count)",
            "N": "tissue (extravasated) neutrophils (cells)",
            "T": "activated T cells (cells)",
            "X": "oxidant / reactive-oxygen species (model concentration units, patch-scaled)",
            "A": "antibody (model concentration units)",
            "B": "antibody-producing B-cell pool (cells)",
            "P": "antigen-presenting cells, APC (cells)",
            "W": "APC-derived T-helper/cytokine mediator (model concentration units)",
            "G": "effector cytokine (granzyme/IFN-gamma-like) mediator (model concentration units)",
            "O": "upstream APC-driven cytokine (model concentration units)",
            "H": "healthy epithelial-cell count (cells, input)",
            "I": "infected epithelial-cell count (cells, input)",
            "M": "macrophage count (cells, input)",
            "K": "NK-cell count (cells, input)",
            "E": "CD8+ effector-cell count (cells, input; unused by the RHS here)",
            "DH": "dead epithelial-cell count (cells, input)",
            "D": "susceptible/dead-pool complement, N_ep - H - I (cells)",
            "DI": "dead-infected cells, N_ep - H - I - DH (cells)",
            "V": "extracellular virus level (field units, input)",
            "F": "interferon/IFN level (field units, input)",
            "C": "chemokine level (field units, input)",
            "L": "lymph/antigen-transport level (field units, input)",
            "Sigma1": "a_11*T + a_12*D, emitted as sig_1 (dimensionless activation signal)",
            "Sigma2": "Sigma1 + a_21*V/(a_22+V), the T-activation drive (dimensionless)",
            "B_ei": "resistance-weighted infected-load CD8+ outflow term (per-MCS rate, input)",
            "G_ki": "resistance-weighted infected-load NK outflow term (per-MCS rate, input)",
            "b_*": "production/recruitment rate constants (per-MCS, patch-scaled)",
            "mu_*": "first-order decay/attrition rate constants (per-MCS)",
            "a_*": "Hill/Michaelis half-saturation constants (same units as their field)",
            "h_m, h_k, h_e, h_o": "Hill exponents (dimensionless)",
            "s_t, s_v, s_l": "time, population, and density cellularization scale factors (dimensionless)",
        },
        assumptions=[
            "Only the 10 systemic species are integrated; the spatialized cell/field species (H,I,M,K,E,L,C,V,DH) and the derived algebraic vars (D, Sigma1, Sigma2, DI, R) are held fixed across the step as inputs.",
            "The three nearby surrogates M_nb/K_nb/E_nb have d/dt = 0 in this process and are carried through unchanged (recruitment maintains them elsewhere).",
            "M/K/E inputs are the uncapped agent-level counts used as given; unlike the per-MCS CPM driver, no '+ state[M_nb]' nearby-surrogate correction is applied here.",
            "scipy solve_ivp LSODA (rtol=1e-6, atol=1e-9) approximates RoadRunner's stiff CVODE integrator; one default 60 s step advances exactly one MCS of the scaled per-MCS model.",
            "Integrated values are clamped at 0 to absorb tiny numerical negatives.",
            "The resistance R = F/(a_rf' + F) is computed in derived_inputs but does not enter the 10 integrated ODEs; the E input is likewise not consumed by the RHS in this build.",
            "Rate constants are pre-resolved and patch-cellularized (s_t, s_v, s_l) by price_ode.resolve_constants, with the T-equation denominator using unit-corrected g_2_T, d_2_T.",
            "CD8+ recruitment is driven by APC (P), not the chemokine field, and has no homeostatic inflow baseline (cd8_inflow = 0 when P <= 0) — an intentional source-faithful asymmetry, not to be symmetrized with macrophage/NK.",
        ],
        references=[
            "Mochan-Keef, Swigon, Ermentrout, Clermont (Price et al.) 2015 — global ODE model of the within-host influenza immune response.",
            "Sego et al. 2022 — CC3D ViralInfectionVTM spatial model and its spatial-coupling global-ODE string (ImmuneModelLib.py:immune_model_string); authority docs/cc3d-reference/sego2022-global-ode.md.",
        ],
    )

    config_schema = {
        # Task 3.2 fix: was "map[float]" -- bigraph_schema treats any
        # underscore-prefixed dict key as reserved schema metadata and
        # strips it (`bigraph_schema.strip_schema_keys` /
        # `is_schema_key`), which a "map[...]" config value goes through
        # when a real `Composite` realizes this process's config (unlike
        # `test_ode_process.py`'s direct `SystemicODEProcess(...)`
        # constructor call, which bypasses that pipeline entirely).
        # `price_ode.resolve_constants` stashes derived helpers under
        # underscore keys (`_s_t`, `_s_v`, `_s_l`) that `GlobalODE.__init__`
        # requires -- through a Composite those got silently dropped,
        # raising `KeyError: '_s_t'`. "tree" is not subject to the same
        # per-key stripping (verified: a "tree"-typed config value keeps
        # underscore keys intact through a full Composite realize), so it
        # is used here instead; no value semantics change, `consts` is
        # still consumed as a plain `dict[str, float]`.
        "consts": "tree",
        "num_epithelial": "integer",
        "dt_seconds": {"_type": "float", "_default": 60.0},
    }

    def initialize(self, config):
        self.consts = dict(config["consts"])
        self.num_epithelial = int(config["num_epithelial"])
        self.dt_seconds = float(config.get("dt_seconds", 60.0))
        self.ode = price_ode.GlobalODE(self.consts, num_epithelial=self.num_epithelial)
        # Same ICs `run_full_model` seeds before its coupling loop
        # (run.py:1284-1285 / :1833-1834): v0=0.0, resist0 defaults to 0.0.
        self.state = price_ode.initial_state(
            self.consts, num_epithelial=self.num_epithelial, v0=0.0
        )

    def inputs(self):
        return {
            "H": "integer",
            "I": "integer",
            "M": "integer",
            "K": "integer",
            "E": "integer",
            "DH": "integer",
            "V": "float",
            "F": "float",
            "C": "float",
            "L": "float",
            "B_ei": "float",
            "G_ki": "float",
        }

    def outputs(self):
        return {
            "ode_state": "overwrite[map[float]]",
            "recruit_drivers": "overwrite[map[float]]",
            "sig_1": "overwrite[float]",
        }

    def update(self, state, interval):
        state = state or {}

        inputs = {k: state.get(k, 0) for k in _ODE_INPUT_INT_KEYS}
        inputs.update({k: float(state.get(k, 0.0)) for k in _ODE_INPUT_FLOAT_KEYS})

        H = inputs["H"]
        I = inputs["I"]
        C = inputs["C"]
        G_ki = inputs["G_ki"]
        B_ei = inputs["B_ei"]

        self.state = self.ode.step(self.state, inputs, dt_seconds=self.dt_seconds)

        P = self.state["P"]
        recruit_drivers = {
            "macro_inflow": recruitment.macrophage_inflow(C, self.consts),
            "nk_inflow": recruitment.nk_inflow(C, self.consts),
            "cd8_inflow": recruitment.cd8_inflow(P, self.consts),
            "macro_outflow": recruitment.macrophage_outflow(self.consts),
            "nk_outflow": recruitment.nk_outflow(G_ki, self.consts),
            "cd8_outflow": recruitment.cd8_outflow(B_ei, self.consts),
        }

        D = self.num_epithelial - H - I
        sig_1 = self.consts["a_11"] * self.state["T"] + self.consts["a_12"] * D

        ode_state = {k: float(self.state[k]) for k in price_ode.INTEGRATED_STATES}

        return {
            "ode_state": ode_state,
            "recruit_drivers": recruit_drivers,
            "sig_1": float(sig_1),
        }
