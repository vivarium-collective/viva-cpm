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

from . import price_ode, recruitment

_ODE_INPUT_INT_KEYS = ("H", "I", "M", "K", "E", "DH")
_ODE_INPUT_FLOAT_KEYS = ("V", "F", "C", "L", "B_ei", "G_ki")


class SystemicODEProcess(Process):
    config_schema = {
        "consts": "map[float]",
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
