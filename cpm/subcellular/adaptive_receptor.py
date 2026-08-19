from process_bigraph import Process
from cpm.subcellular.receptor import ReceptorSubcell


class AdaptiveReceptorSubcell(Process):
    """Receptor occupancy with slow adaptation (fold-change detection).

    theta = Hill occupancy of (ligand + background); m is a slow per-cell
    setpoint, dm/dt = epsilon*(theta - m); the response that drives chemotaxis
    is r = max(0, theta - m). epsilon=0 reduces this to pure Hill occupancy
    (response == theta) — the hill_occupancy rung and the adaptation knockout.
    A cell activates (naive->activated) when r >= activate_occupancy.
    """
    config_schema = {
        "kd": {"_type": "float", "_default": 2.9},
        "hill": {"_type": "float", "_default": 2.0},
        "conc_scale": {"_type": "float", "_default": 0.02},
        "activate_occupancy": {"_type": "float", "_default": 0.5},
        "epsilon": {"_type": "float", "_default": 0.1},
        "background": {"_type": "float", "_default": 0.0},
        "naive_type": {"_type": "integer", "_default": 2},
        "activated_type": {"_type": "integer", "_default": 3},
    }

    def __init__(self, config=None, core=None):
        # bigraph_schema's core.fill() treats an explicit 0.0/0 as "empty" for
        # Float/Integer leaves and silently substitutes the schema default (see
        # bigraph_schema.methods.is_empty). That would make epsilon=0.0 — the
        # documented pure-Hill knockout — unreachable. Stash the caller's raw
        # config so explicit zeros survive; fall back to the filled config
        # (which still supplies real defaults) for anything left unset.
        self._raw_config = dict(config) if config else {}
        super().__init__(config, core=core)

    def initialize(self, config):
        def pick(key):
            return self._raw_config[key] if key in self._raw_config else config[key]

        kd, hill, conc_scale = pick("kd"), pick("hill"), pick("conc_scale")
        activate_occupancy = pick("activate_occupancy")
        naive_type, activated_type = pick("naive_type"), pick("activated_type")

        # reuse the canonical Hill formula (no fourth copy of c^h/(kd^h+c^h))
        self._hill = ReceptorSubcell(
            {"kd": kd, "hill": hill, "conc_scale": conc_scale,
             "activate_occupancy": activate_occupancy,
             "naive_type": naive_type, "activated_type": activated_type},
            core=self.core)
        self.epsilon = float(pick("epsilon")); self.background = float(pick("background"))
        self.activate_occupancy = float(activate_occupancy)
        self.naive_type = int(naive_type); self.activated_type = int(activated_type)

    def inputs(self):  return {"ligand": "float", "m_prev": "float"}
    def outputs(self): return {"fate": "overwrite[integer]", "m": "overwrite[float]"}

    def signal(self, ligand, m_prev):
        theta = self._hill.occupancy(float(ligand) + self.background)
        m_new = float(m_prev) + self.epsilon * (theta - float(m_prev))
        response = max(0.0, theta - float(m_prev))
        return theta, m_new, response

    def update(self, state, interval):
        s = state or {}
        _, m_new, response = self.signal(s.get("ligand", 0.0), s.get("m_prev", 0.0))
        fate = self.activated_type if response >= self.activate_occupancy else self.naive_type
        return {"fate": fate, "m": m_new}
