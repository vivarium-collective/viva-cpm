from process_bigraph import Process


class ReceptorSubcell(Process):
    """Minimal per-cell ligand-receptor occupancy model. Reads the local
    chemokine concentration, computes fractional receptor occupancy via a Hill
    relation, and emits a responder cell-type (fate): activated when occupancy
    crosses threshold, else naive. Deterministic.

        theta = c^h / (Kd^h + c^h),   c = ligand * conc_scale

    With hill=1, conc_scale=1 and activate_occupancy=0.5 the activation threshold
    sits exactly at c = Kd, so the cited Kd is directly visible in a dose-response.
    """

    config_schema = {
        "kd": {"_type": "float", "_default": 1.0},
        "hill": {"_type": "float", "_default": 1.0},
        "conc_scale": {"_type": "float", "_default": 1.0},
        "activate_occupancy": {"_type": "float", "_default": 0.5},
        "naive_type": {"_type": "integer", "_default": 2},
        "activated_type": {"_type": "integer", "_default": 3},
    }

    def initialize(self, config):
        self.kd = float(config["kd"])
        self.hill = float(config["hill"])
        self.conc_scale = float(config["conc_scale"])
        self.activate_occupancy = float(config["activate_occupancy"])
        self.naive_type = int(config["naive_type"])
        self.activated_type = int(config["activated_type"])

    def inputs(self):
        return {"ligand": "float"}

    def outputs(self):
        return {"fate": "overwrite[integer]"}

    def occupancy(self, ligand):
        c = max(0.0, float(ligand)) * self.conc_scale
        if c <= 0.0:
            return 0.0
        ch = c ** self.hill
        return ch / (self.kd ** self.hill + ch)

    def update(self, state, interval):
        theta = self.occupancy(float((state or {}).get("ligand", 0.0)))
        fate = self.activated_type if theta >= self.activate_occupancy else self.naive_type
        return {"fate": fate}
