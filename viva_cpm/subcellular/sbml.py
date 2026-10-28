from process_bigraph import Process
from viva_tellurium.processes import TelluriumProcess


class SBMLSubcell(Process):
    """Per-cell SBML ODE via viva-tellurium. The environmental ligand is pushed
    into a held floating species each step; a chosen species is published as
    the cell's scalar ``state`` (e.g. stemness)."""

    config_schema = {
        "model": {"_type": "string", "_default": ""},
        "ligand_species": {"_type": "string", "_default": "Wnt"},
        "state_species": {"_type": "string", "_default": "S"},
        "ligand_scale": {"_type": "float", "_default": 1.0},
    }

    def initialize(self, config):
        self.ligand_species = config["ligand_species"]
        self.state_species = config["state_species"]
        self.scale = float(config["ligand_scale"])
        # Pin the legacy per-species map port. viva-tellurium >= 0.1.2 defaults
        # bulk_species=True, which makes the `species` port a numpy array AND expects
        # the species input as an array — so this process's map-shaped ligand override
        # and string-keyed state read both break (IndexError, then a dropped override).
        # bulk_species=False is the bit-identical back-compat map port. Older
        # viva-tellurium (no bulk_species key) ignores the extra config harmlessly.
        self._tp = TelluriumProcess(
            {"model": config["model"], "bulk_species": False}, core=self.core)

    def inputs(self):
        return {"ligand": "float"}

    def outputs(self):
        return {"state": "overwrite[float]"}

    def update(self, state, interval):
        ligand = float((state or {}).get("ligand", 0.0)) * self.scale
        out = self._tp.update({"species": {self.ligand_species: ligand}}, interval)
        species = out["species"]
        # viva-tellurium >= 0.1.2 bulk mode returns `species` as a numpy array ordered by
        # getFloatingSpeciesIds; 0.1.0 returned a {id: value} map. Read the chosen state
        # species by key when it's a map, else by its position in the id order.
        if isinstance(species, dict):
            value = species[self.state_species]
        else:
            value = species[self._tp._species_ids.index(self.state_species)]
        return {"state": float(value)}
