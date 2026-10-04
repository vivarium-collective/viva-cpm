from process_bigraph import Process

from bigraph_schema.contract import ProcessContract

from viva_cpm.schema import load_world


class CPMProcess(Process):
    """process-bigraph wrapper around the Rust CPM engine, with a per-cell
    coupling surface so subcellular processes can read the local environment
    and drive differentiation via a cell-type switch."""

    contract = ProcessContract(
        summary=(
            "Advances a Cellular Potts (Glazier–Graner–Hogeweg) model by running "
            "`mcs_per_update` Monte Carlo sweeps of the wrapped Rust engine each process "
            "step. Each sweep attempts random boundary copies between adjacent lattice "
            "sites and accepts or rejects them by the Metropolis rule on the change ΔH in a "
            "Hamiltonian that sums cell–cell adhesion, a per-cell volume constraint, and a "
            "per-cell surface-area constraint (plus any optional terms the spec opts in). "
            "Before stepping it applies incoming cell-type switches (`fates`), and after "
            "stepping it emits per-cell geometry (volumes, types, centre-of-mass positions) "
            "and local-environment readouts (mean field value at each cell, count of "
            "secretory neighbours) so subcellular processes can couple to the tissue."
        ),
        description=(
            "Core cellular-dynamics engine of the workspace. The morphology lives in a Rust "
            "`World` built from a declarative CPM `spec` (see viva_cpm.schema.load_world); this "
            "Process owns the step loop and the coupling surface to process-bigraph.\n\n"
            "Energy terms compose additively in the Metropolis Hamiltonian; each term's "
            "`*_type` list opts specific cell types in, so tissues mix behaviours. The core "
            "terms are:\n"
            "  - Adhesion/contact: sum over neighbouring site pairs of different cells of the "
            "contact energy J(tau_a, tau_b) between their types (the `contact` adhesion matrix).\n"
            "  - Volume constraint: lambda_V * (V - V_target)^2 per cell.\n"
            "  - Surface constraint: lambda_S * (S - S_target)^2 per cell.\n"
            "Optional spec terms (chemotaxis on reaction-diffusion fields, connectivity, "
            "basement membrane, junctions, elongation length, external force) add further "
            "contributions to ΔH but are not computed in this wrapper.\n\n"
            "Each update: (1) write incoming fates into the engine as cell-type switches; "
            "(2) run `mcs_per_update` Monte Carlo sweeps (each sweep also advances the "
            "reaction-diffusion fields when present); (3) read back per-cell geometry and "
            "the per-cell environment readouts used for subcellular coupling."
        ),
        inputs={
            "fates": (
                "map[integer] keyed by string cell id -> new cell-type tag. Lets a subcellular "
                "process drive differentiation by wiring a single cell's fate at "
                "`[fates, str(cid)]`; each non-zero value for a positive cell id is applied via "
                "world.set_cell_type before the sweeps. The composite pre-seeds a key per wired "
                "cell so per-key overwrites land (an absent map key would be dropped)."
            ),
        },
        outputs={
            "volumes": (
                "overwrite[list] of current cell volumes V (lattice-site counts), indexed by "
                "cell id (index 0 is the medium); emitted as a whole list each step."
            ),
            "types": (
                "overwrite[list] of each cell's current type tag tau, indexed by cell id "
                "(index 0 is the medium)."
            ),
            "positions": (
                "overwrite[list] of each cell's centre-of-mass [x, y, z] in lattice units, "
                "indexed by cell id."
            ),
            "field_at_cell": (
                "overwrite[map[float]] keyed by string cell id -> mean value of field 0 over "
                "the cell's occupied sites (empty when n_fields == 0). The overwrite wrapper "
                "replaces the whole store atomically so a subcell can read `[field_at_cell, str(cid)]`."
            ),
            "neighbor_secretory": (
                "overwrite[map[integer]] keyed by string cell id -> number of distinct "
                "face-adjacent neighbour cells whose type is in `secretory_types` (medium and "
                "same-cell contacts excluded)."
            ),
        },
        config={
            "spec": (
                "Declarative CPM tree passed to load_world: potts block (dims, boundary, "
                "neighbor_order, temperature T, seed), cell/seed-label placement with per-cell "
                "target_volume/lambda_volume/target_surface/lambda_surface, the contact adhesion "
                "matrix J, optional reaction-diffusion fields, and optional extension terms."
            ),
            "mcs_per_update": (
                "integer (default 10): number of Monte Carlo sweeps (MCS) the engine runs per "
                "process update — the number of simulated CPM time units advanced each step."
            ),
            "n_fields": (
                "integer (default 0): number of reaction-diffusion fields; when > 0 the mean of "
                "field 0 at each cell is emitted on field_at_cell."
            ),
            "secretory_types": (
                "list (default []) of cell-type tags treated as secretory when tallying "
                "neighbor_secretory face-adjacency counts."
            ),
        },
        math=[
            "H = \\sum_{(i,j),\\,\\sigma_i \\neq \\sigma_j} J(\\tau_{\\sigma_i}, \\tau_{\\sigma_j})"
            " \\;+\\; \\sum_{\\sigma} \\lambda_V (V_\\sigma - V^{*}_\\sigma)^2"
            " \\;+\\; \\sum_{\\sigma} \\lambda_S (S_\\sigma - S^{*}_\\sigma)^2",
            "\\Delta H = H_{\\text{after copy}} - H_{\\text{before copy}}",
            "P(\\text{accept}) = \\begin{cases} 1 & \\Delta H \\le 0 \\\\ \\exp(-\\Delta H / T) & \\Delta H > 0 \\end{cases}",
            "1\\ \\text{MCS} = N\\ \\text{copy attempts},\\quad N = \\text{number of lattice sites}",
        ],
        symbols={
            "sigma": "cell index (spin) labelling the cell occupying a lattice site (dimensionless; 0 = medium)",
            "tau_sigma": "cell-type tag of cell sigma, selecting adhesion and constraint parameters (dimensionless)",
            "J": "contact (adhesion) energy between two abutting cell types from the contact matrix (energy units)",
            "H": "total CPM Hamiltonian, sum of adhesion, volume, and surface energies (energy units)",
            "ΔH": "change in H produced by a single trial site-copy, used in the Metropolis rule (energy units)",
            "V_sigma": "current volume of cell sigma (number of lattice sites, voxels)",
            "V_target": "target volume V* of a cell set by its type (voxels)",
            "lambda_V": "volume-constraint stiffness, penalty per squared volume deviation (energy / voxel^2)",
            "S_sigma": "current surface area of cell sigma (number of boundary faces, sites)",
            "S_target": "target surface area S* of a cell set by its type (sites)",
            "lambda_S": "surface-constraint stiffness, penalty per squared surface deviation (energy / site^2)",
            "T": "CPM temperature / fluctuation amplitude in the Metropolis acceptance exp(-ΔH/T) (energy units)",
            "MCS": "Monte Carlo sweep: N site-copy attempts where N is the number of lattice sites (time unit)",
            "mcs_per_update": "Monte Carlo sweeps advanced per process step (sweeps)",
        },
        assumptions=[
            "Dynamics (Hamiltonian evaluation, trial copies, Metropolis acceptance, field PDE "
            "integration) live in the wrapped Rust engine; this Process only drives the step "
            "loop and the port coupling.",
            "Energy terms are additive and each term's `*_type` list opts specific cell types "
            "in, so only the core adhesion/volume/surface terms are always present.",
            "Cell ids are contiguous positive integers with 0 reserved for the medium; per-cell "
            "output lists are indexed by id and the id-0 medium slot is excluded from the "
            "per-cell maps.",
            "fates are applied as hard type switches before stepping; a value of 0 or a "
            "non-positive cell id is ignored.",
            "neighbor_secretory counts only face-adjacent (6-neighbour) distinct cell pairs; "
            "medium and same-cell contacts do not count.",
            "field_at_cell reports only field index 0 and is empty unless n_fields > 0.",
        ],
        references=[
            "Graner & Glazier (1992), Simulation of biological cell sorting using a two-dimensional "
            "extended Potts model, Phys. Rev. Lett. 69:2013.",
            "Glazier & Graner (1993), Simulation of the differential adhesion driven rearrangement "
            "of biological cells, Phys. Rev. E 47:2128.",
            "Hogeweg (2000), Evolving mechanisms of morphogenesis: on the interplay between "
            "differential adhesion and cell differentiation, J. Theor. Biol. 203:317.",
        ],
    )

    config_schema = {
        "spec": "tree",
        "mcs_per_update": {"_type": "integer", "_default": 10},
        "n_fields": {"_type": "integer", "_default": 0},
        "secretory_types": {"_type": "list", "_default": []},
    }

    def initialize(self, config):
        self.world = load_world(self.config["spec"])
        self.mcs = int(self.config["mcs_per_update"])
        self.n_fields = int(self.config["n_fields"])
        self.secretory = set(int(t) for t in self.config.get("secretory_types") or [])
        self.dims = self.world.dims()

    def inputs(self):
        # fates is a map keyed by string cell id so subcellular processes can
        # wire a single cell's fate (``[fates, str(cid)]``); integer-index paths
        # into a plain list are rejected by the type system. The composite
        # pre-initialises the fates store with a key per wired cell so per-key
        # overwrite writes land (an absent map key would be dropped).
        return {"fates": "map[integer]"}

    def outputs(self):
        return {
            "volumes": "overwrite[list]",
            "types": "overwrite[list]",
            "positions": "overwrite[list]",
            # per-cell readouts are maps keyed by string cell id, so a subcell
            # can wire ``[<port>, str(cid)]`` to its own cell. ``overwrite`` wraps
            # the map so emitting the whole dict atomically replaces the store
            # (a bare ``map[..]`` apply ignores keys absent from the prior map).
            "field_at_cell": "overwrite[map[float]]",
            "neighbor_secretory": "overwrite[map[integer]]",
        }

    def _neighbor_secretory_counts(self, types):
        """Count, per cell, face-adjacent cells whose type is secretory."""
        nx, ny, nz = self.dims
        lab = self.world.snapshot()
        n = len(types)
        counts = [0] * n
        seen = [set() for _ in range(n)]
        def owner(x, y, z):
            return lab[x + y * nx + z * nx * ny]
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    a = owner(x, y, z)
                    if a == 0:
                        continue
                    for dx, dy, dz in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                        xx, yy, zz = x + dx, y + dy, z + dz
                        if xx >= nx or yy >= ny or zz >= nz:
                            continue
                        b = owner(xx, yy, zz)
                        if b == a or b == 0:
                            continue
                        if types[b] in self.secretory and b not in seen[a]:
                            counts[a] += 1; seen[a].add(b)
                        if types[a] in self.secretory and a not in seen[b]:
                            counts[b] += 1; seen[b].add(a)
        return counts

    def update(self, state, interval):
        fates = (state or {}).get("fates") or {}
        for cid_key, t in fates.items():
            cid = int(cid_key)
            if t and cid > 0:
                self.world.set_cell_type(cid, int(t))
        self.world.step(self.mcs)

        types = list(self.world.cell_types())
        n = len(types)
        field_at = {}
        if self.n_fields > 0:
            for cid in range(1, n):
                field_at[str(cid)] = self.world.field_mean_at_cell(0, cid)
        neigh = self._neighbor_secretory_counts(types)
        return {
            "volumes": list(self.world.cell_volumes()),
            "types": types,
            "positions": [list(c) for c in self.world.cell_coms()],
            "field_at_cell": field_at,
            "neighbor_secretory": {str(cid): neigh[cid] for cid in range(1, n)},
        }
