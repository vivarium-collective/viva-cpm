# Sego et al. 2022 — CC3D source acquisition notes

Paper: T.J. Sego, Ericka D. Mochan, G. Bard Ermentrout, James A. Glazier, "A multiscale
multicellular spatiotemporal model of local influenza infection and immune response,"
*J. Theor. Biol.* 532 (2022) 110918. DOI `10.1016/j.jtbi.2021.110918`. PMID 34592264,
PMCID PMC8478073. Local copy: `workspace/references/papers/sego-2022-influenza.pdf`.

The paper states (Sec. 2, "Models and methods"): "The package to run the model in CompuCell3D
can be found in Supplementary Materials, in the package Source Code." This is Appendix A /
the article's Supplementary Materials, not a separate GitHub repo maintained under
`github.com/tjsego` — see "What was tried" below.

## Result: RETRIEVED — full CC3D source with exact per-parameter constants

The supplementary "Source Code" package (`mmc3.zip`, 20.7 KB) was downloaded and unpacked.
It contains the complete, exact CC3D project used to generate the paper's results: CC3DML,
all Python steppables, and the ODE-model (Antimony/SBML) generator functions. Every constant
in the model inputs files is a literal numeric value with an inline description string — no
transcription or approximation from the paper's prose/tables is needed for anything covered
by this package. Per the `docs/cc3d-reference/demo-parameters.md` convention, this source
package (not just validated behavior) is the primary authority for Increment 0 parameter
transcription; the paper's Tables 1-4 remain the authority for the ODE-level formulas,
scenario definitions, and any constant this package doesn't set explicitly.

### Exact URLs that worked

- Elsevier supplementary-content CDN (direct file, no auth, no JS/PoW challenge):
  `https://ars.els-cdn.com/content/image/1-s2.0-S0022519321003374-mmc3.zip`
  (200 OK, `content-length: 21177`, matches ScienceDirect's listed "20.7 KB"; the other two
  supplementary files are at the same path pattern: `...-mmc1.docx`, `...-mmc2.docx`.)
- Europe PMC supplementary-files bundle (all three `mmc*` files zipped together):
  `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8478073/supplementaryFiles`
  (200 OK, returns a ~1 MB zip containing mmc1.docx, mmc2.docx, mmc3.zip).
- Canonical article page (for reference, not machine-fetchable — see below):
  `https://doi.org/10.1016/j.jtbi.2021.110918` → redirects to
  `https://linkinghub.elsevier.com/retrieve/pii/S0022519321003374` →
  `https://www.sciencedirect.com/science/article/pii/S0022519321003374`.
- PMC article (open access, PMCID PMC8478073):
  `https://pmc.ncbi.nlm.nih.gov/articles/PMC8478073/`, with supplementary files linked at
  `https://pmc.ncbi.nlm.nih.gov/articles/instance/8478073/bin/mmc{1,2,3}.{docx,docx,zip}`.

### What did NOT work directly

- `WebFetch` on the ScienceDirect article page itself returned HTTP 403 (bot-blocked).
- The PMC `bin/mmc3.zip` links return HTTP 200 but serve an HTML "Preparing to download..."
  interstitial page (a client-side proof-of-work challenge, `cloudpmc-viewer-pow`) instead of
  the file — this requires executing JS in a real browser, so `curl`/`WebFetch` cannot follow
  it. **Do not rely on the PMC `bin/` URL for automated re-fetching; use the Elsevier CDN URL
  or the Europe PMC bundle above instead**, both of which are plain HTTP GETs that returned
  the real file with a normal browser `User-Agent` header.
- `mmc3.zip` is itself a thin ZIP wrapper containing a single 7z archive, `Source Code.7z`
  (21,041 bytes) — needs a 7z extractor (macOS/BSD `unzip`/`tar` can't read `.7z`; `py7zr`
  via `pip install py7zr` handled it fine).

### Package contents (unpacked `Source Code.7z`)

```
Source Code/
├── ViralInfectionVTM.cc3d                       # CC3D project file (CC3D version 4.2.2)
├── Simulation/
│   ├── ViralInfectionVTM.xml                    # CC3DML: lattice, Potts, cell types,
│   │                                             #   Contact (adhesion J), diffusion fields
│   ├── ViralInfectionVTM.py                      # CC3D run/steppable-registration script
│   ├── ViralInfectionVTMSteppables.py            # 15 steppable classes (1637 lines): init,
│   │                                             #   viral internalization/secretion/death,
│   │                                             #   contact killing, chemotaxis, IFN/
│   │                                             #   chemokine/IL-10 secretion, recovery, etc.
│   ├── ViralInfectionVTMModelInputs.py           # ALL numeric constants for the spatial
│   │                                             #   model, each with a __param_desc__ string
│   │                                             #   (unit conversions, lattice/volume
│   │                                             #   constants, MCS-to-time mapping, etc.)
│   └── ViralInfectionVTMLib.py                   # shared helper functions
├── ImmuneModel/
│   ├── ImmuneModelInputs.py                      # constants for the immune-model layer
│   │                                             #   (decay rates carried over verbatim from
│   │                                             #   the calibrated ODE model, local immune
│   │                                             #   population ratios, initial viral load
│   │                                             #   options, etc.), imports the Simulation
│   │                                             #   inputs above as a base
│   ├── ImmuneModelLib.py                         # `immune_model_string(...)` /
│   │                                             #   `immune_model_string_ode(...)`: generate
│   │                                             #   the ODE sub-model as an **Antimony/SBML**
│   │                                             #   string (run via CC3D's SBML solver) —
│   │                                             #   this is the calibrated Price et al. 2015
│   │                                             #   ODE model, cellularized
│   └── ImmuneModelSteppableBasePy.py
└── nCoVToolkit/
    ├── nCoVSteppableBase.py
    └── nCoVUtils.py                              # small utilities (e.g. rate/probability
                                                    #   conversions)
```

- CC3DML header comments attribute the base simulation structure to an earlier, related
  paper: Sego et al., "A modular framework for multiscale multicellular spatial modeling of
  viral infection, immune response and drug therapy timing and efficacy in epithelial
  tissues," bioRxiv 2020.04.27.064139 (published as Sego et al. 2020, *PLOS Comput. Biol.*
  16(12):e1008451) — the Sego 2022 influenza package is a parameterization/extension of that
  framework, not a from-scratch model. Useful for cross-checking steppable logic if a
  constant or mechanism is ambiguous.
- Confirmed examples of exact retrievable constants (used directly, no transcription
  approximation): `s_to_mcs = 60` (1 min/MCS, matches task 0.2's Δt = 1 min/step),
  `um_to_lat_width = 2.0` (2 µm/lattice site, matches task 0.2's 2 µm lattice), Potts
  `Temperature = 10.0`, `Dimensions x=175 y=175 z=2` (default demo lattice; task 0.2's 0.3/1.0
  mm patches are separate scenario configs, likely set via `um_to_lat_width`/dimensions
  overrides — check `ViralInfectionVTMModelInputs.py`/`ImmuneModelInputs.py` directly for the
  scenario-specific values when transcribing), full `Contact` energy matrix in
  `ViralInfectionVTM.xml`, and virus/type-I-IFN/IL-10 decay constants in
  `ImmuneModelInputs.py` (e.g. `virus_decay_ODE = 0.412015488642712`,
  `t1ifn_decay_ODE = 112.230629642229`) taken verbatim from the calibrated ODE model.

### Retrieval method used (for reproducing / re-fetching)

```bash
curl -sL -A "Mozilla/5.0" \
  "https://ars.els-cdn.com/content/image/1-s2.0-S0022519321003374-mmc3.zip" \
  -o mmc3.zip
unzip mmc3.zip                    # -> "Source Code.7z"
pip install py7zr                 # macOS/BSD unzip/tar can't read .7z
python3 -c "import py7zr; py7zr.SevenZipFile('Source Code.7z').extractall('unpacked')"
```

`shasum -a 256` of the two intermediate archives, for reference (not committed — see below):
- `mmc3.zip` (Elsevier CDN download): `6a17288df4dd1b1267d3272f3dc19c676b3c0dc06d52247a2b11894208f4760a`
- `Source Code.7z` (unwrapped): `ddfc0131c6c430f2d079f8c7d2b0bf1650ce5a9bc445eb83d51d1ae79a44ac56`

### Not vendored into this repo

The extracted source tree was inspected in a scratch directory only and is **not** committed
here — this task's deliverable is the acquisition record, and the package is reliably
re-fetchable via the Elsevier CDN URL above. Task 0.2 (parameter transcription) and later
tasks that need to read exact constants from these files should re-run the retrieval command
above (or ask for the source tree to be vendored under `docs/cc3d-reference/` at that point,
if repeated re-fetching becomes a friction point).

## What was tried for (b) GitHub/nanoHUB and (c) Price et al. 2015

- Searched `github.com/tjsego` and "Sego cellularization influenza CompuCell3D": no
  dedicated GitHub repo for this specific paper's model was found. T.J. Sego's public repos
  found (e.g. `tjsego/simservice`) are unrelated infrastructure libraries, not this model.
  The general `CompuCell3D/CompuCell3D` repo is the simulator engine itself, not this
  paper's parameterized project. Since (a) fully succeeded, (b)/(c) were not needed as a
  primary source, but nanoHUB was not separately searched given (a)'s success — flagging
  in case a nanoHUB-hosted interactive version is independently useful later (e.g. for
  behavioral cross-validation) but it is not required for parameter retrieval.
- Price et al. 2015 (*J. Theor. Biol.* 374:83-93, the underlying non-spatial ODE model) was
  not separately fetched, since the Sego 2022 source package's `ImmuneModelInputs.py`/
  `ImmuneModelLib.py` already carry the calibrated ODE constants and the Antimony/SBML
  model string verbatim (comment: "Taken from calibrated ODE model parameters"). If a
  constant is needed that isn't in the CC3D package, fall back to Price et al. 2015 directly.

## Fidelity convention for downstream tasks

Because the exact source is retrievable, prefer literal values from
`ViralInfectionVTMModelInputs.py` / `ImmuneModelInputs.py` / `ViralInfectionVTM.xml` over the
paper's Tables 1-4 wherever they overlap (the source is the ground truth the paper's tables
summarize). Where the source and paper disagree, note the discrepancy explicitly rather than
silently preferring one. For anything genuinely absent from both the source and the paper,
follow the `demo-parameters.md` convention: mark it "CC3D-typical" and make the *validated
behavior* (not the exact constant) the fidelity criterion.
