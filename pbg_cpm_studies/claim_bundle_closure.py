"""Derive the claim-bundle closure verdict for the chemotactic-recruitment slice.

This computes the four scientific-closure predicates from the biological-claim
bundle (``claim.yaml`` / ``mechanism.yaml`` / ``realization.yaml``) plus its
member studies, and writes ``closure.generated.yaml`` next to the bundle. The
verdict is DERIVED, never hand-set — mirroring the study gate verdict's
``evaluated_by: code`` discipline (see ``gg1993/validate.py``). It is the
concrete, in-repo realization of the Vivarium Semantic Modeling Layer's closure
card for this investigation's reference slice.

Predicates
----------
well_formed
    The three bundle files are schema-shaped and their cross-references resolve
    (realization -> claim/mechanism ids; mechanism.linked_claims -> claim id).
claim_complete
    No unresolved hole is ``required_for`` a downstream claim. The
    ``hole:lambda-calibration`` hole is required_for a quantitative prediction,
    so this slice is deliberately ``false`` — demonstrated qualitatively, not
    calibrated.
backend_realizable
    The realization names a backend, binds the mechanism roles, and every
    ``generated-studies`` slug exists on disk.
execution_ground
    Every generated study has a completed, gated run (``gate_status: passed``).
    (The workbench computes the finer "every readout is a real emit-plan leaf"
    check; this is the bundle-level approximation.)

Pure stdlib + PyYAML — no Rust engine import, so it stays fast and decoupled.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import yaml

# repo layout: <repo>/pbg_cpm_studies/claim_bundle_closure.py ; workspace is a sibling
_REPO_ROOT = Path(__file__).resolve().parents[1]
_DEFAULT_BUNDLE = (
    _REPO_ROOT
    / "workspace/investigations/chemotactic-recruitment/claim-bundle"
)
_DEFAULT_STUDIES = _REPO_ROOT / "workspace/studies"


def _load(path: Path) -> dict[str, Any]:
    with path.open() as fh:
        data = yaml.safe_load(fh)
    return data if isinstance(data, dict) else {}


def _slug(ref: str) -> str:
    """`study:recruitment-baseline` -> `recruitment-baseline`."""
    return ref.split(":", 1)[-1].strip()


def derive_closure(
    bundle_dir: Path | str = _DEFAULT_BUNDLE,
    studies_dir: Path | str = _DEFAULT_STUDIES,
) -> dict[str, Any]:
    """Compute the closure verdict for a claim bundle. Returns the closure dict."""
    bundle_dir = Path(bundle_dir)
    studies_dir = Path(studies_dir)

    claim = _load(bundle_dir / "claim.yaml")
    mechanism = _load(bundle_dir / "mechanism.yaml")
    realization = _load(bundle_dir / "realization.yaml")

    problems: list[str] = []

    # --- well_formed -------------------------------------------------------
    if not claim.get("id"):
        problems.append("claim.yaml missing id")
    if not isinstance(claim.get("assertion"), dict):
        problems.append("claim.yaml missing assertion")
    if not mechanism.get("id"):
        problems.append("mechanism.yaml missing id")
    if not isinstance(mechanism.get("roles"), dict):
        problems.append("mechanism.yaml missing roles")
    if realization.get("claim") != claim.get("id"):
        problems.append("realization.claim does not resolve to claim.id")
    if realization.get("mechanism") != mechanism.get("id"):
        problems.append("realization.mechanism does not resolve to mechanism.id")
    if claim.get("id") not in (mechanism.get("linked_claims") or []):
        problems.append("mechanism.linked_claims does not reference the claim")
    well_formed = not problems

    # --- claim_complete ----------------------------------------------------
    holes = list(realization.get("unresolved") or [])
    required_holes = [h for h in holes if h.get("required_for")]
    claim_complete = len(required_holes) == 0

    # --- backend_realizable ------------------------------------------------
    gen_slugs = [_slug(s) for s in (realization.get("generated-studies") or [])]
    missing_studies = [
        s for s in gen_slugs if not (studies_dir / s / "study.yaml").is_file()
    ]
    backend_realizable = bool(
        realization.get("backend")
        and realization.get("bindings")
        and gen_slugs
        and not missing_studies
    )

    # --- execution_ground --------------------------------------------------
    ungated: list[str] = []
    for s in gen_slugs:
        spec = _load(studies_dir / s / "study.yaml")
        if str(spec.get("gate_status") or "").strip().lower() != "passed":
            ungated.append(s)
    execution_ground = bool(gen_slugs) and not ungated

    closure: dict[str, Any] = {
        "claim": claim.get("id"),
        "well_formed": well_formed,
        "claim_complete": claim_complete,
        "backend_realizable": backend_realizable,
        "execution_ground": execution_ground,
        "evaluated_by": "code",
        "unresolved": [
            {
                "id": h.get("id"),
                "required_for": h.get("required_for") or [],
                "consequence": h.get("consequence-if-unresolved"),
            }
            for h in required_holes
        ],
    }
    if problems:
        closure["problems"] = problems
    if missing_studies:
        closure["missing_studies"] = missing_studies
    if ungated:
        closure["ungated_studies"] = ungated
    return closure


def write_closure(
    bundle_dir: Path | str = _DEFAULT_BUNDLE,
    studies_dir: Path | str = _DEFAULT_STUDIES,
) -> Path:
    """Derive the closure and write ``closure.generated.yaml`` into the bundle."""
    bundle_dir = Path(bundle_dir)
    closure = derive_closure(bundle_dir, studies_dir)
    out = bundle_dir / "closure.generated.yaml"
    header = (
        "# GENERATED by pbg_cpm_studies.claim_bundle_closure — do not hand-edit.\n"
        "# Derived closure verdict for the claim bundle (evaluated_by: code).\n"
        "# Regenerate: python pbg_cpm_studies/claim_bundle_closure.py\n"
    )
    with out.open("w") as fh:
        fh.write(header)
        yaml.safe_dump({"closure": closure}, fh, sort_keys=False, allow_unicode=True)
    return out


def main() -> int:
    out = write_closure()
    closure = derive_closure()
    print(f"wrote {out.relative_to(_REPO_ROOT)}")
    for k in ("well_formed", "claim_complete", "backend_realizable", "execution_ground"):
        print(f"  {k}: {closure[k]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
