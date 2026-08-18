#!/usr/bin/env python
"""Render the receptor-recruitment study figures to each study's ``viz/`` dir.

The three receptor visualizations in
``pbg_cpm_studies.visualizations.receptor_studies`` are self-contained (their
Plotly data is baked in from the analysed 5-seed runs, so rendering needs no
live simulation). This script wraps each fragment in a standalone HTML page —
matching the ``chemotaxis/scene.py`` convention that produced
``recruitment-baseline/viz/scene.html`` — and writes them into both receptor
studies so they carry rendered output (study-audit L3 outputs-present).

Run: ``python scripts/render_receptor_viz.py``
"""
from __future__ import annotations

from pathlib import Path

from pbg_cpm_studies.visualizations import receptor_studies as V

# study slug -> the receptor investigation's shared figure set (baseline vs
# blocked are two conditions of one comparison, so both carry the full set;
# the recruitment curve overlays both arms with CI ribbons).
STUDIES = ("recruitment-receptor-baseline", "recruitment-receptor-blocked")

FIGURES = (
    ("scene.html", V.ReceptorRecruitmentScene),
    ("occupancy-law.html", V.ReceptorOccupancyLaw),
    ("recruitment-curve.html", V.ReceptorRecruitmentCurve),
)


def _standalone(fragment: str, title: str) -> str:
    """Wrap a Plotly card fragment in a minimal dark standalone document."""
    return (
        '<!doctype html><html><head><meta charset="utf-8">'
        f"<title>{title}</title>"
        '<style>body{margin:0;background:#0a1020;padding:16px 0}</style>'
        f"</head><body>{fragment}</body></html>"
    )


def main() -> int:
    ws = Path(__file__).resolve().parent.parent
    for slug in STUDIES:
        viz = ws / "workspace" / "studies" / slug / "viz"
        viz.mkdir(parents=True, exist_ok=True)
        for filename, fn in FIGURES:
            title = f"{slug} — {fn.__name__}"
            (viz / filename).write_text(_standalone(fn(), title), encoding="utf-8")
            print(f"wrote {viz.relative_to(ws)}/{filename}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
