#!/usr/bin/env python
"""Render the Increment 3-8 influenza study figures to each study's ``viz/`` dir.

The six backfilled influenza visualizations in
``pbg_cpm_studies.visualizations.influenza_studies`` are self-contained (their
Plotly data is baked in from real engine runs at each study's documented
scenario, so rendering needs no live simulation). This script wraps each
fragment in a standalone HTML page — matching the convention that produced
``epithelial-sheet-baseline/viz/confluent-sheet.html`` — and writes them into
each study so the dashboard's ``embed_visualizations`` URL resolves.

Run: ``python scripts/render_influenza_studies_viz.py``
"""
from __future__ import annotations

from pathlib import Path

from pbg_cpm_studies.visualizations import influenza_studies as V

# study slug -> (embed filename, zero-arg accessor)
STUDIES = {
    "ifn-resistance": ("ifn-resistance.html", V.InfluenzaIfnResistance),
    "epithelial-fate": ("epithelial-fate.html", V.InfluenzaEpithelialFate),
    "macrophage-response": ("macrophage-response.html", V.InfluenzaMacrophageResponse),
    "signaling-fields": ("signaling-fields.html", V.InfluenzaSignalingFields),
    "cytotoxic-killing": ("cytotoxic-killing.html", V.InfluenzaCytotoxicKilling),
    "global-coupling": ("global-coupling.html", V.InfluenzaGlobalCoupling),
}


def _standalone(fragment: str, title: str) -> str:
    """Wrap a Plotly card fragment in a minimal light standalone document
    (same shell as the other influenza study embeds)."""
    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f"<title>{title}</title></head>"
        '<body style="margin:0;padding:20px;background:#f6f6f4;">'
        f"{fragment}</body></html>"
    )


def main() -> int:
    ws = Path(__file__).resolve().parent.parent
    for slug, (filename, fn) in STUDIES.items():
        viz = ws / "workspace" / "studies" / slug / "viz"
        viz.mkdir(parents=True, exist_ok=True)
        title = filename[: -len(".html")]
        (viz / filename).write_text(_standalone(fn(), title), encoding="utf-8")
        print(f"wrote {viz.relative_to(ws)}/{filename}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
