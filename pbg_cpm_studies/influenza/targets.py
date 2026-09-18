"""Load digitized figure targets (acceptance bands) for the reproduction studies."""
from __future__ import annotations
import json
from pathlib import Path

_DIR = Path(__file__).with_name("targets")


def load_target(name: str) -> dict:
    return json.loads((_DIR / f"{name}.json").read_text())
