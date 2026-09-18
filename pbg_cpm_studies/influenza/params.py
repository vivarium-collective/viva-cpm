"""Load the canonical Sego-2022 parameter set from params.yaml."""
from __future__ import annotations
from functools import lru_cache
from pathlib import Path
import yaml

_PARAMS_PATH = Path(__file__).with_name("params.yaml")


@lru_cache(maxsize=1)
def load_params() -> dict:
    return yaml.safe_load(_PARAMS_PATH.read_text())


PARAMS = load_params()
