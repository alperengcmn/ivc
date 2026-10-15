"""Load project configuration (config/*.yaml)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = ROOT / "config"


def load_yaml(path: Path) -> dict[str, Any]:
    import yaml  # analysis dependency; validators do not need it

    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def parameters(path: Path | None = None) -> dict[str, Any]:
    return load_yaml(path or CONFIG_DIR / "parameters.yaml")
