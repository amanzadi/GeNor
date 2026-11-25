"""Locations and settings.

Override with environment variables:

- ``GENOR_HOME``: base directory for everything below (useful on clusters).
- ``GENOR_DATA_DIR``: reference data cache (large; LDSC, LD panels).
- ``GENOR_CONFIG_DIR``: backend registrations (``backends.json``).
- ``GENOR_JOBS_DIR``: background job records.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from platformdirs import user_config_dir, user_data_dir, user_state_dir

APP = "genor"


def _base(env: str, default: str) -> Path:
    if v := os.environ.get(env):
        return Path(v).expanduser()
    if home := os.environ.get("GENOR_HOME"):
        return Path(home).expanduser() / Path(default).name
    return Path(default)


def data_dir() -> Path:
    p = _base("GENOR_DATA_DIR", user_data_dir(APP))
    p.mkdir(parents=True, exist_ok=True)
    return p


def config_dir() -> Path:
    p = _base("GENOR_CONFIG_DIR", user_config_dir(APP))
    p.mkdir(parents=True, exist_ok=True)
    return p


def jobs_dir() -> Path:
    p = _base("GENOR_JOBS_DIR", str(Path(user_state_dir(APP)) / "jobs"))
    p.mkdir(parents=True, exist_ok=True)
    return p


def bin_dir() -> Path:
    p = data_dir() / "bin"
    p.mkdir(parents=True, exist_ok=True)
    return p


def _settings_path() -> Path:
    return config_dir() / "settings.json"


def load_settings() -> dict[str, Any]:
    p = _settings_path()
    if p.exists():
        return json.loads(p.read_text())
    return {}


def save_settings(settings: dict[str, Any]) -> None:
    _settings_path().write_text(json.dumps(settings, indent=2, sort_keys=True))
