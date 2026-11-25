# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Run directories and manifests.

Every analysis writes into a run directory:

    <run_dir>/
      manifest.json   inputs (with sha256), parameters, versions, refdata, backends,
                      genome build, ancestry, citations, warnings, outputs
      result.json     the same short JSON summary the CLI prints
      ...             method outputs (TSV / Parquet / PNG / Markdown)

The manifest is what makes a result reproducible and what the report builder
reads to cite tools and data correctly.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import math
import platform
import sys
import tomllib
from importlib import resources
from pathlib import Path
from typing import Any

from genor import RESEARCH_USE_NOTICE, __version__
from genor.errors import GenorError, Warnings

_HASH_LIMIT = 4 * 1024**3  # don't hash inputs larger than 4 GB (record size + mtime instead)


def now() -> str:
    return _dt.datetime.now(_dt.UTC).isoformat(timespec="seconds")


def sha256_file(path: Path, limit: int = _HASH_LIMIT) -> str | None:
    if path.stat().st_size > limit:
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_citations() -> dict[str, dict[str, str]]:
    text = resources.files("genor").joinpath("citations.toml").read_text()
    return tomllib.loads(text)


def default_run_dir(command: str) -> Path:
    stamp = _dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    return Path("genor_runs") / f"{stamp}-{command.replace('.', '-')}"


def _jsonable(v: Any) -> Any:
    """Plain JSON types; non-finite floats (NaN, inf) become None so strict parsers accept it."""
    if isinstance(v, Path):
        return str(v)
    if isinstance(v, dict):
        return {k: _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if hasattr(v, "item") and callable(v.item):  # numpy scalar
        try:
            v = v.item()
        except Exception:  # noqa: BLE001
            return str(v)
    if isinstance(v, float) and not math.isfinite(v):
        return None
    return v


class Run:
    """Context manager that owns one run directory and its manifest."""

    def __init__(self, command: str, out: Path | str | None, params: dict[str, Any] | None = None):
        self.command = command
        self.dir = Path(out) if out else default_run_dir(command)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.params = {k: _jsonable(v) for k, v in (params or {}).items()}
        self.warnings = Warnings()
        self.summary: dict[str, Any] = {}
        self.outputs: dict[str, str] = {}
        self.inputs: list[dict[str, Any]] = []
        self.refdata: list[dict[str, Any]] = []
        self.backends: list[dict[str, Any]] = []
        self.citations: list[str] = ["genor"]
        self.meta: dict[str, Any] = {}
        self.next_steps: list[str] = []
        self.ask_user: list[dict[str, Any]] = []
        self.started = now()
        self.status = "running"
        self.error: dict[str, Any] | None = None

    # -- recording -------------------------------------------------------
    def add_input(self, path: Path | str, role: str, **extra: Any) -> None:
        p = Path(path)
        rec: dict[str, Any] = {"role": role, "path": str(p.resolve())}
        if p.is_file():
            rec["size"] = p.stat().st_size
            rec["sha256"] = sha256_file(p)
        rec.update({k: _jsonable(v) for k, v in extra.items()})
        self.inputs.append(rec)

    def add_refdata(self, name: str, **info: Any) -> None:
        self.refdata.append({"name": name, **_jsonable(info)})

    def add_backend(self, name: str, **info: Any) -> None:
        self.backends.append({"name": name, **_jsonable(info)})

    def cite(self, *keys: str) -> None:
        for k in keys:
            if k not in self.citations:
                self.citations.append(k)

    def warn(self, code: str, message: str, severity: str = "warning", **details: Any) -> None:
        self.warnings.add(code, message, severity, **_jsonable(details))

    def output(self, key: str, path: Path | str) -> Path:
        p = Path(path)
        self.outputs[key] = str(p)
        return p

    def path(self, name: str) -> Path:
        return self.dir / name

    def ask(self, question: str, options: list[str] | None = None, why: str = "") -> None:
        """Record a decision the agent must put to the human before relying on the result."""
        self.ask_user.append({"question": question, "options": options or [], "why": why})

    # -- serialisation ---------------------------------------------------
    def manifest(self) -> dict[str, Any]:
        cites = load_citations()
        return {
            "genor_version": __version__,
            "command": self.command,
            "argv": sys.argv,
            "status": self.status,
            "started": self.started,
            "finished": now() if self.status != "running" else None,
            "platform": {
                "python": sys.version.split()[0],
                "system": platform.system(),
                "machine": platform.machine(),
            },
            "params": self.params,
            "inputs": self.inputs,
            "refdata": self.refdata,
            "backends": self.backends,
            "meta": _jsonable(self.meta),
            "summary": _jsonable(self.summary),
            "warnings": self.warnings.to_list(),
            "ask_user": self.ask_user,
            "outputs": self.outputs,
            "citations": [
                {"key": k, **cites.get(k, {"text": k})} for k in self.citations
            ],
            "error": self.error,
            "notice": RESEARCH_USE_NOTICE,
        }

    def result(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "status": self.status,
            "command": self.command,
            "run_dir": str(self.dir),
            "summary": _jsonable(self.summary),
            "warnings": self.warnings.to_list(),
            "files": self.outputs,
        }
        if self.ask_user:
            out["ask_user"] = self.ask_user
        if self.next_steps:
            out["next_steps"] = self.next_steps
        if self.error:
            out["error"] = self.error
        return _jsonable(out)

    def write(self) -> None:
        (self.dir / "manifest.json").write_text(json.dumps(_jsonable(self.manifest()), indent=2, default=str,
                                                           allow_nan=False))
        (self.dir / "result.json").write_text(json.dumps(self.result(), indent=2, default=str, allow_nan=False))

    def __enter__(self) -> Run:
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        if exc is None:
            self.status = "ok"
        else:
            self.status = "error"
            if isinstance(exc, GenorError):
                self.error = exc.to_dict()
            else:
                self.error = {"code": "internal_error", "message": f"{type(exc).__name__}: {exc}"}
        self.write()
        return False


def read_manifest(run_dir: Path | str) -> dict[str, Any]:
    p = Path(run_dir) / "manifest.json"
    if not p.exists():
        raise GenorError(f"No manifest.json in {run_dir}", code="not_a_run_dir",
                         hint="Pass a directory written by a genor command.")
    return json.loads(p.read_text())
