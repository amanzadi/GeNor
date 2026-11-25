# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""External tool resolution.

Resolution order for a backend:
1. a path the user registered with ``genor backend add <name> <path>``
2. a binary GeNor downloaded earlier into its data dir
3. the executable on ``PATH``
"""

from __future__ import annotations

import json
import os
import platform
import shutil
import stat
import subprocess
import tomllib
import zipfile
from importlib import resources
from pathlib import Path
from typing import Any

from genor.config import bin_dir, config_dir
from genor.errors import GenorError


def registry() -> dict[str, dict[str, Any]]:
    return tomllib.loads(resources.files("genor.backends").joinpath("registry.toml").read_text())


def platform_key() -> str:
    system = platform.system().lower()
    machine = platform.machine().lower()
    machine = {"amd64": "x86_64", "aarch64": "arm64"}.get(machine, machine)
    return f"{system}-{machine}"


def _has_avx2() -> bool:
    try:
        return "avx2" in Path("/proc/cpuinfo").read_text()
    except OSError:
        return False


def _registrations_path() -> Path:
    return config_dir() / "backends.json"


def registrations() -> dict[str, str]:
    p = _registrations_path()
    return json.loads(p.read_text()) if p.exists() else {}


def _save_registrations(d: dict[str, str]) -> None:
    _registrations_path().write_text(json.dumps(d, indent=2, sort_keys=True))


def entry(name: str) -> dict[str, Any]:
    reg = registry()
    if name not in reg:
        raise GenorError(f"Unknown backend '{name}'", code="unknown_backend",
                         hint=f"Known backends: {', '.join(sorted(reg))}")
    return reg[name]


def add(name: str, path: str | Path) -> dict[str, Any]:
    e = entry(name)
    p = Path(path).expanduser().resolve()
    if p.is_dir():
        cand = p / e["executable"]
        if not cand.exists():
            raise GenorError(f"{p} does not contain '{e['executable']}'", code="backend_not_found")
        p = cand
    if not p.exists():
        raise GenorError(f"{p} does not exist", code="backend_not_found")
    if p.name != e["executable"]:
        raise GenorError(f"{p} is not '{e['executable']}'", code="backend_not_found",
                         hint=f"Pass the {name} executable named '{e['executable']}' or the directory containing it.")
    regs = registrations()
    regs[name] = str(p)
    _save_registrations(regs)
    return {"name": name, "path": str(p), "version": _version(name, p)}


def remove(name: str) -> None:
    regs = registrations()
    regs.pop(name, None)
    _save_registrations(regs)


def _downloaded(name: str) -> Path | None:
    e = entry(name)
    p = bin_dir() / name / e["executable"]
    return p if p.exists() else None


def _platform_entry(e: dict[str, Any]) -> dict[str, Any] | None:
    plats = e.get("platforms", {})
    key = platform_key()
    if key == "linux-x86_64" and _has_avx2() and f"{key}-avx2" in plats:
        return plats[f"{key}-avx2"]
    return plats.get(key)


def find(name: str) -> Path | None:
    e = entry(name)
    if e["mode"] == "python":
        import importlib.util

        spec = importlib.util.find_spec(e["module"])
        return Path(spec.origin).parent if spec and spec.origin else None
    regs = registrations()
    if name in regs and Path(regs[name]).exists():
        return Path(regs[name])
    if p := _downloaded(name):
        return p
    if w := shutil.which(e["executable"]):
        return Path(w)
    return None


def _version(name: str, path: Path) -> str | None:
    e = entry(name)
    if e["mode"] == "python":
        from importlib.metadata import PackageNotFoundError, version

        try:
            return version(e["module"])
        except PackageNotFoundError:
            return None
    args = e.get("version_args")
    if not args:
        return None
    try:
        r = subprocess.run([str(path), *args], capture_output=True, text=True, timeout=20)
        text = (r.stdout or r.stderr).strip().splitlines()
        for line in text:
            if line.strip():
                return line.strip()[:120]
    except (OSError, subprocess.SubprocessError):
        return None
    return None


def require(name: str) -> Path:
    p = find(name)
    if p:
        return p
    e = entry(name)
    if e["mode"] == "download" and _platform_entry(e):
        hint = f"Run: genor backend install {name}  (downloads a pinned {name} binary)"
    elif e["mode"] == "download":
        hint = f"No pinned {name} build for {platform_key()}; install it yourself and run: genor backend add {name} <path>"
    else:
        hint = e.get("install", f"genor backend add {name} <path>")
    raise GenorError(f"Backend '{name}' is not available.", code="backend_missing", hint=hint,
                     backend=name)


def install(name: str) -> dict[str, Any]:
    """Download a pinned, redistributable binary."""
    import pooch

    e = entry(name)
    if e["mode"] != "download":
        raise GenorError(f"'{name}' cannot be downloaded by GeNor ({e['licence']}).",
                         code="not_downloadable", hint=e.get("install"))
    pe = _platform_entry(e)
    if not pe:
        raise GenorError(f"No pinned {name} build for {platform_key()}.", code="no_platform_build",
                         hint=f"Install it yourself and run: genor backend add {name} <path>")
    dest = bin_dir() / name
    dest.mkdir(parents=True, exist_ok=True)
    archive = Path(pooch.retrieve(pe["url"], known_hash=f"sha256:{pe['sha256']}", path=dest,
                                  fname=Path(pe["url"]).name, progressbar=False))
    with zipfile.ZipFile(archive) as zf:
        zf.extractall(dest / "unpacked")
    member = pe["member"]
    if member.startswith("auto:"):
        target_name = member.split(":", 1)[1]
        hits = [p for p in (dest / "unpacked").rglob(target_name) if "__MACOSX" not in p.parts]
        if not hits:
            raise GenorError(f"{target_name} not found in {archive.name}", code="backend_archive")
        src = hits[0]
    else:
        src = dest / "unpacked" / member
    exe = dest / e["executable"]
    if exe.exists() or exe.is_symlink():
        exe.unlink()
    # keep the unpacked tree (some builds ship dylibs next to the binary) and symlink
    exe.symlink_to(src.resolve())
    src.chmod(src.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    archive.unlink()
    return {"name": name, "path": str(exe), "version": _version(name, exe)}


def status() -> list[dict[str, Any]]:
    out = []
    for name, e in registry().items():
        p = find(name)
        rec: dict[str, Any] = {
            "name": name,
            "available": p is not None,
            "mode": e["mode"],
            "redistributable": e.get("redistributable", False),
            "provides": e.get("provides", []),
            "description": e["description"],
        }
        if p:
            rec["path"] = str(p)
            rec["version"] = _version(name, p)
        else:
            if e["mode"] == "download":
                rec["how_to_enable"] = (f"genor backend install {name}" if _platform_entry(e)
                                        else f"no pinned build for {platform_key()}; genor backend add {name} <path>")
            else:
                rec["how_to_enable"] = e.get("install", f"genor backend add {name} <path>")
        out.append(rec)
    return out


def run(name: str, args: list[str], cwd: Path | None = None, log: Path | None = None,
        timeout: int | None = None) -> subprocess.CompletedProcess[str]:
    exe = require(name)
    env = dict(os.environ)
    r = subprocess.run([str(exe), *args], cwd=cwd, capture_output=True, text=True,
                       timeout=timeout, env=env)
    if log:
        log.write_text(f"$ {exe} {' '.join(args)}\n\n--- stdout ---\n{r.stdout}\n--- stderr ---\n{r.stderr}")
    if r.returncode != 0:
        tail = (r.stderr or r.stdout).strip().splitlines()[-15:]
        raise GenorError(f"{name} exited with code {r.returncode}", code="backend_failed",
                         hint="\n".join(tail), log=str(log) if log else None)
    return r
