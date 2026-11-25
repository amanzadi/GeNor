"""Reference data registry and lazy, checksum-pinned downloads.

Nothing is downloaded implicitly. Methods call :func:`require`, which raises a
``refdata_missing`` error naming the dataset and its size; the agent then asks
the user and runs ``genor refdata fetch <name>``.
"""

from __future__ import annotations

import json
import shutil
import tarfile
import tomllib
import zipfile
from importlib import resources
from pathlib import Path
from typing import Any

from genor.config import data_dir
from genor.errors import GenorError

LARGE_BYTES = 100 * 1024**2  # downloads above this need explicit confirmation


def registry() -> dict[str, dict[str, Any]]:
    text = resources.files("genor.refdata").joinpath("registry.toml").read_text()
    return tomllib.loads(text)


def root() -> Path:
    p = data_dir() / "refdata"
    p.mkdir(parents=True, exist_ok=True)
    return p


def dataset_dir(name: str) -> Path:
    return root() / name


def human_size(n: int | float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return str(n)


def _marker(name: str) -> Path:
    return dataset_dir(name) / ".genor-complete.json"


def is_installed(name: str) -> bool:
    return _marker(name).exists()


def entry(name: str) -> dict[str, Any]:
    reg = registry()
    if name not in reg:
        raise GenorError(
            f"Unknown reference dataset '{name}'",
            code="unknown_refdata",
            hint=f"Known datasets: {', '.join(sorted(reg))}",
        )
    return reg[name]


def status() -> list[dict[str, Any]]:
    out = []
    for name, e in registry().items():
        rec = {
            "name": name,
            "installed": is_installed(name),
            "size": human_size(e["size"]),
            "build": e.get("build"),
            "ancestry": e.get("ancestry"),
            "pinned": bool(e.get("hash")),
            "provides": e.get("provides", []),
            "description": e.get("description", ""),
        }
        if rec["installed"]:
            rec["path"] = str(dataset_dir(name))
        out.append(rec)
    return out


def fetch(name: str, yes: bool = False, allow_unpinned: bool = False) -> dict[str, Any]:
    """Download (and unpack) a dataset. Idempotent."""
    import pooch

    e = entry(name)
    target = dataset_dir(name)
    if is_installed(name):
        return {"name": name, "path": str(target), "downloaded": False}
    if e["size"] > LARGE_BYTES and not yes:
        raise GenorError(
            f"'{name}' is {human_size(e['size'])}. Large downloads need confirmation.",
            code="confirm_download",
            hint=f"Ask the user, then run: genor refdata fetch {name} --yes",
            size=human_size(e["size"]),
        )
    known_hash = e.get("hash") or None
    if not known_hash and not allow_unpinned:
        raise GenorError(
            f"'{name}' has no pinned checksum in GeNor's registry.",
            code="unpinned_refdata",
            hint=f"If the user accepts an unverified download: genor refdata fetch {name} --allow-unpinned --yes",
        )
    target.mkdir(parents=True, exist_ok=True)
    fname = e.get("filename") or f"{name}.{ {'tgz': 'tgz', 'zip': 'zip'}.get(e.get('unpack', ''), 'dat') }"
    path = pooch.retrieve(
        url=e["url"],
        known_hash=known_hash,
        fname=fname,
        path=target,
        progressbar=False,
    )
    path = Path(path)
    seen_hash = pooch.file_hash(str(path), alg="sha256")
    unpack = e.get("unpack", "")
    if unpack == "tgz":
        with tarfile.open(path) as tf:
            if hasattr(tarfile, "data_filter"):
                tf.extractall(target, filter="data")
            else:  # Python < 3.11.4
                root_ = target.resolve()
                for m in tf.getmembers():
                    if not (target / m.name).resolve().is_relative_to(root_) or m.issym() or m.islnk():
                        raise GenorError(f"Unsafe path in archive: {m.name}", code="refdata_corrupt")
                tf.extractall(target)
        path.unlink()
    elif unpack == "zip":
        with zipfile.ZipFile(path) as zf:
            zf.extractall(target)
        path.unlink()
    _marker(name).write_text(json.dumps({"url": e["url"], "hash": known_hash, "sha256_seen": seen_hash}))
    return {"name": name, "path": str(target), "downloaded": True, "sha256": seen_hash}


def remove(name: str) -> None:
    entry(name)
    shutil.rmtree(dataset_dir(name), ignore_errors=True)


def require(name: str) -> Path:
    """Return the dataset directory or raise an actionable error."""
    e = entry(name)
    if not is_installed(name):
        raise GenorError(
            f"Reference dataset '{name}' is not downloaded ({human_size(e['size'])}).",
            code="refdata_missing",
            hint=(
                f"Ask the user before downloading, then run: genor refdata fetch {name}"
                + (" --yes" if e["size"] > LARGE_BYTES else "")
            ),
            dataset=name,
            size=human_size(e["size"]),
        )
    return dataset_dir(name)


def info_for_manifest(name: str) -> dict[str, Any]:
    e = entry(name)
    rec = {"url": e["url"], "hash": e.get("hash"), "build": e.get("build"),
           "ancestry": e.get("ancestry"), "licence": e.get("licence")}
    m = _marker(name)
    if m.exists():
        rec["sha256_seen"] = json.loads(m.read_text()).get("sha256_seen")
    return rec


def ldscore_prefix(name: str) -> str:
    """Locate the per-chromosome LD score prefix inside a dataset (``<prefix>{chr}.l2.ldscore.gz``)."""
    d = require(name)
    e = entry(name)
    pre = e.get("prefix", "auto")
    if pre != "auto":
        return str(d / pre)
    hits = sorted(d.rglob("*22.l2.ldscore.gz"))
    if not hits:
        raise GenorError(f"No LD score files found in {d}", code="refdata_corrupt",
                         hint=f"genor refdata remove {name} && genor refdata fetch {name}")
    return str(hits[0]).removesuffix("22.l2.ldscore.gz")


def plink_prefix(name: str) -> str:
    """Locate the per-chromosome PLINK prefix (``<prefix>{chr}.bed``)."""
    d = require(name)
    hits = sorted(d.rglob("*.22.bed"))
    if not hits:
        raise GenorError(f"No PLINK .bed files found in {d}", code="refdata_corrupt",
                         hint=f"genor refdata remove {name} && genor refdata fetch {name}")
    return str(hits[0]).removesuffix("22.bed")


def file_in(name: str) -> Path:
    d = require(name)
    e = entry(name)
    return d / e["filename"]
