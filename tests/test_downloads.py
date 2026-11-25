# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Offline tests for refdata.fetch and backends.install (pooch.retrieve is faked; no network)."""

import io
import shutil
import tarfile
import zipfile
from pathlib import Path

import pooch
import pytest

from genor import backends, refdata
from genor.errors import GenorError

URL = "https://example.invalid/pkg.bin"
CALLS: list[str] = []


@pytest.fixture
def served(tmp_path, monkeypatch):
    """Fake pooch.retrieve: copy a local file for URL, enforcing known_hash like pooch does."""
    src = tmp_path / "served.bin"
    CALLS.clear()

    def retrieve(url, known_hash=None, fname=None, path=None, progressbar=False, **kw):
        CALLS.append(url)
        out = Path(path) / (fname or "download")
        shutil.copy(src, out)
        if known_hash:
            alg, _, want = known_hash.rpartition(":")
            if pooch.file_hash(str(out), alg=alg or "sha256") != want:
                out.unlink()
                raise ValueError(f"hash mismatch for {url}")
        return str(out)

    monkeypatch.setattr(pooch, "retrieve", retrieve)
    return src


def _tgz(path, members):
    with tarfile.open(path, "w:gz") as tf:
        for name, data in members.items():
            ti = tarfile.TarInfo(name)
            ti.size = len(data)
            tf.addfile(ti, io.BytesIO(data))


def _ref(monkeypatch, src, size=10, hash_="auto"):
    h = pooch.file_hash(str(src), alg="sha256") if hash_ == "auto" else hash_
    e = {"url": URL, "size": size, "hash": h, "unpack": "tgz", "filename": "pkg.tgz"}
    monkeypatch.setattr(refdata, "registry", lambda: {"toy": e})


def test_refdata_fetch_unpacks_and_is_idempotent(served, monkeypatch):
    _tgz(served, {"a/x.txt": b"hello"})
    _ref(monkeypatch, served)
    r = refdata.fetch("toy")
    assert r["downloaded"] is True
    assert (refdata.dataset_dir("toy") / "a" / "x.txt").read_bytes() == b"hello"
    assert (refdata.dataset_dir("toy") / ".genor-complete.json").exists()
    assert refdata.is_installed("toy")
    assert not (refdata.dataset_dir("toy") / "pkg.tgz").exists()
    assert refdata.fetch("toy")["downloaded"] is False
    assert len(CALLS) == 1


def test_refdata_wrong_hash_leaves_no_marker(served, monkeypatch):
    _tgz(served, {"x.txt": b"hello"})
    _ref(monkeypatch, served, hash_="sha256:" + "0" * 64)
    with pytest.raises(ValueError):
        refdata.fetch("toy")
    assert not refdata.is_installed("toy")


def test_refdata_unpinned_needs_flag(served, monkeypatch):
    _tgz(served, {"x.txt": b"hello"})
    _ref(monkeypatch, served, hash_="")
    with pytest.raises(GenorError) as ei:
        refdata.fetch("toy")
    assert ei.value.code == "unpinned_refdata"
    assert not CALLS
    assert refdata.fetch("toy", allow_unpinned=True)["downloaded"] is True


def test_refdata_large_needs_confirmation(served, monkeypatch):
    _tgz(served, {"x.txt": b"hello"})
    _ref(monkeypatch, served, size=refdata.LARGE_BYTES + 1)
    with pytest.raises(GenorError) as ei:
        refdata.fetch("toy")
    assert ei.value.code == "confirm_download"
    assert not CALLS and not refdata.is_installed("toy")
    assert refdata.fetch("toy", yes=True)["downloaded"] is True


def test_refdata_tar_path_traversal_contained(served, monkeypatch):
    _tgz(served, {"../evil": b"pwned", "ok.txt": b"fine"})
    _ref(monkeypatch, served)
    try:
        refdata.fetch("toy")
    except Exception:
        pass  # refusing the archive is acceptable
    assert not (refdata.root() / "evil").exists()
    assert not (refdata.root().parent / "evil").exists()


def test_backend_install_symlinks_executable(served, monkeypatch):
    with zipfile.ZipFile(served, "w") as zf:
        zf.writestr("pkg/fakebin", "#!/bin/sh\necho fakebin 1.0\n")
    sha = pooch.file_hash(str(served), alg="sha256")
    e = {"mode": "download", "executable": "fakebin", "licence": "MIT", "description": "toy",
         "version_args": ["--version"],
         "platforms": {backends.platform_key(): {"url": URL, "sha256": sha, "member": "auto:fakebin"}}}
    monkeypatch.setattr(backends, "registry", lambda: {"fakebin": e})
    monkeypatch.setattr(backends.shutil, "which", lambda _n: None)
    r = backends.install("fakebin")
    exe = Path(r["path"])
    assert exe.exists() and exe.resolve().stat().st_mode & 0o111
    assert backends.find("fakebin") == exe
