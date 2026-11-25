# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Detached background jobs.

A job is ``python -m genor <args> --json`` started in its own session, with stdout,
stderr and a ``job.json`` record in ``<jobs_dir>/<id>/``. MCP tools return the id at
once (Codex's default MCP tool timeout is 60 s); agents poll ``genor job status``.
"""

from __future__ import annotations

import json
import os
import re
import signal
import subprocess
import sys
import time
import uuid
from pathlib import Path
from typing import Any

from genor.config import jobs_dir
from genor.errors import GenorError
from genor.runs import now

_JOB_ID = re.compile(r"^\d{8}-\d{6}-[0-9a-f]{6}$")


def _dir(job_id: str) -> Path:
    d = jobs_dir() / job_id
    if not _JOB_ID.fullmatch(job_id) or not d.exists():
        raise GenorError(f"No job {job_id}", code="job_not_found", hint="genor job list")
    return d


def submit(args: list[str], cwd: Path | None = None, label: str | None = None) -> dict[str, Any]:
    job_id = time.strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:6]
    d = jobs_dir() / job_id
    d.mkdir(parents=True)
    cmd = [sys.executable, "-m", "genor", *args]
    if "--json" not in args:
        cmd.append("--json")
    rec = {"id": job_id, "args": args, "label": label or " ".join(args[:2]),
           "cwd": str((cwd or Path.cwd()).resolve()), "submitted": now(), "status": "running"}
    with open(d / "stdout.log", "w") as out, open(d / "stderr.log", "w") as err:
        p = subprocess.Popen(cmd, cwd=rec["cwd"], stdout=out, stderr=err, stdin=subprocess.DEVNULL,
                             start_new_session=True)
    rec["pid"] = p.pid
    (d / "job.json").write_text(json.dumps(rec, indent=2))
    # completion is detected lazily in status(): process gone + result JSON on stdout
    return {"job_id": job_id, "status": "running", "pid": p.pid}


def _alive(pid: int) -> bool:
    try:
        r = os.waitpid(pid, os.WNOHANG)
        if r != (0, 0):
            return False
    except ChildProcessError:
        pass
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    # zombie check on Linux
    stat = Path(f"/proc/{pid}/stat")
    if stat.exists() and stat.read_text().split()[2] == "Z":
        return False
    return True


def _result_from_stdout(d: Path) -> dict[str, Any] | None:
    text = (d / "stdout.log").read_text().strip()
    if not text:
        return None
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        # last JSON object in the stream
        start = text.rfind("\n{")
        try:
            return json.loads(text[start + 1:] if start >= 0 else text)
        except json.JSONDecodeError:
            return None


def status(job_id: str) -> dict[str, Any]:
    d = _dir(job_id)
    rec = json.loads((d / "job.json").read_text())
    # a result on stdout means the job finished even if its PID now belongs to another process
    if rec["status"] == "running" and (not _alive(rec["pid"]) or _result_from_stdout(d) is not None):
        res = _result_from_stdout(d)
        if res is None:
            rec["status"] = "failed"
        else:
            rec["status"] = "done" if res.get("status") == "ok" else "failed"
            rec["run_dir"] = res.get("run_dir")
        rec["finished"] = now()
        (d / "job.json").write_text(json.dumps(rec, indent=2))
    out = {k: rec.get(k) for k in ("id", "label", "status", "submitted", "finished", "run_dir")}
    if out["status"] == "running":
        tail = (d / "stderr.log").read_text().strip().splitlines()[-3:]
        if tail:
            out["last_log_lines"] = tail
    return out


def result(job_id: str) -> dict[str, Any]:
    st = status(job_id)
    if st["status"] == "running":
        return {"job_id": job_id, "status": "running", "hint": "Poll again later."}
    res = _result_from_stdout(_dir(job_id))
    if res is None:
        err = (_dir(job_id) / "stderr.log").read_text().strip().splitlines()[-20:]
        return {"job_id": job_id, "status": "failed", "stderr_tail": err}
    return {"job_id": job_id, **res}


def logs(job_id: str, n: int = 50) -> dict[str, Any]:
    d = _dir(job_id)
    return {"job_id": job_id,
            "stderr": (d / "stderr.log").read_text().splitlines()[-n:],
            "stdout": (d / "stdout.log").read_text().splitlines()[-n:]}


def cancel(job_id: str) -> dict[str, Any]:
    d = _dir(job_id)
    status(job_id)  # marks a finished job done/failed, so a reused PID is never signalled
    rec = json.loads((d / "job.json").read_text())
    if rec["status"] == "running":
        os.killpg(os.getpgid(rec["pid"]), signal.SIGTERM)
        rec["status"] = "cancelled"
        rec["finished"] = now()
        (d / "job.json").write_text(json.dumps(rec, indent=2))
    return {"job_id": job_id, "status": rec["status"]}


def list_jobs(limit: int = 20) -> list[dict[str, Any]]:
    ids = sorted((p.name for p in jobs_dir().iterdir() if (p / "job.json").exists()), reverse=True)
    return [status(i) for i in ids[:limit]]


def wait(job_id: str, timeout: float) -> dict[str, Any]:
    t0 = time.monotonic()
    while time.monotonic() - t0 < timeout:
        st = status(job_id)
        if st["status"] != "running":
            return result(job_id)
        time.sleep(0.5)
    return {"job_id": job_id, "status": "running",
            "hint": f"Still running after {timeout:.0f}s. Poll with job_status / genor job status {job_id}."}
