# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Skill files: locate, serve as guidance, and install for each agent.

Skills live in the repository's ``skills/`` directory (the Claude Code plugin reads them
from there) and are copied into the wheel as ``genor/skills``.

Install targets:
- Claude Code: ``.claude/skills`` (project) or ``~/.claude/skills`` (user)
- Codex, opencode: ``.agents/skills`` (project) or ``~/.agents/skills`` (user)
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path
from typing import Any

from genor.errors import GenorError

TARGETS = {
    "claude": (".claude/skills", "~/.claude/skills"),
    "codex": (".agents/skills", "~/.agents/skills"),
    "opencode": (".agents/skills", "~/.agents/skills"),
}


def skills_root() -> Path:
    here = Path(__file__).resolve().parent
    for cand in (here / "skills", here.parents[1] / "skills"):
        if (cand / "genor" / "SKILL.md").exists():
            return cand
    raise GenorError("Skill files not found in this installation.", code="skills_missing")


def _front(text: str) -> dict[str, str]:
    m = re.match(r"^---\n(.*?)\n---\n", text, re.S)
    out: dict[str, str] = {}
    if m:
        for line in m.group(1).splitlines():
            if ":" in line:
                k, v = line.split(":", 1)
                out[k.strip()] = v.strip().strip('"')
    return out


def list_skills() -> list[dict[str, Any]]:
    root = skills_root()
    out = []
    for d in sorted(p for p in root.iterdir() if (p / "SKILL.md").exists()):
        fm = _front((d / "SKILL.md").read_text())
        refs = sorted(str(p.relative_to(d)) for p in (d / "references").glob("*.md")) \
            if (d / "references").exists() else []
        out.append({"name": d.name, "description": fm.get("description", ""), "references": refs})
    return out


def guidance(topic: str | None = None) -> dict[str, Any]:
    """Return skill text for agents without native skills support (MCP get_guidance)."""
    root = skills_root()
    if not topic:
        return {"status": "ok", "topics": list_skills(),
                "hint": "Call again with topic='<name>' or '<name>/references/<file>.md'. "
                        "Start with topic='genor'."}
    path = root / topic
    if path.is_dir():
        path = path / "SKILL.md"
    elif not path.suffix:
        path = path.with_suffix(".md")
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError as e:
        raise GenorError("Invalid topic", code="bad_topic") from e
    if not path.exists():
        raise GenorError(f"No guidance topic '{topic}'", code="bad_topic",
                         hint=", ".join(s["name"] for s in list_skills()))
    return {"status": "ok", "topic": topic, "text": path.read_text()}


def install(agent: str, scope: str = "project", dest: Path | None = None, force: bool = False) -> dict[str, Any]:
    if agent not in TARGETS:
        raise GenorError(f"Unknown agent '{agent}'", code="bad_param", hint="claude | codex | opencode")
    proj, user = TARGETS[agent]
    target = dest or (Path(user).expanduser() if scope == "user" else Path.cwd() / proj)
    target.mkdir(parents=True, exist_ok=True)
    root = skills_root()
    installed, skipped = [], []
    for d in sorted(p for p in root.iterdir() if (p / "SKILL.md").exists()):
        t = target / d.name
        if t.exists() and not force:
            skipped.append(d.name)
            continue
        if t.exists():
            shutil.rmtree(t)
        shutil.copytree(d, t)
        installed.append(d.name)
    return {"status": "ok", "agent": agent, "path": str(target), "installed": installed,
            "skipped_existing": skipped,
            "hint": "Use --force to overwrite." if skipped else None}
