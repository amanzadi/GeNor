# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Combine one or more run directories into a Markdown (and HTML) report.

The report states what was run, key numbers, every warning, the decisions that were
left to the human, figures, and the citations each run recorded.
"""

from __future__ import annotations

import html
import json
import re
from pathlib import Path
from typing import Any

from genor import RESEARCH_USE_NOTICE
from genor.runs import read_manifest


def _flat(d: Any, prefix: str = "", depth: int = 0) -> list[tuple[str, Any]]:
    out: list[tuple[str, Any]] = []
    if isinstance(d, dict) and depth < 2:
        for k, v in d.items():
            out += _flat(v, f"{prefix}{k}." if prefix or True else k, depth + 1)
    elif isinstance(d, (list, dict)):
        s = json.dumps(d, default=str)
        out.append((prefix.rstrip("."), s if len(s) < 200 else s[:197] + "..."))
    else:
        out.append((prefix.rstrip("."), d))
    return out


def _fmt(v: Any) -> str:
    if isinstance(v, float):
        return f"{v:.4g}"
    return str(v)


def section(run_dir: Path, out_dir: Path) -> tuple[str, list[dict[str, Any]]]:
    m = read_manifest(run_dir)
    rel = lambda p: Path(p).resolve().relative_to(out_dir.resolve()) if Path(p).resolve().is_relative_to(out_dir.resolve()) else Path(p).resolve()  # noqa: E731
    lines = [f"## {m['command']}  ", f"Run directory: `{run_dir}` · status **{m['status']}** · "
             f"GeNor {m['genor_version']} · {m.get('finished') or m['started']}", ""]
    if m.get("inputs"):
        lines.append("**Inputs**")
        for i in m["inputs"]:
            lines.append(f"- {i['role']}: `{i['path']}`" + (f" (sha256 {i['sha256'][:12]}…)" if i.get("sha256") else ""))
        lines.append("")
    if m.get("refdata") or m.get("backends"):
        lines.append("**Reference data and tools**")
        for r in m.get("refdata", []):
            lines.append(f"- {r['name']}: {r.get('build', '')} {r.get('ancestry', '')} ({r.get('licence', '')})")
        for b in m.get("backends", []):
            lines.append(f"- {b['name']} {b.get('version') or ''}")
        lines.append("")
    if m.get("params"):
        lines.append("**Parameters:** " + ", ".join(f"`{k}={_fmt(v)}`" for k, v in m["params"].items()
                                                    if v is not None))
        lines.append("")
    summ = [(k, v) for k, v in _flat(m.get("summary", {})) if not k.endswith("_path")]
    if summ:
        lines += ["| Result | Value |", "|---|---|"]
        lines += [f"| {k} | {_fmt(v)} |" for k, v in summ[:60]]
        lines.append("")
    if m.get("warnings"):
        lines.append("**Warnings**")
        for w in m["warnings"]:
            lines.append(f"- [{w['severity']}] `{w['code']}`: {w['message']}")
        lines.append("")
    if m.get("ask_user"):
        lines.append("**Decisions for the human**")
        for a in m["ask_user"]:
            lines.append(f"- {a['question']}" + (f" ({a['why']})" if a.get("why") else ""))
        lines.append("")
    for key, path in m.get("outputs", {}).items():
        p = Path(path)
        if p.suffix == ".png" and p.exists():
            lines.append(f"![{key}]({rel(p)})")
    if (run_dir / "strobe_mr.md").exists():
        lines += ["", (run_dir / "strobe_mr.md").read_text().replace("# STROBE-MR", "### STROBE-MR")]
    if m.get("error"):
        lines.append(f"**Error:** `{m['error'].get('code')}`: {m['error'].get('message')}")
    return "\n".join(lines) + "\n", m.get("citations", [])


def build(run_dirs: list[Path], out: Path, title: str = "GeNor analysis report") -> dict[str, Any]:
    out.parent.mkdir(parents=True, exist_ok=True)
    parts = [f"# {title}", "", f"> {RESEARCH_USE_NOTICE}", ""]
    cites: dict[str, str] = {}
    for d in run_dirs:
        sec, c = section(Path(d), out.parent)
        parts.append(sec)
        for x in c:
            cites.setdefault(x["key"], x.get("text", x["key"]))
    parts += ["## References", ""] + [f"{i}. {t}" for i, t in enumerate(cites.values(), 1)]
    md = "\n".join(parts) + "\n"
    out.write_text(md)
    html_path = out.with_suffix(".html")
    html_path.write_text(to_html(md, title))
    return {"markdown": str(out), "html": str(html_path), "n_runs": len(run_dirs), "n_citations": len(cites)}


def to_html(md: str, title: str) -> str:
    """Minimal Markdown -> HTML (headings, tables, lists, images, code, bold)."""
    body: list[str] = []
    in_table = in_list = False
    for line in md.splitlines():
        esc = html.escape(line)
        esc = re.sub(r"`([^`]+)`", r"<code>\1</code>", esc)
        esc = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", esc)
        esc = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)", r'<img alt="\1" src="\2" style="max-width:100%">', esc)
        if line.startswith("|"):
            cells = [c.strip() for c in esc.strip("|").split("|")]
            if set(line.replace("|", "").strip()) <= {"-", " "}:
                continue
            if not in_table:
                body.append("<table>")
                in_table = True
                body.append("<tr>" + "".join(f"<th>{c}</th>" for c in cells) + "</tr>")
            else:
                body.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")
            continue
        if in_table:
            body.append("</table>")
            in_table = False
        if line.startswith(("- ", "1. ")) or re.match(r"^\d+\. ", line):
            if not in_list:
                body.append("<ul>")
                in_list = True
            item = re.sub(r"^(- |\d+\. )", "", esc)
            body.append(f"<li>{item}</li>")
            continue
        if in_list:
            body.append("</ul>")
            in_list = False
        if m := re.match(r"^(#{1,4}) (.*)", line):
            n = len(m.group(1))
            body.append(f"<h{n}>{html.escape(m.group(2))}</h{n}>")
        elif line.startswith("&gt;") or line.startswith(">"):
            body.append(f"<blockquote>{esc.lstrip('&gt;').strip()}</blockquote>")
        elif line.strip():
            body.append(f"<p>{esc}</p>")
    if in_table:
        body.append("</table>")
    if in_list:
        body.append("</ul>")
    css = ("body{font-family:system-ui,sans-serif;max-width:980px;margin:2rem auto;padding:0 1rem;"
           "line-height:1.5;color:#1d1d1f}table{border-collapse:collapse;margin:1rem 0;font-size:.9rem}"
           "td,th{border:1px solid #ddd;padding:.3rem .6rem;text-align:left;vertical-align:top}"
           "code{background:#f3f3f3;padding:0 .2rem}blockquote{border-left:4px solid #c33;margin:0;"
           "padding:.3rem 1rem;background:#fff6f6}")
    return (f"<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title>"
            f"<style>{css}</style></head><body>{''.join(body)}</body></html>")
