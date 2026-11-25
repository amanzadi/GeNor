# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi

"""Structured errors and warnings.

Every failure an agent can hit carries a stable ``code`` and a ``hint`` that says
what to do next, so an agent can branch on the code instead of parsing prose.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


class GenorError(Exception):
    """An error the user or agent can act on."""

    def __init__(self, message: str, code: str = "error", hint: str | None = None, **details: Any):
        super().__init__(message)
        self.message = message
        self.code = code
        self.hint = hint
        self.details = details

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {"code": self.code, "message": self.message}
        if self.hint:
            out["hint"] = self.hint
        if self.details:
            out["details"] = self.details
        return out


class IndividualLevelDataError(GenorError):
    """Raised when a study input looks like individual-level genotype data."""


class NeedsUserDecision(GenorError):
    """Raised when an analysis cannot continue without a choice that belongs to the human.

    Agents must stop and ask the user; they must not guess.
    """

    def __init__(self, message: str, question: str, options: list[str] | None = None, **details: Any):
        super().__init__(message, code="needs_user_decision", hint=question, **details)
        self.question = question
        self.options = options or []

    def to_dict(self) -> dict[str, Any]:
        out = super().to_dict()
        out["ask_user"] = {"question": self.question, "options": self.options}
        return out


@dataclass
class Warning_:
    code: str
    message: str
    severity: str = "warning"  # info | warning | serious
    details: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        if not d["details"]:
            d.pop("details")
        return d


class Warnings:
    """Collects structured warnings during a command."""

    def __init__(self) -> None:
        self.items: list[Warning_] = []

    def add(self, code: str, message: str, severity: str = "warning", **details: Any) -> None:
        self.items.append(Warning_(code, message, severity, details))

    def to_list(self) -> list[dict[str, Any]]:
        return [w.to_dict() for w in self.items]

    def __len__(self) -> int:
        return len(self.items)
