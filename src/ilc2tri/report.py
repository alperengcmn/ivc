"""Markdown findings collector and report output shared by the validators."""

from __future__ import annotations

import sys
from collections.abc import Iterable
from pathlib import Path


def codes(items: Iterable[str]) -> str:
    return ", ".join(f"`{x}`" for x in items)


class Findings:
    """Ordered list of Markdown bullets with an error counter."""

    def __init__(self) -> None:
        self.lines: list[str] = []
        self.errors = 0

    @property
    def passed(self) -> bool:
        return self.errors == 0

    def error(self, message: str) -> None:
        self.errors += 1
        self.lines.append(f"- **ERROR:** {message}")

    def info(self, message: str) -> None:
        self.lines.append(f"- {message}")

    def note(self, message: str) -> None:
        self.lines.append(f"- **INFO:** {message}")

    def unresolved(self, message: str) -> None:
        self.lines.append(f"- **UNRESOLVED:** {message}")

    def heading(self, text: str) -> None:
        self.lines += ["", f"### {text}"]


def status_line(label: str, passed: bool) -> str:
    """E.g. 'STRUCTURAL METADATA VALIDATION: PASS' — the label scopes what passed."""
    return f"{label}: {'PASS' if passed else 'FAIL'}"


def emit(report: str, path: Path | None, status: str, passed: bool) -> int:
    """Write or print the report. Exit codes: 0 pass, 1 fail, 2 I/O error."""
    if path is None:
        print(report, end="")
    else:
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(report, encoding="utf-8")
        except OSError as exc:
            print(f"Error: could not write report to {path}: {exc}", file=sys.stderr)
            return 2
        print(f"{status} — report saved to {path}")
    return 0 if passed else 1
