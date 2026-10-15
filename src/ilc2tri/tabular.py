"""Strict reader for small delimited metadata tables (CSV / TSV)."""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

SUPPORTED_SUFFIXES = frozenset({".csv", ".tsv", ".txt"})
DELIMITER_CHOICES = ("auto", "comma", "tab")
_EXPLICIT = {"comma": ",", "tab": "\t"}


@dataclass
class Table:
    path: Path
    headers: list[str]
    rows: list[dict[str, str]]
    issues: list[str] = field(default_factory=list)  # structural problems, reported as errors

    def duplicate_headers(self) -> list[str]:
        return sorted(h for h, n in Counter(self.headers).items() if h and n > 1)

    def column(self, name: str) -> list[str]:
        return [cell(row, name) for row in self.rows]


def cell(row: dict[str, str], column: str) -> str:
    return row.get(column, "").strip()


def _resolve_delimiter(path: Path, option: str, sample: str) -> str:
    # An explicit --delimiter always wins; then the file extension; then sniffing.
    if option in _EXPLICIT:
        return _EXPLICIT[option]
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return ","
    if suffix == ".tsv":
        return "\t"
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t;").delimiter
    except csv.Error as exc:
        raise ValueError("Could not detect delimiter; pass --delimiter comma or tab.") from exc


def read_table(path: Path, delimiter: str = "auto") -> Table:
    """Read a header + rows table, recording (not hiding) malformed rows."""
    if delimiter not in DELIMITER_CHOICES:
        raise ValueError(f"Unknown delimiter option: {delimiter!r}")
    if not path.is_file():
        raise FileNotFoundError(f"Metadata file not found: {path}")
    if path.suffix.lower() not in SUPPORTED_SUFFIXES:
        raise ValueError("Unsupported file extension; use .csv, .tsv or .txt.")

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(8192)
        handle.seek(0)
        reader = csv.reader(
            handle, delimiter=_resolve_delimiter(path, delimiter, sample), strict=True
        )

        raw_headers = next(reader, None)
        if raw_headers is None:
            raise ValueError("The metadata file is empty.")
        headers = [h.strip() for h in raw_headers]
        if not any(headers):
            raise ValueError("The header row is blank.")

        rows: list[dict[str, str]] = []
        issues: list[str] = []
        for values in reader:
            line = reader.line_num  # physical line, correct even with quoted newlines
            if not values or all(not v.strip() for v in values):
                issues.append(f"Line {line}: blank row.")
                continue
            if len(values) != len(headers):
                issues.append(f"Line {line}: {len(values)} fields; expected {len(headers)}.")
            row: dict[str, str] = {}
            for index, name in enumerate(headers):
                # setdefault keeps the first occurrence of a duplicated header;
                # duplicates are reported separately and fail validation.
                row.setdefault(name, values[index].strip() if index < len(values) else "")
            rows.append(row)

    return Table(path=path, headers=headers, rows=rows, issues=issues)
