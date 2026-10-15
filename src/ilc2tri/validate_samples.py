"""Script 01 — structural validation of GSE295590 sample-level metadata.

Checks that a sample table (one row per GEO sample) is well formed and agrees
with the three deposited GEO samples. It never infers mouse identity, pooling
or biological replicates: those stay UNRESOLVED unless supported by evidence.

Required columns (configurable): geo_accession, infection_group.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from collections.abc import Sequence
from pathlib import Path

from .gse295590 import BY_GEO, GROUPS, REPLICATE_SAFEGUARDS, SERIES
from .report import Findings, codes, emit, status_line
from .tabular import DELIMITER_CHOICES, Table, read_table

STATUS_LABEL = "STRUCTURAL METADATA VALIDATION"
DEFAULT_REPLICATE_COLUMNS = ("mouse_id", "biological_replicate")
GSM = re.compile(r"^GSM\d+$", re.IGNORECASE)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("metadata", type=Path, help="Sample table (.csv, .tsv or .txt)")
    p.add_argument("--report", type=Path, help="Write the Markdown report here (default: stdout)")
    p.add_argument("--geo-column", default="geo_accession")
    p.add_argument("--group-column", default="infection_group")
    p.add_argument(
        "--replicate-columns",
        nargs="*",
        default=list(DEFAULT_REPLICATE_COLUMNS),
        help="Columns inspected for replicate metadata (default: mouse_id biological_replicate)",
    )
    p.add_argument("--delimiter", choices=DELIMITER_CHOICES, default="auto")
    return p


def validate(
    table: Table,
    geo_column: str = "geo_accession",
    group_column: str = "infection_group",
    replicate_columns: Sequence[str] = DEFAULT_REPLICATE_COLUMNS,
) -> Findings:
    f = Findings()
    headers, rows = table.headers, table.rows

    # --- structure -------------------------------------------------------
    if dups := table.duplicate_headers():
        f.error(f"Duplicate header(s): {codes(dups)}.")
    if missing := [c for c in (geo_column, group_column) if c not in headers]:
        f.error(f"Required column(s) missing: {codes(missing)}.")
    for issue in table.issues:
        f.error(issue)
    if not rows:
        f.error("No metadata records were found.")

    f.info(f"Records read: **{len(rows)}**")
    f.info(f"Columns: {codes(headers)}")
    f.info(f"Required columns: {codes([geo_column, group_column])}")

    # --- GEO identifiers -------------------------------------------------
    if geo_column in headers:
        ids = [x.upper() for x in table.column(geo_column)]
        if blanks := [str(i) for i, x in enumerate(ids, 1) if not x]:
            f.error(f"Missing GEO sample identifier in data row(s): {', '.join(blanks)}.")
        if dup_ids := sorted(x for x, n in Counter(x for x in ids if x).items() if n > 1):
            f.error(f"Duplicate GEO sample identifier(s): {codes(dup_ids)}.")
        if bad := sorted({x for x in ids if x and not GSM.fullmatch(x)}):
            f.error(f"Identifier(s) not in GEO `GSM` + digits format: {codes(bad)}.")

        observed = {x for x in ids if x}
        absent = [g for g in BY_GEO if g not in observed]
        extra = sorted(observed - BY_GEO.keys())
        if absent:
            f.note(f"{SERIES} sample(s) absent from this file: {codes(absent)}.")
        if extra:
            f.note(f"Sample(s) outside the {SERIES} sample list: {codes(extra)}.")
        if not absent and not extra:
            f.info(f"GEO sample set matches the three deposited {SERIES} samples.")

    # --- infection groups ------------------------------------------------
    if group_column in headers:
        groups = [g.casefold() for g in table.column(group_column)]
        if blanks := [str(i) for i, g in enumerate(groups, 1) if not g]:
            f.error(f"Missing infection-group label in data row(s): {', '.join(blanks)}.")
        if invalid := sorted({g for g in groups if g and g not in GROUPS}):
            f.error(
                f"Unrecognised infection-group label(s): {codes(invalid)}. "
                f"Allowed: `primary`, `secondary` (case-insensitive)."
            )
        conflicts = []
        if geo_column in headers:
            for geo, group in zip(
                (x.upper() for x in table.column(geo_column)), groups, strict=True
            ):
                expected = BY_GEO.get(geo)
                if expected and group in GROUPS and group != expected.group:
                    conflicts.append(
                        f"{geo}: supplied `{group}`, deposited title indicates `{expected.group}`"
                    )
        if conflicts:
            f.error(
                "Infection group conflicts with the deposited sample: " + "; ".join(conflicts) + "."
            )
        if not blanks and not invalid and not conflicts:
            f.info(
                "Infection-group labels: populated, valid, and consistent with deposited sample titles."
            )

    # --- completeness ----------------------------------------------------
    gaps = [(h, sum(not v for v in table.column(h))) for h in headers if h]
    if gaps := [(h, n) for h, n in gaps if n]:
        f.info("Blank cells by column: " + "; ".join(f"`{h}`: {n}" for h, n in gaps) + ".")
    else:
        f.info("Blank cells: none in the supplied columns.")

    # --- replicate metadata (reported, never inferred) --------------------
    f.heading("Replicate metadata")
    wanted = list(dict.fromkeys(replicate_columns))
    present = [c for c in wanted if c in headers]
    if not present:
        f.unresolved(f"No replicate-identifying column found (checked: {codes(wanted) or 'none'}).")
    filled_any = False
    for c in present:
        filled = sum(bool(v) for v in table.column(c))
        filled_any = filled_any or filled > 0
        f.info(f"Replicate field `{c}`: {filled}/{len(rows)} record(s) populated.")
    if present and not filled_any:
        f.unresolved(
            "Replicate fields exist but are empty: no replicate information has been recorded yet."
        )
    elif filled_any:
        f.unresolved(
            "Populated replicate fields do not by themselves verify independent biological "
            "replicates or the mouse-to-library mapping. Check against the original study "
            "records; no assignments were inferred."
        )
    for line in REPLICATE_SAFEGUARDS:
        f.info(line)

    f.heading("Result")
    f.info(f"Validation errors: **{f.errors}**")
    return f


def render(source: Path, f: Findings) -> str:
    return "\n".join(
        [
            f"# Sample metadata validation — {SERIES}",
            "",
            f"- Input: `{source}`",
            f"- **{status_line(STATUS_LABEL, f.passed)}**",
            "",
            "> **Scope.** This status covers the *structure* of the sample metadata only "
            "(identifiers, group labels, completeness). It does **not** establish the biological "
            "replicate structure, the mouse-to-library mapping, or pre-library pooling.",
            "",
            "## Findings",
            *f.lines,
            "",
        ]
    )


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        table = read_table(args.metadata, args.delimiter)
    except (OSError, UnicodeError, csv.Error, ValueError) as exc:
        print(f"Error: could not read metadata: {exc}", file=sys.stderr)
        return 2
    f = validate(table, args.geo_column, args.group_column, args.replicate_columns)
    return emit(
        render(args.metadata, f), args.report, status_line(STATUS_LABEL, f.passed), f.passed
    )


if __name__ == "__main__":
    raise SystemExit(main())
