"""Script 02 — GEO → BioSample → SRA Experiment → SRR crosswalk validation.

Input: one row per SRR run with columns
    geo_accession, biosample, sra_experiment, srr_accession   (+ optional condition)

Validates the *accession relationships* only. The number of SRR runs per
experiment is an observation from the current SRA run table — SRR runs are
sequencing runs, never biological replicates. Mouse-to-library mapping and
pre-library pooling are not inferred.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from collections.abc import Sequence
from pathlib import Path

from .gse295590 import GROUPS, REPLICATE_SAFEGUARDS, SAMPLES, SERIES
from .report import Findings, codes, emit, status_line
from .tabular import DELIMITER_CHOICES, Table, cell, read_table

STATUS_LABEL = "GEO-SRA CROSSWALK VALIDATION"
REQUIRED = ("geo_accession", "biosample", "sra_experiment", "srr_accession")
PATTERNS = {
    "geo_accession": re.compile(r"^GSM\d+$", re.IGNORECASE),
    "biosample": re.compile(r"^SAMN\d+$", re.IGNORECASE),
    "sra_experiment": re.compile(r"^SRX\d+$", re.IGNORECASE),
    "srr_accession": re.compile(r"^SRR\d+$", re.IGNORECASE),
}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("metadata", type=Path, help="Run-level crosswalk (.csv, .tsv or .txt)")
    p.add_argument("--report", type=Path, help="Write the Markdown report here (default: stdout)")
    p.add_argument("--delimiter", choices=DELIMITER_CHOICES, default="auto")
    return p


def _safeguards(f: Findings) -> None:
    f.heading("Interpretation safeguards")
    f.info("SRR run counts describe sequencing runs, not biological replicates.")
    for line in REPLICATE_SAFEGUARDS:
        f.info(line)
    f.info(
        "This validation covers accession relationships only; it does not establish the experimental unit."
    )


def validate(table: Table) -> Findings:
    f = Findings()
    headers, rows = table.headers, table.rows

    # --- structure -------------------------------------------------------
    if dups := table.duplicate_headers():
        f.error(f"Duplicate header(s): {codes(dups)}.")
    missing = [c for c in REQUIRED if c not in headers]
    if missing:
        f.error(f"Missing required column(s): {codes(missing)}.")
    for issue in table.issues:
        f.error(issue)
    if not rows:
        f.error("No data rows found.")

    f.info(f"Run-level records: **{len(rows)}**")
    f.info(f"Required columns: {codes(REQUIRED)}")
    if missing:
        _safeguards(f)
        f.heading("Result")
        f.info(f"Validation errors: **{f.errors}**")
        return f

    # --- identifier formats ---------------------------------------------
    for column, pattern in PATTERNS.items():
        values = table.column(column)
        if blanks := [str(i) for i, v in enumerate(values, 1) if not v]:
            f.error(f"Blank `{column}` in data row(s): {', '.join(blanks)}.")
        if bad := sorted({v for v in values if v and not pattern.fullmatch(v)}):
            f.error(f"Invalid `{column}` identifier(s): {codes(bad)}.")

    norm = [
        {c: cell(r, c).upper() for c in REQUIRED} | {"condition": cell(r, "condition").casefold()}
        for r in rows
    ]

    if dup_srr := sorted(
        x
        for x, n in Counter(r["srr_accession"] for r in norm if r["srr_accession"]).items()
        if n > 1
    ):
        f.error(f"Duplicate SRR run ID(s): {codes(dup_srr)}.")

    observed_geo = {r["geo_accession"] for r in norm if r["geo_accession"]}
    expected_geo = {s.geo for s in SAMPLES}
    if absent := sorted(expected_geo - observed_geo):
        f.error(f"Expected GEO accession(s) missing: {codes(absent)}.")
    if extra := sorted(observed_geo - expected_geo):
        f.error(f"Unexpected GEO accession(s): {codes(extra)}.")

    invalid = sorted(
        {r["condition"] for r in norm if r["condition"] and r["condition"] not in GROUPS}
    )
    if "condition" in headers and invalid:
        f.error(
            f"Unrecognised condition label(s): {codes(invalid)}. Allowed: `primary`, `secondary`."
        )

    # --- per deposited sample -------------------------------------------
    f.heading("Per-sample crosswalk")
    for s in SAMPLES:
        group = [r for r in norm if r["geo_accession"] == s.geo]
        if not group:
            continue
        biosamples = {r["biosample"] for r in group if r["biosample"]}
        experiments = {r["sra_experiment"] for r in group if r["sra_experiment"]}
        runs = {r["srr_accession"] for r in group if r["srr_accession"]}
        expected_runs = set(s.srr_runs)

        if biosamples != {s.biosample}:
            f.error(
                f"{s.geo} BioSample mapping is {codes(sorted(biosamples))}; expected `{s.biosample}`."
            )
        if experiments != {s.sra_experiment}:
            f.error(
                f"{s.geo} SRA Experiment mapping is {codes(sorted(experiments))}; "
                f"expected `{s.sra_experiment}`."
            )
        conditions = {r["condition"] for r in group if r["condition"]}
        if conditions and conditions != {s.group}:
            f.error(f"{s.geo} condition is {codes(sorted(conditions))}; expected `{s.group}`.")
        if missing_runs := sorted(expected_runs - runs):
            f.error(
                f"{s.geo}: SRR run(s) listed in the current SRA run table are missing: {codes(missing_runs)}."
            )
        if extra_runs := sorted(runs - expected_runs):
            f.error(f"{s.geo}: SRR run(s) not in the current SRA run table: {codes(extra_runs)}.")

        f.info(
            f"`{s.geo}` → `{s.biosample}` → `{s.sra_experiment}`: {len(runs)} distinct SRR run(s) "
            f"observed; the current SRA run table lists {len(expected_runs)} for this experiment "
            f"(an observation about sequencing runs, not a replicate count)."
        )

    # --- one-to-one checks ----------------------------------------------
    for label, key in (("BioSample", "biosample"), ("SRA Experiment", "sra_experiment")):
        linked: dict[str, set[str]] = defaultdict(set)
        for r in norm:
            if r[key] and r["geo_accession"]:
                linked[r[key]].add(r["geo_accession"])
        if shared := {k: v for k, v in linked.items() if len(v) > 1}:
            f.error(
                f"{label} ID linked to multiple GEO samples: "
                + "; ".join(f"`{k}` → {', '.join(sorted(v))}" for k, v in sorted(shared.items()))
                + "."
            )

    runs_by_exp: dict[str, set[str]] = defaultdict(set)
    for r in norm:
        if r["sra_experiment"] and r["srr_accession"]:
            runs_by_exp[r["sra_experiment"]].add(r["srr_accession"])
    f.heading("SRR runs by SRA Experiment")
    for exp in sorted(runs_by_exp):
        f.info(
            f"`{exp}`: {len(runs_by_exp[exp])} distinct SRR run(s) — {codes(sorted(runs_by_exp[exp]))}."
        )

    _safeguards(f)
    f.heading("Result")
    f.info(f"Validation errors: **{f.errors}**")
    return f


def render(source: Path, f: Findings) -> str:
    return "\n".join(
        [
            f"# GEO–SRA crosswalk validation — {SERIES}",
            "",
            f"- Input: `{source}`",
            f"- **{status_line(STATUS_LABEL, f.passed)}**",
            "",
            "> **Scope.** This status covers the accession relationships "
            "GEO → BioSample → SRA Experiment → SRR only. It does **not** validate the underlying "
            "biological replicate structure.",
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
    f = validate(table)
    return emit(
        render(args.metadata, f), args.report, status_line(STATUS_LABEL, f.passed), f.passed
    )


if __name__ == "__main__":
    raise SystemExit(main())
