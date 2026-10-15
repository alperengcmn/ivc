"""Reference accession structure for GSE295590 — the single source of truth.

Every validator compares supplied metadata against these records. They describe
the *deposited accession structure* (GEO → BioSample → SRA Experiment → SRR) as
it appears in GEO and the SRA run table. They say nothing about which mice
contributed to which library: that mapping is UNRESOLVED (docs/GSE295590_design.md).
"""

from __future__ import annotations

from dataclasses import dataclass

SERIES = "GSE295590"
GROUPS = frozenset({"primary", "secondary"})


@dataclass(frozen=True)
class DepositedSample:
    geo: str
    title: str
    group: str
    biosample: str
    sra_experiment: str
    srr_runs: tuple[str, ...]  # sequencing runs — NOT biological replicates


SAMPLES: tuple[DepositedSample, ...] = (
    DepositedSample(
        "GSM8953514",
        "Primary-Inf-ILCs_S1",
        "primary",
        "SAMN48137559",
        "SRX28549631",
        ("SRR33304601", "SRR33304602"),
    ),
    DepositedSample(
        "GSM8953515",
        "Primary-Inf-ILCs_S2",
        "primary",
        "SAMN48137558",
        "SRX28549632",
        ("SRR33304599", "SRR33304600"),
    ),
    DepositedSample(
        "GSM8953516",
        "Secondary-Inf-ILCs_S3",
        "secondary",
        "SAMN48137557",
        "SRX28549633",
        ("SRR33304597", "SRR33304598"),
    ),
)

BY_GEO: dict[str, DepositedSample] = {s.geo: s for s in SAMPLES}

# Wording reused verbatim by every report so the caveat never drifts.
REPLICATE_SAFEGUARDS: tuple[str, ...] = (
    "Mouse-to-library mapping: **UNKNOWN / not established**.",
    "Pre-library pooling status: **UNKNOWN / not established**.",
    "Independent biological replicate structure: **UNKNOWN / not inferred**.",
)
