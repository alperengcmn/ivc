"""Script 03 — Phase 3 single-cell quality control for GSE295590.

Loads the per-sample 10x matrices, computes per-cell QC metrics, scores
doublets, and writes a QC summary, figures and a Markdown report.

Two modes, chosen by config/parameters.yaml → phase3_qc:

* **report-only** (any threshold is `null`, the default): no cell is removed.
  Cells are *flagged* with data-driven MAD outlier rules so thresholds can be
  agreed with the data in view. The Phase 3 gate is reported as NOT EVALUABLE.
* **filter** (all thresholds set): thresholds are applied, a filtered object
  is written, and the gate is evaluated.

Low cell recovery (GSM8953515) and elevated mitochondrial fraction
(GSM8953516) are flagged for investigation, never excluded automatically.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .gse295590 import BY_GEO, SERIES

THRESHOLD_KEYS = (
    "min_cells_per_sample",
    "max_pct_mito",
    "min_genes_per_cell",
    "max_genes_per_cell",
)
MT_PREFIX = "mt-"  # mouse
RIBO_PATTERN = r"^Rp[sl]\d"
# Outlier rules used only for *suggestions* in report-only mode (Germain et al. 2020; Heumos et al. 2023).
MAD_COUNTS = 5
MAD_MITO = 3
UNFILTERED_BARCODE_LIMIT = 50_000  # more barcodes than this suggests a raw (unfiltered) matrix
LOW_RECOVERY_FRACTION = 0.25  # flag samples with < 25 % of the median cell count
GSM_RE = re.compile(r"(GSM\d+)")
GROUP_COLORS = {
    "primary": "#2563eb",
    "secondary": "#d97706",
}  # blue / amber, colour-blind safe pair
SCRUBLET_MIN_CELLS = 500  # scrublet estimates are unreliable on very small libraries
SCRUBLET_SILENT_MIN_CELLS = (
    2_000  # 0 % predicted doublets in a library this large = auto-threshold failure
)
DEPTH_RATIO_FLAG = 3.0  # flag libraries sequenced > 3x shallower than the deepest one
HIGH_MITO_PCT = 20.0  # descriptive cut-off for the high-mito tail
HIGH_MITO_TAIL_FLAG = 5.0  # flag samples where > 5 % of cells exceed HIGH_MITO_PCT
COUNT_FLOOR_SPREAD = 1.2  # all per-sample minima within 20 % of each other -> shared pre-filter


@dataclass
class SampleFiles:
    geo: str
    prefix: str
    directory: Path


# ---------------------------------------------------------------- discovery


def discover(raw_dir: Path) -> list[SampleFiles]:
    """Find 10x triplets (<prefix>matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)."""
    found: list[SampleFiles] = []
    for matrix in sorted(raw_dir.rglob("*matrix.mtx.gz")):
        prefix = matrix.name[: -len("matrix.mtx.gz")]
        missing = [
            n
            for n in ("barcodes.tsv.gz", "features.tsv.gz")
            if not (matrix.parent / f"{prefix}{n}").is_file()
        ]
        if missing:
            raise FileNotFoundError(
                f"{matrix}: companion file(s) missing: {', '.join(prefix + m for m in missing)}"
            )
        match = GSM_RE.search(prefix) or GSM_RE.search(str(matrix.parent))
        if not match:
            raise ValueError(f"Cannot find a GSM accession in {matrix}")
        found.append(SampleFiles(match.group(1).upper(), prefix, matrix.parent))
    geos = [s.geo for s in found]
    if dups := sorted({g for g in geos if geos.count(g) > 1}):
        raise ValueError(f"More than one matrix per sample: {', '.join(dups)}")
    return found


# ---------------------------------------------------------------- core


def load(samples: Sequence[SampleFiles]):
    import anndata as ad
    import scanpy as sc

    adatas = {}
    for s in samples:
        a = sc.read_10x_mtx(
            s.directory, prefix=s.prefix, var_names="gene_symbols", make_unique=True
        )
        meta = BY_GEO.get(s.geo)
        a.obs["sample"] = s.geo
        a.obs["sample_title"] = meta.title if meta else "unknown"
        a.obs["infection_group"] = meta.group if meta else "unknown"
        adatas[s.geo] = a
    adata = ad.concat(adatas, label="sample_key", index_unique="-", merge="same")
    adata.obs = adata.obs.drop(columns="sample_key")
    restore_categories(adata)
    return adata


CATEGORICAL_OBS = ("sample", "sample_title", "infection_group")


def restore_categories(adata) -> None:
    """(Re)apply categorical dtypes. Some scanpy versions drop them inside scrublet."""
    for col in CATEGORICAL_OBS:
        if col in adata.obs:
            adata.obs[col] = adata.obs[col].astype(str).astype("category")


def sample_order(obs) -> list[str]:
    """Sample IDs in a stable order, whatever the column dtype."""
    col = obs["sample"]
    if hasattr(col, "cat"):
        return [str(c) for c in col.cat.categories if (col == c).any()]
    return sorted(col.astype(str).unique())


def compute_metrics(adata, seed: int = 0, doublets: bool = True) -> str:
    """Add QC metrics to adata.obs. Returns a doublet-scoring status message."""
    import scanpy as sc

    adata.var["mt"] = adata.var_names.str.startswith(MT_PREFIX)
    adata.var["ribo"] = adata.var_names.str.match(RIBO_PATTERN)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=True, inplace=True
    )

    if not doublets:
        return "skipped (--no-doublets)"
    cells_before = adata.obs_names.copy()
    try:
        sc.pp.scrublet(adata, batch_key="sample", random_state=seed)
    except Exception as exc:  # report the failure, never hide it
        return f"FAILED — {type(exc).__name__}: {exc}"
    finally:
        restore_categories(adata)
    if not adata.obs_names.equals(cells_before):
        # Never silently misattach metrics to the wrong cells.
        raise RuntimeError("scrublet changed the cell order; refusing to continue.")
    return "scrublet, per sample"


def mad_outlier(values, n_mads: float, upper_only: bool = False):
    import numpy as np

    med = np.median(values)
    mad = np.median(np.abs(values - med))
    high = values > med + n_mads * mad
    return high if upper_only else high | (values < med - n_mads * mad)


def flag_outliers(adata) -> None:
    """Per-sample MAD flags. Flags only — nothing is removed here."""
    import numpy as np

    flags = {k: np.zeros(adata.n_obs, dtype=bool) for k in ("outlier_counts", "outlier_mito")}
    for sample in sample_order(adata.obs):
        idx = np.flatnonzero((adata.obs["sample"] == sample).to_numpy())
        obs = adata.obs.iloc[idx]
        flags["outlier_counts"][idx] = mad_outlier(
            obs["log1p_total_counts"].to_numpy(), MAD_COUNTS
        ) | mad_outlier(obs["log1p_n_genes_by_counts"].to_numpy(), MAD_COUNTS)
        flags["outlier_mito"][idx] = mad_outlier(
            obs["pct_counts_mt"].to_numpy(), MAD_MITO, upper_only=True
        )
    for k, v in flags.items():
        adata.obs[k] = v


def summarise(adata):
    import pandas as pd

    obs = adata.obs
    g = obs.groupby("sample", observed=True)
    summary = pd.DataFrame(
        {
            "sample_title": g["sample_title"].first(),
            "infection_group": g["infection_group"].first(),
            "n_cells": g.size(),
            "min_total_counts": g["total_counts"].min(),
            "median_total_counts": g["total_counts"].median(),
            "median_genes": g["n_genes_by_counts"].median(),
            "median_pct_mt": g["pct_counts_mt"].median(),
            "p90_pct_mt": g["pct_counts_mt"].quantile(0.9),
            "pct_cells_mt_gt20": g["pct_counts_mt"].apply(
                lambda s: (s > HIGH_MITO_PCT).mean() * 100
            ),
            "median_pct_ribo": g["pct_counts_ribo"].median(),
            "pct_flag_counts": g["outlier_counts"].mean() * 100,
            "pct_flag_mito": g["outlier_mito"].mean() * 100,
        }
    )
    if "predicted_doublet" in obs:
        summary["pct_predicted_doublet"] = g["predicted_doublet"].mean() * 100
    # Suggested mito ceiling per sample = median + 3 MAD (for discussion only).
    summary["suggested_max_pct_mt"] = g["pct_counts_mt"].apply(
        lambda s: s.median() + MAD_MITO * (s - s.median()).abs().median()
    )
    median_cells = summary["n_cells"].median()
    summary["cells_vs_median"] = summary["n_cells"] / median_cells
    return summary.round(3)


def dataset_notes(summary) -> list[str]:
    """Observations that concern the dataset as a whole rather than one sample."""
    notes = []
    floors = summary["min_total_counts"]
    if len(floors) > 1 and floors.min() > 0 and floors.max() / floors.min() <= COUNT_FLOOR_SPREAD:
        notes.append(
            f"Every library starts at ≈ {floors.min():,.0f} UMI per cell — the deposited matrices appear to be "
            "**already filtered** (e.g. a ≥ 500-UMI cut-off by the authors). Barcodes below that floor cannot be "
            "recovered from these files."
        )
    depth = summary["median_total_counts"]
    if len(depth) > 1 and depth.min() > 0 and depth.max() / depth.min() > DEPTH_RATIO_FLAG:
        shallow = depth.idxmin()
        notes.append(
            f"**Sequencing depth differs {depth.max() / depth.min():.1f}×** between libraries "
            f"(`{shallow}` median {depth.min():,.0f} UMI vs {depth.max():,.0f}). Depth is a technical confounder of "
            "any between-sample comparison: normalisation (Phase 4) must account for it, and Phase 11 should include "
            "a depth-matched (downsampled) sensitivity analysis."
        )
    return notes


def sample_notes(summary) -> list[str]:
    notes = []
    mt_median_all = summary["median_pct_mt"].median()
    for geo, row in summary.iterrows():
        if row["n_cells"] > UNFILTERED_BARCODE_LIMIT:
            notes.append(
                f"`{geo}`: {int(row['n_cells']):,} barcodes — the matrix may be **unfiltered** (raw droplets)."
            )
        if row["cells_vs_median"] < LOW_RECOVERY_FRACTION:
            notes.append(
                f"`{geo}`: low cell recovery ({int(row['n_cells']):,} cells, {row['cells_vs_median']:.0%} of the median). "
                "Investigate library recovery / input; do **not** use it to infer mouse numbers."
            )
        if "pct_predicted_doublet" in row and row["n_cells"] < SCRUBLET_MIN_CELLS:
            notes.append(
                f"`{geo}`: doublet rate ({row['pct_predicted_doublet']:.1f} %) is unreliable with "
                f"< {SCRUBLET_MIN_CELLS} cells; do not filter this sample on predicted doublets alone."
            )
        if (
            "pct_predicted_doublet" in row
            and row["n_cells"] >= SCRUBLET_SILENT_MIN_CELLS
            and row["pct_predicted_doublet"] == 0
        ):
            notes.append(
                f"`{geo}`: scrublet predicted **0 %** doublets in {int(row['n_cells']):,} cells — its automatic threshold "
                "most likely failed (10x libraries of this size typically carry several % doublets). Keep "
                "`doublet_score` and assess doublets at cluster level in Phase 5 instead."
            )
        if mt_median_all > 0 and row["median_pct_mt"] > 1.5 * mt_median_all:
            notes.append(
                f"`{geo}`: elevated mitochondrial fraction (median {row['median_pct_mt']:.2f} % vs "
                f"{mt_median_all:.2f} % across samples). Investigate before excluding."
            )
        if row["pct_cells_mt_gt20"] > HIGH_MITO_TAIL_FLAG:
            notes.append(
                f"`{geo}`: {row['pct_cells_mt_gt20']:.1f} % of cells have > {HIGH_MITO_PCT:.0f} % mitochondrial "
                "reads — a distinct high-mito / low-complexity tail typical of damaged cells. Inspect "
                "`qc_counts_vs_genes.png`."
            )
    return notes


def apply_thresholds(adata, t: dict[str, Any]):
    keep = (
        (adata.obs["pct_counts_mt"] <= t["max_pct_mito"])
        & (adata.obs["n_genes_by_counts"] >= t["min_genes_per_cell"])
        & (adata.obs["n_genes_by_counts"] <= t["max_genes_per_cell"])
    )
    if t.get("remove_predicted_doublets") and "predicted_doublet" in adata.obs:
        keep &= ~adata.obs["predicted_doublet"].astype(bool)
    return adata[keep.to_numpy()].copy()


# ---------------------------------------------------------------- figures


def figures(adata, out_dir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    obs = adata.obs
    samples = sample_order(obs)
    paths = []

    metrics = [
        ("total_counts", "UMI counts per cell", True),
        ("n_genes_by_counts", "Genes per cell", True),
        ("pct_counts_mt", "Mitochondrial %", False),
        ("pct_counts_ribo", "Ribosomal %", False),
    ]
    fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 4))
    for ax, (col, label, log) in zip(axes, metrics, strict=True):
        data = [obs.loc[obs["sample"] == s, col].to_numpy() for s in samples]
        parts = ax.violinplot(data, showmedians=True)
        for body, s in zip(parts["bodies"], samples, strict=True):
            group = obs.loc[obs["sample"] == s, "infection_group"].iloc[0]
            body.set_facecolor(GROUP_COLORS.get(group, "#888888"))
            body.set_alpha(0.6)
        ax.set_xticks(
            range(1, len(samples) + 1),
            [f"{s}\n({obs.loc[obs['sample'] == s, 'infection_group'].iloc[0]})" for s in samples],
            fontsize=8,
        )
        ax.set_title(label)
        if log:
            ax.set_yscale("log")
    fig.tight_layout()
    paths.append(out_dir / "qc_violins.png")
    fig.savefig(paths[-1], dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, len(samples), figsize=(4.5 * len(samples), 4), squeeze=False)
    for ax, s in zip(axes[0], samples, strict=True):
        sub = obs[obs["sample"] == s]
        sca = ax.scatter(
            sub["total_counts"],
            sub["n_genes_by_counts"],
            c=sub["pct_counts_mt"],
            s=3,
            cmap="viridis",
            rasterized=True,
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"{s} (n = {len(sub):,})")
        ax.set_xlabel("UMI counts")
        ax.set_ylabel("Genes")
        fig.colorbar(sca, ax=ax, label="Mito %")
    fig.tight_layout()
    paths.append(out_dir / "qc_counts_vs_genes.png")
    fig.savefig(paths[-1], dpi=150)
    plt.close(fig)
    return paths


# ---------------------------------------------------------------- report


def _table(df) -> str:
    cols = ["sample", *df.columns]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for idx, row in df.iterrows():
        cells = [str(idx)] + [
            f"{v:,.2f}" if isinstance(v, float) else f"{v:,}" if isinstance(v, int) else str(v)
            for v in row
        ]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def render(
    summary,
    mode: str,
    thresholds: dict[str, Any],
    doublet_status: str,
    notes: list[str],
    figure_paths: list[Path],
    retained=None,
    global_notes: list[str] | None = None,
) -> str:
    unset = [k for k in THRESHOLD_KEYS if thresholds.get(k) is None]
    if mode == "report-only":
        gate = (
            "**PHASE 3 GATE: NOT EVALUABLE** — threshold(s) not yet agreed: "
            + ", ".join(f"`{k}`" for k in unset)
            + ". No cells were removed."
        )
    else:
        floor = thresholds["min_cells_per_sample"]
        ok = bool((retained >= floor).all())
        gate = (
            f"**PHASE 3 GATE: {'PASS' if ok else 'FAIL'}** — {int(retained.sum()):,} cells retained; "
            f"every sample needs ≥ {floor:,} cells after filtering."
        )
    t = summary[
        [
            "sample_title",
            "infection_group",
            "n_cells",
            "cells_vs_median",
            "median_total_counts",
            "median_genes",
            "median_pct_mt",
            "p90_pct_mt",
            "pct_cells_mt_gt20",
            "median_pct_ribo",
        ]
    ]
    flags = summary[
        [
            c
            for c in (
                "pct_flag_counts",
                "pct_flag_mito",
                "pct_predicted_doublet",
                "suggested_max_pct_mt",
            )
            if c in summary
        ]
    ]
    return "\n".join(
        [
            f"# Phase 3 — single-cell QC · {SERIES}",
            "",
            f"- Mode: **{mode}**",
            f"- {gate}",
            f"- Doublet scoring: {doublet_status}",
            "",
            "> Cell counts are not biological sample sizes. These metrics describe libraries; the replicate "
            "structure of GSE295590 remains UNRESOLVED (docs/GSE295590_design.md).",
            "",
            "## Per-sample metrics",
            _table(t),
            "",
            "## Flags (descriptive — no cells removed in report-only mode)",
            f"Counts/genes outliers: > {MAD_COUNTS} MAD on log scale, per sample. Mito outliers: > {MAD_MITO} MAD above "
            "the sample median. `suggested_max_pct_mt` = median + 3 MAD — a starting point for the threshold discussion, "
            "not an agreed value. On narrow mito distributions MAD-based cut-offs are very stringent and can remove "
            "metabolically active (OXPHOS-high) cells, so treat them as an upper bound on stringency.",
            "",
            _table(flags),
            "",
            *(
                ["## Retained after filtering", _table(retained.to_frame("retained_cells")), ""]
                if retained is not None
                else []
            ),
            "## Dataset-level observations",
            *([f"- {n}" for n in global_notes or []] or ["- None."]),
            "",
            "## Samples to investigate",
            *([f"- {n}" for n in notes] or ["- None flagged."]),
            "",
            "## Figures",
            *[f"- `{p}`" for p in figure_paths],
            "",
            "## Next step",
            "Agree `phase3_qc` thresholds in `config/parameters.yaml` (log the decision in `docs/decision_log.md`), "
            "then re-run `make qc` to apply them and evaluate the gate.",
            "",
        ]
    )


# ---------------------------------------------------------------- cli


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--raw-dir", type=Path, default=Path("data/raw") / SERIES)
    p.add_argument("--out-dir", type=Path, default=Path("results/qc"))
    p.add_argument("--report", type=Path, default=Path("reports/qc/03_qc_report.md"))
    p.add_argument(
        "--params",
        type=Path,
        default=None,
        help="parameters.yaml (default: config/parameters.yaml)",
    )
    p.add_argument("--no-doublets", action="store_true", help="skip scrublet doublet scoring")
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        import scanpy  # noqa: F401
    except ImportError:
        print(
            "Error: Script 03 needs the analysis environment (conda env create -f environment.yml).",
            file=sys.stderr,
        )
        return 2

    from .config import parameters

    params = parameters(args.params)
    seed = int(params.get("random_seed", 0))
    thresholds = params.get("phase3_qc", {}) or {}
    mode = "filter" if all(thresholds.get(k) is not None for k in THRESHOLD_KEYS) else "report-only"

    if not args.raw_dir.is_dir():
        print(f"Error: {args.raw_dir} not found. Run `make fetch` first.", file=sys.stderr)
        return 2
    try:
        samples = discover(args.raw_dir)
    except (FileNotFoundError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    if not samples:
        print(f"Error: no *matrix.mtx.gz files under {args.raw_dir}.", file=sys.stderr)
        return 2
    print(f"[qc] {len(samples)} sample(s): {', '.join(s.geo for s in samples)} · mode: {mode}")

    adata = load(samples)
    doublet_status = compute_metrics(adata, seed=seed, doublets=not args.no_doublets)
    flag_outliers(adata)
    summary = summarise(adata)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.out_dir / "qc_summary.csv")
    adata.obs.to_csv(args.out_dir / "cell_qc_metrics.csv.gz")
    fig_paths = figures(adata, args.out_dir / "qc_figures")

    retained = None
    if mode == "filter":
        filtered = apply_thresholds(adata, thresholds)
        retained = filtered.obs["sample"].value_counts().reindex(summary.index, fill_value=0)
        filtered.write_h5ad(args.out_dir / f"{SERIES}_qc_filtered.h5ad", compression="gzip")
    else:
        adata.write_h5ad(args.out_dir / f"{SERIES}_qc_flagged_unfiltered.h5ad", compression="gzip")

    args.report.parent.mkdir(parents=True, exist_ok=True)
    report = render(
        summary,
        mode,
        thresholds,
        doublet_status,
        sample_notes(summary),
        fig_paths,
        retained,
        dataset_notes(summary),
    )
    args.report.write_text(report, encoding="utf-8")
    print(f"[qc] report → {args.report}")
    print(f"[qc] outputs → {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
