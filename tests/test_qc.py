"""Phase 3 QC on a synthetic dataset shaped like GSE295590.

S2 is deliberately small and S3 deliberately mitochondria-rich, so the report
must flag both — without removing any cells in report-only mode.
"""

import gzip

import numpy as np
import pytest

pytest.importorskip("scanpy")
pytest.importorskip("yaml")

from ilc2tri import qc

GENES = (
    [f"mt-Gene{i}" for i in range(13)]
    + [f"Rpl{i}" for i in range(10)]
    + [f"Gene{i}" for i in range(377)]
)
SPEC = {  # geo: (title, n_cells, mito boost)
    "GSM8953514": ("Primary-Inf-ILCs_S1", 600, 1.0),
    "GSM8953515": ("Primary-Inf-ILCs_S2", 60, 1.0),
    "GSM8953516": ("Secondary-Inf-ILCs_S3", 580, 6.0),
}


def write_10x(directory, prefix, n_cells, mito_boost, rng):
    rates = rng.gamma(0.6, 1.0, len(GENES))
    rates[:13] *= mito_boost
    counts = rng.poisson(rates * rng.uniform(1.5, 4.0, (n_cells, 1)))
    rows, cols = np.nonzero(counts.T)  # genes x cells
    with gzip.open(directory / f"{prefix}matrix.mtx.gz", "wt") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n%\n")
        fh.write(f"{len(GENES)} {n_cells} {len(rows)}\n")
        for r, c in zip(rows, cols, strict=True):
            fh.write(f"{r + 1} {c + 1} {counts[c, r]}\n")
    with gzip.open(directory / f"{prefix}barcodes.tsv.gz", "wt") as fh:
        fh.writelines(f"CELL{i:05d}-1\n" for i in range(n_cells))
    with gzip.open(directory / f"{prefix}features.tsv.gz", "wt") as fh:
        fh.writelines(f"ENSMUSG{i:011d}\t{g}\tGene Expression\n" for i, g in enumerate(GENES))


@pytest.fixture
def raw_dir(tmp_path):
    rng = np.random.default_rng(1)
    d = tmp_path / "raw"
    d.mkdir()
    for geo, (title, n, boost) in SPEC.items():
        write_10x(d, f"{geo}_{title}_", n, boost, rng)
    return d


@pytest.fixture
def params(tmp_path):
    def _write(**overrides):
        t = {
            "min_cells_per_sample": None,
            "max_pct_mito": None,
            "min_genes_per_cell": None,
            "max_genes_per_cell": None,
        }
        t.update(overrides)
        lines = ["random_seed: 0", "phase3_qc:"] + [
            f"  {k}: {'null' if v is None else v}" for k, v in t.items()
        ]
        path = tmp_path / "parameters.yaml"
        path.write_text("\n".join(lines) + "\n")
        return path

    return _write


def test_discover_matches_geo_prefixes(raw_dir):
    found = qc.discover(raw_dir)
    assert [s.geo for s in found] == list(SPEC)


def test_discover_reports_missing_companion(raw_dir):
    (raw_dir / "GSM8953515_Primary-Inf-ILCs_S2_features.tsv.gz").unlink()
    with pytest.raises(FileNotFoundError):
        qc.discover(raw_dir)


def test_report_only_flags_without_removing(raw_dir, params, tmp_path):
    out, report = tmp_path / "qc", tmp_path / "reports" / "qc.md"
    rc = qc.main(
        [
            "--raw-dir",
            str(raw_dir),
            "--out-dir",
            str(out),
            "--report",
            str(report),
            "--params",
            str(params()),
            "--no-doublets",
        ]
    )
    assert rc == 0
    text = report.read_text()
    assert "NOT EVALUABLE" in text
    assert "low cell recovery" in text and "GSM8953515" in text
    assert "elevated mitochondrial fraction" in text and "GSM8953516" in text

    import anndata as ad
    import pandas as pd

    summary = pd.read_csv(out / "qc_summary.csv", index_col=0)
    assert summary["n_cells"].to_dict() == {g: n for g, (_, n, _) in SPEC.items()}
    adata = ad.read_h5ad(out / "GSE295590_qc_flagged_unfiltered.h5ad")
    assert adata.n_obs == sum(n for _, n, _ in SPEC.values())  # nothing removed
    assert {"outlier_counts", "outlier_mito", "pct_counts_mt"} <= set(adata.obs.columns)
    assert (out / "qc_figures" / "qc_violins.png").stat().st_size > 0


def test_filter_mode_applies_thresholds(raw_dir, params, tmp_path):
    out, report = tmp_path / "qc", tmp_path / "qc.md"
    p = params(
        min_cells_per_sample=100, max_pct_mito=20, min_genes_per_cell=50, max_genes_per_cell=400
    )
    assert (
        qc.main(
            [
                "--raw-dir",
                str(raw_dir),
                "--out-dir",
                str(out),
                "--report",
                str(report),
                "--params",
                str(p),
                "--no-doublets",
            ]
        )
        == 0
    )
    text = report.read_text()
    assert "PHASE 3 GATE: FAIL" in text  # S2 has 60 cells < 100

    import anndata as ad

    filtered = ad.read_h5ad(out / "GSE295590_qc_filtered.h5ad")
    assert filtered.obs["pct_counts_mt"].max() <= 20


def test_missing_raw_dir_exits_2(tmp_path, params):
    assert qc.main(["--raw-dir", str(tmp_path / "nope"), "--params", str(params())]) == 2
