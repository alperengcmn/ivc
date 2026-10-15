"""Script 09 — Phase 9: epigenomic integration with GSE172258 ATAC-seq (pre-registered, D-014).

Differential accessibility (pydeseq2) per comparison; peaks are linked to the nearest RefSeq
TSS (mm10) within 50 kb; a gene's accessibility change is the mean Wald statistic of its peaks.
The frozen signature is then tested: do `recall_up` genes gain accessibility, and
`recall_down` genes lose it, relative to all other genes?

Motif enrichment is not performed (needs a genome sequence and motif database; D-014).
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import shutil
import sys
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

MAX_DIST = 50_000
MIN_COUNT = 10  # peak kept if summed counts >= MIN_COUNT * n_samples
P_SUPPORT = 0.05
SIGNATURE = Path("signatures/GSE295590_recall_ILC2_signature_v1.tsv")


@dataclass(frozen=True)
class AtacContrast:
    name: str
    label: str
    file: str
    test: str  # column prefix of the memory/sensitised group
    ref: str  # column prefix of the comparator
    recall: bool  # True = recall comparison (counts towards the verdict)


CONTRASTS = (
    AtacContrast(
        "rag_recall",
        "Rag1⁻/⁻ Alt/Alt vs Sal/Alt",
        "GSE172258_Alam_ILC2_memory_ATAC_counts.txt.gz",
        "Allergen_recall",
        "Saline_recall",
        True,
    ),
    AtacContrast(
        "wt_recall",
        "WT Alt/Alt vs Sal/Alt",
        "GSE172258_WT_Alam_ILC2_memory_ATAC_counts.txt.gz",
        "Allergen_Allergen",
        "Saline_Allergen",
        True,
    ),
    AtacContrast(
        "rag_resting",
        "Rag1⁻/⁻ Alt/– vs Sal/– (resting memory)",
        "GSE172258_Alam_ILC2_memory_ATAC_counts.txt.gz",
        "Allergen_norecall",
        "Saline_norecall",
        False,
    ),
)


def load_tss(refgene: Path) -> dict[str, tuple[list[int], list[str]]]:
    """chrom -> (sorted TSS positions, gene names) from UCSC refGene."""
    per_chrom: dict[str, list[tuple[int, str]]] = defaultdict(list)
    with gzip.open(refgene, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom, strand, start, end, gene = f[2], f[3], int(f[4]), int(f[5]), f[12]
            if "_" in chrom:
                continue
            per_chrom[chrom].append((start if strand == "+" else end, gene))
    out = {}
    for chrom, items in per_chrom.items():
        items = sorted(set(items))
        out[chrom] = ([p for p, _ in items], [g for _, g in items])
    return out


def nearest_gene(peak: str, tss) -> tuple[str | None, int | None]:
    chrom, rng = peak.split(":")
    s, e = (int(x) for x in rng.split("-"))
    mid = (s + e) // 2
    if chrom not in tss:
        return None, None
    pos, genes = tss[chrom]
    i = bisect.bisect_left(pos, mid)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(pos):
            d = abs(pos[j] - mid)
            if best is None or d < best[1]:
                best = (genes[j], d)
    if best is None or best[1] > MAX_DIST:
        return None, None
    return best


def differential(counts, c: AtacContrast, seed: int = 0):
    import pandas as pd
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    cols = [x for x in counts.columns if x.startswith((c.test + "_", c.ref + "_"))]
    sub = counts[cols]
    sub = sub[sub.sum(axis=1) >= MIN_COUNT * len(cols)]
    meta = pd.DataFrame(
        {"group": ["test" if x.startswith(c.test + "_") else "ref" for x in cols]}, index=cols
    )
    dds = DeseqDataSet(
        counts=sub.T.round().astype(int),
        metadata=meta,
        design="~group",
        inference=DefaultInference(n_cpus=4),
        quiet=True,
    )
    dds.deseq2()
    st = DeseqStats(
        dds, contrast=["group", "test", "ref"], inference=DefaultInference(n_cpus=4), quiet=True
    )
    st.summary()
    res = st.results_df.copy()
    res.attrs["n_test"] = int((meta["group"] == "test").sum())
    res.attrs["n_ref"] = int((meta["group"] == "ref").sum())
    return res


def gene_level(res, tss):

    genes, dists = zip(*(nearest_gene(p, tss) for p in res.index), strict=True)
    res = res.assign(gene=genes, distance=dists).dropna(subset=["gene", "stat"])
    g = res.groupby("gene").agg(
        mean_stat=("stat", "mean"),
        n_peaks=("stat", "size"),
        n_up=("padj", lambda s: 0),
        max_abs_lfc=("log2FoldChange", lambda s: s.abs().max()),
    )
    g["n_up"] = (
        res[(res["padj"] < 0.05) & (res["log2FoldChange"] > 0)]
        .groupby("gene")
        .size()
        .reindex(g.index)
        .fillna(0)
    )
    g["n_down"] = (
        res[(res["padj"] < 0.05) & (res["log2FoldChange"] < 0)]
        .groupby("gene")
        .size()
        .reindex(g.index)
        .fillna(0)
    )
    return g, res


def signature_test(gl, sig):
    from scipy.stats import mannwhitneyu

    out = {}
    for name, alt in (("recall_up", "greater"), ("recall_down", "less")):
        genes = [g for g in sig.loc[sig["set"] == name, "gene"] if g in gl.index]
        rest = gl.drop(index=genes)
        p = (
            mannwhitneyu(gl.loc[genes, "mean_stat"], rest["mean_stat"], alternative=alt).pvalue
            if genes
            else 1.0
        )
        out[name] = {
            "n_genes": len(genes),
            "median_stat_signature": float(gl.loc[genes, "mean_stat"].median())
            if genes
            else float("nan"),
            "median_stat_background": float(rest["mean_stat"].median()),
            "p_one_sided": float(p),
        }
    return out


def figure(gls: dict, sig, out_dir: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(gls), figsize=(5 * len(gls), 4), squeeze=False)
    for ax, (_name, (label, gl)) in zip(axes[0], gls.items(), strict=True):
        up = gl.loc[gl.index.isin(sig.loc[sig["set"] == "recall_up", "gene"]), "mean_stat"]
        dn = gl.loc[gl.index.isin(sig.loc[sig["set"] == "recall_down", "gene"]), "mean_stat"]
        bg = gl.loc[~gl.index.isin(sig["gene"]), "mean_stat"]
        parts = ax.violinplot([bg.clip(-8, 8), up.clip(-8, 8), dn.clip(-8, 8)], showmedians=True)
        for body, col in zip(parts["bodies"], ["#9ca3af", "#d97706", "#2563eb"], strict=True):
            body.set_facecolor(col)
            body.set_alpha(0.6)
        ax.axhline(0, c="k", lw=0.5)
        ax.set_xticks(
            [1, 2, 3],
            [
                f"all genes\n(n={len(bg):,})",
                f"recall_up\n(n={len(up)})",
                f"recall_down\n(n={len(dn)})",
            ],
            fontsize=8,
        )
        ax.set_ylabel("gene accessibility change (mean Wald stat)")
        ax.set_title(label, fontsize=10)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "atac_signature_accessibility.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--atac-dir", type=Path, default=Path("data/raw/GSE172258"))
    p.add_argument("--refgene", type=Path, default=Path("references/ucsc_mm10_refGene.txt.gz"))
    p.add_argument("--signature", type=Path, default=SIGNATURE)
    p.add_argument("--out-dir", type=Path, default=Path("results/regulatory/atac"))
    p.add_argument("--report", type=Path, default=Path("reports/phase9/09_chromatin_report.md"))
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    import pandas as pd

    missing = [p for p in (args.refgene, args.signature) if not p.is_file()]
    files = {c.file for c in CONTRASTS}
    missing += [args.atac_dir / f for f in files if not (args.atac_dir / f).is_file()]
    if missing:
        status = "NOT PERFORMED — missing input(s): " + ", ".join(map(str, missing))
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(f"# Phase 9 — chromatin\n\n**Epigenomic integration: {status}**\n")
        print(f"[atac] {status}", file=sys.stderr)
        return 2

    sig = pd.read_csv(args.signature, sep="\t")
    tss = load_tss(args.refgene)
    tables = {f: pd.read_csv(args.atac_dir / f, sep="\t", index_col=0) for f in files}
    args.out_dir.mkdir(parents=True, exist_ok=True)

    rows, gls = [], {}
    for c in CONTRASTS:
        res = differential(tables[c.file], c)
        gl, peaks = gene_level(res, tss)
        peaks.to_csv(args.out_dir / f"da_peaks_{c.name}.csv.gz")
        gl.to_csv(args.out_dir / f"gene_accessibility_{c.name}.csv")
        gls[c.name] = (c.label, gl)
        t = signature_test(gl, sig)
        n_da = int((res["padj"] < 0.05).sum())
        rows.append(
            {
                "comparison": c.label,
                "recall": c.recall,
                "n": f"{res.attrs['n_test']} v {res.attrs['n_ref']}",
                "peaks_tested": len(res),
                "DA_peaks_FDR5": n_da,
                "up_genes": t["recall_up"]["n_genes"],
                "p_up_gain": t["recall_up"]["p_one_sided"],
                "down_genes": t["recall_down"]["n_genes"],
                "p_down_loss": t["recall_down"]["p_one_sided"],
            }
        )
        print(
            f"[atac] {c.name}: {n_da} DA peaks; p_up={t['recall_up']['p_one_sided']:.3g} "
            f"p_down={t['recall_down']['p_one_sided']:.3g}"
        )
    summary = pd.DataFrame(rows).set_index("comparison")
    summary.to_csv(args.out_dir / "signature_accessibility_summary.csv")
    recall = summary[summary["recall"]]
    n_ok = int((recall["p_up_gain"] < P_SUPPORT).sum())
    verdict = {
        2: "REGULATORY SUPPORT",
        1: "PARTIAL REGULATORY SUPPORT",
        0: "NO REGULATORY SUPPORT",
    }[min(n_ok, 2)]

    # per-gene view for the signature in the two recall comparisons
    per_gene = sig[["set", "gene"]].set_index("gene")
    for c in CONTRASTS:
        per_gene[f"stat_{c.name}"] = gls[c.name][1]["mean_stat"].reindex(per_gene.index)
    per_gene.to_csv(args.out_dir / "signature_gene_accessibility.csv")
    consistent = per_gene[
        (
            (per_gene["set"] == "recall_up")
            & (per_gene["stat_rag_recall"] > 2)
            & (per_gene["stat_wt_recall"] > 2)
        )
        | (
            (per_gene["set"] == "recall_down")
            & (per_gene["stat_rag_recall"] < -2)
            & (per_gene["stat_wt_recall"] < -2)
        )
    ]

    fig = figure(gls, sig, args.out_dir / "figures")
    (args.report.parent / "figures").mkdir(parents=True, exist_ok=True)
    shutil.copy2(fig, args.report.parent / "figures" / fig.name)

    def md(df):
        cols = [df.index.name or "", *map(str, df.columns)]
        out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
        for i, r in df.iterrows():
            out.append(
                "| "
                + " | ".join([str(i)] + [f"{v:.3g}" if isinstance(v, float) else str(v) for v in r])
                + " |"
            )
        return "\n".join(out)

    report = "\n".join(
        [
            "# Phase 9 — epigenomic integration (GSE172258 ATAC-seq)",
            "",
            "- Status: **performed** (compatible chromatin data identified in Phase 1).",
            f"- **Verdict (pre-registered, D-014): {verdict}** — `recall_up` genes gain accessibility "
            f"(one-sided p < {P_SUPPORT}) in {n_ok} of 2 recall comparisons.",
            f"- Peaks linked to the nearest RefSeq TSS (mm10) within {MAX_DIST // 1000} kb; gene statistic = mean "
            "DESeq2 Wald statistic of its peaks; pydeseq2 with default settings.",
            "- Motif enrichment: **not performed** — no genome sequence / motif database in this project (D-014).",
            "- Cross-model caveat: Alternaria allergen memory (Rag1⁻/⁻ and WT) vs *N. brasiliensis* recall in "
            "the discovery data — agreement indicates a shared program, not identical biology.",
            "",
            "## Signature vs accessibility",
            md(summary),
            "",
            f"## Signature genes with consistent accessibility in both recall comparisons (|stat| > 2): {len(consistent)}",
            ", ".join(f"{g} ({r['set'].replace('recall_', '')})" for g, r in consistent.iterrows())
            or "none",
            "",
            "## Figures",
            f"![{fig.stem}](figures/{fig.name})",
            "",
        ]
    )
    args.report.write_text(report, encoding="utf-8")
    print(f"[atac] verdict: {verdict} · report → {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
