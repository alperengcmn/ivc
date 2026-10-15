"""Script 10 — Phase 10: cross-dataset validation of the frozen signature (pre-registered, D-014).

Rules are in docs/validation_plan.md and were committed before this module was run.

Core memory/recall datasets get a verdict (Replicated / Partially supported / Not supported /
Not testable). The acute-influenza dataset annotates which signature genes are shared with
acute activation. Mechanistic datasets (BATF, iron, Blimp-1/stimulation) are enrichment context.
"""

from __future__ import annotations

import argparse
import gzip
import io
import shutil
import sys
from collections.abc import Sequence
from pathlib import Path

RAW = Path("data/raw")
SIGNATURE = Path("signatures/GSE295590_recall_ILC2_signature_v1.tsv")
MIN_GENES = 20
P = 0.05
MIN_DETECT = 0.05  # gene measurable if detected in >= 5 % of units (count data)
MIN_WELL_UMI = 1000  # plate wells below this are failed libraries (applied before any score)


# ----------------------------------------------------------------------------- metrics


def signature_sets(path: Path = SIGNATURE):
    import pandas as pd

    sig = pd.read_csv(path, sep="\t")
    return list(sig.loc[sig["set"] == "recall_up", "gene"]), list(
        sig.loc[sig["set"] == "recall_down", "gene"]
    )


def score_units(logexpr, up, down):
    """logexpr: DataFrame genes x units (log scale). z per gene across units; up − down."""
    z = logexpr.sub(logexpr.mean(axis=1), axis=0).div(logexpr.std(axis=1).replace(0, 1), axis=0)
    u = [g for g in up if g in z.index]
    d = [g for g in down if g in z.index]
    return z.loc[u].mean() - z.loc[d].mean(), len(u) + len(d)


def gene_agreement(lfc, up, down):
    from scipy.stats import binomtest

    u = lfc.reindex(up).dropna()
    d = lfc.reindex(down).dropna()
    agree = int((u > 0).sum() + (d < 0).sum())
    n = len(u) + len(d)
    p = binomtest(agree, n, 0.5, alternative="greater").pvalue if n else 1.0
    p_opp = binomtest(agree, n, 0.5, alternative="less").pvalue if n else 1.0
    return {
        "n_genes": n,
        "agree": agree,
        "fraction": agree / n if n else float("nan"),
        "p": p,
        "p_opposite": p_opp,
        "up_agree": f"{int((u > 0).sum())}/{len(u)}",
        "down_agree": f"{int((d < 0).sum())}/{len(d)}",
    }


def score_test(scores, test_mask, ref_mask):
    from scipy.stats import mannwhitneyu

    a, b = scores[test_mask], scores[ref_mask]
    if len(a) < 3 or len(b) < 3:
        return {
            "n_test": len(a),
            "n_ref": len(b),
            "delta": float(a.mean() - b.mean()),
            "p": None,
            "p_opposite": None,
        }
    return {
        "n_test": len(a),
        "n_ref": len(b),
        "delta": float(a.mean() - b.mean()),
        "p": float(mannwhitneyu(a, b, alternative="greater").pvalue),
        "p_opposite": float(mannwhitneyu(a, b, alternative="less").pvalue),
    }


def verdict(score: dict, genes: dict) -> str:
    """Pre-registered rule (docs/validation_plan.md)."""
    if genes["n_genes"] < MIN_GENES:
        return "Not testable with available data"
    opposite = genes["p_opposite"] < P or (
        score["p_opposite"] is not None and score["p_opposite"] < P
    )
    if opposite:
        return "Not supported"
    s_ok = score["p"] is not None and score["p"] < P
    g_ok = genes["p"] < P
    if score["p"] is None:
        return "Partially supported" if g_ok else "Not supported"
    if s_ok and g_ok:
        return "Replicated"
    if s_ok or g_ok:
        return "Partially supported"
    return "Not supported"


# ----------------------------------------------------------------------------- datasets


def gse172258_rna(up, down):
    import numpy as np
    import pandas as pd

    path = RAW / "GSE172258" / "GSE172258_expression_matrix_mgi.txt.gz"
    m = pd.read_csv(path, sep="\t", index_col=0)
    m.index = [str(i).split(" ")[-1] for i in m.index]
    m = m.groupby(level=0).sum()
    grp = pd.Series(m.columns, index=m.columns).str.extract(r"[._](\d[AB])_(\d+)cells")
    grp.columns = ["group", "cells"]
    keep = grp["group"].notna() & (m.sum() >= MIN_WELL_UMI)
    m, grp = m.loc[:, keep], grp[keep]
    det = (m > 0).mean(axis=1)
    m = m[det >= MIN_DETECT]
    logx = np.log1p(m / m.sum() * 1e4)
    scores, _ = score_units(logx, up, down)
    out = []
    for test, ref, label in (
        ("1B", "2B", "after recall: Alt-sensitised (1B) vs saline (2B)"),
        ("1A", "2A", "before recall: Alt-sensitised (1A) vs saline (2A)"),
    ):
        t, r = (grp["group"] == test).to_numpy(), (grp["group"] == ref).to_numpy()
        lfc = logx.loc[:, t].mean(axis=1) - logx.loc[:, r].mean(axis=1)
        s, g = score_test(scores, t, r), gene_agreement(lfc, up, down)
        strata = {}
        for cells in ("20", "40"):
            cm = (grp["cells"] == cells).to_numpy()
            strata[cells] = score_test(scores, t & cm, r & cm)["delta"]
        out.append(
            {
                "dataset": "GSE172258 RNA",
                "comparison": label,
                "unit": "20/40-cell well",
                "score": s,
                "genes": g,
                "verdict": verdict(s, g),
                "note": f"Δscore 20-cell wells {strata['20']:.2f}, "
                f"40-cell wells {strata['40']:.2f}; {int(keep.sum())} wells ≥ {MIN_WELL_UMI} UMI",
                "scores": scores,
                "groups": grp["group"],
            }
        )
    return out


def _read_series_matrix(path: Path):
    import pandas as pd

    meta, lines, in_table = {}, [], False
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                lines.append(line)
            elif line.startswith("!Sample_"):
                key = line.split("\t")[0]
                vals = [v.strip().strip('"') for v in line.rstrip("\n").split("\t")[1:]]
                meta.setdefault(key, []).append(vals)
    table = pd.read_csv(io.StringIO("".join(lines)), sep="\t", index_col=0)
    return meta, table


def gse81700(up, down, platform: Path):
    import pandas as pd

    if not platform.is_file():
        return [
            {
                "dataset": "GSE81700",
                "comparison": "lung ILC2 memory-like vs naive",
                "unit": "mouse",
                "score": {
                    "n_test": 0,
                    "n_ref": 0,
                    "delta": float("nan"),
                    "p": None,
                    "p_opposite": None,
                },
                "genes": {
                    "n_genes": 0,
                    "agree": 0,
                    "fraction": float("nan"),
                    "p": 1.0,
                    "p_opposite": 1.0,
                    "up_agree": "-",
                    "down_agree": "-",
                },
                "verdict": "Not testable with available data",
                "note": "platform annotation (references/GPL16570_symbols.tsv) not available",
            }
        ]
    meta, x = _read_series_matrix(RAW / "GSE81700" / "GSE81700_series_matrix.txt.gz")
    sym = pd.read_csv(platform, sep="\t", dtype=str).set_index("ID")["gene_symbol"]
    x.index = x.index.astype(str)
    x = x.join(sym.rename("gene"), how="inner").dropna(subset=["gene"])
    x = x.groupby("gene").mean()
    titles = meta["!Sample_title"][0]
    tissue = next(v for v in meta["!Sample_characteristics_ch1"] if v[0].startswith("tissue"))
    tp = next(v for v in meta["!Sample_characteristics_ch1"] if v[0].startswith("time point"))
    lung = [i for i, t in enumerate(tissue) if t == "tissue: lung"]
    x = x.iloc[:, lung]
    tp = pd.Series([tp[i].split(": ")[1] for i in lung], index=x.columns)
    x = x[x.mean(axis=1) > x.mean(axis=1).quantile(0.25)]  # drop the least-expressed quartile
    scores, _ = score_units(x, up, down)
    out = []
    for test in ("4 months", "week 2", "day 4"):
        t, r = (tp == test).to_numpy(), (tp == "naive").to_numpy()
        lfc = x.loc[:, t].mean(axis=1) - x.loc[:, r].mean(axis=1)
        s, g = score_test(scores, t, r), gene_agreement(lfc, up, down)
        out.append(
            {
                "dataset": "GSE81700",
                "comparison": f"lung ILC2 {test} after IL-33 vs naive",
                "unit": "mouse",
                "score": s,
                "genes": g,
                "verdict": verdict(s, g),
                "note": f"{len(titles)} arrays; lung only",
                "scores": scores,
                "groups": tp,
            }
        )
    return out


def _read_mtx(prefix: Path):
    import pandas as pd
    import scipy.io
    import scipy.sparse as sp

    X = sp.csr_matrix(scipy.io.mmread(str(prefix) + "matrix.mtx.gz")).T.tocsr()
    genes = pd.read_csv(str(prefix) + "genes.tsv.gz", sep="\t", header=None)
    bc = pd.read_csv(str(prefix) + "barcodes.tsv.gz", sep="\t", header=None)[0]
    return X, genes[1].astype(str).to_numpy(), bc.to_numpy()


def gse149857_acute(up, down, markers):
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy.stats import spearmanr
    from statsmodels.stats.multitest import multipletests

    parts = {}
    for day, prefix in (
        ("day5", "GSM4514296_scRNAseq_Day5_Lung_ILC_"),
        ("day10", "GSM4514297_scRNAseq_Day10_Lung_ILC_"),
    ):
        X, genes, bc = _read_mtx(RAW / "GSE149857" / "RAW" / prefix)
        a = ad.AnnData(
            X.astype(np.float32), obs=pd.DataFrame(index=bc), var=pd.DataFrame(index=genes)
        )
        a.var_names_make_unique()
        a.obs["day"] = day
        parts[day] = a
    sig_genes = set(up) | set(down)
    rows = {}
    n_ilc2 = {}
    for day, a in parts.items():
        sc.pp.filter_cells(a, min_genes=200)
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        for name, genes in markers["lineage"].items():
            g = [x for x in genes if x in a.var_names]
            if len(g) >= 2:
                sc.tl.score_genes(a, g, score_name=f"lin_{name}", random_state=0)
        lin = [c for c in a.obs if c.startswith("lin_")]
        best = a.obs[lin].idxmax(axis=1)
        ilc = a[(best == "lin_ILC2") & (a.obs["lin_ILC2"] > 0.1)].copy()
        n_ilc2[day] = ilc.n_obs
        axis_genes = [
            g
            for k in ("type2_effector", "activation", "proliferation")
            for g in markers["state"][k]
            if g in ilc.var_names and g not in sig_genes
        ]
        sc.tl.score_genes(ilc, axis_genes, score_name="activation_axis", random_state=0)
        ax = ilc.obs["activation_axis"].to_numpy()
        res = {}
        for g in up + down:
            if g in ilc.var_names:
                v = ilc[:, g].X.toarray().ravel()
                if (v > 0).mean() >= MIN_DETECT:
                    rho, p = spearmanr(v, ax)
                    res[g] = (rho, p)
        df = pd.DataFrame(res, index=["rho", "p"]).T
        df["fdr"] = multipletests(df["p"], method="fdr_bh")[1]
        rows[day] = df
    t = rows["day5"].add_suffix("_d5").join(rows["day10"].add_suffix("_d10"), how="outer")
    exp = pd.Series({g: 1 for g in up} | {g: -1 for g in down})
    t["expected_sign"] = exp.reindex(t.index)
    shared = True
    for d in ("d5", "d10"):
        shared = (
            shared
            & (np.sign(t[f"rho_{d}"]) == t["expected_sign"])
            & (t[f"rho_{d}"].abs() >= 0.1)
            & (t[f"fdr_{d}"] < 0.05)
        )
    t["acute_shared"] = shared.fillna(False)
    opposite = True
    for d in ("d5", "d10"):
        opposite = (
            opposite & (np.sign(t[f"rho_{d}"]) == -t["expected_sign"]) & (t[f"fdr_{d}"] < 0.05)
        )
    t["acute_opposite"] = opposite.fillna(False)
    t["set"] = np.where(t["expected_sign"] > 0, "recall_up", "recall_down")
    return t, n_ilc2


def _deseq(counts, meta, design, contrast):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    counts = counts.loc[:, counts.sum() >= 10 * counts.shape[0]]
    dds = DeseqDataSet(
        counts=counts.round().astype(int),
        metadata=meta,
        design=design,
        inference=DefaultInference(n_cpus=4),
        quiet=True,
    )
    dds.deseq2()
    st = DeseqStats(dds, contrast=contrast, inference=DefaultInference(n_cpus=4), quiet=True)
    st.summary()
    return st.results_df


def fisher_context(res, genes, lfc_sign: int | None = None):
    from scipy.stats import fisher_exact

    res = res.dropna(subset=["padj"])
    hit = res["padj"] < 0.05
    if lfc_sign is not None:
        hit &= (res["log2FoldChange"] * lfc_sign) > 0
    in_set = res.index.isin(genes)
    a, b = int((hit & in_set).sum()), int((~hit & in_set).sum())
    c, d = int((hit & ~in_set).sum()), int((~hit & ~in_set).sum())
    odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return {"in_set": f"{a}/{a + b}", "background": f"{c}/{c + d}", "odds_ratio": odds, "p": p}


def context(up, down):
    import pandas as pd

    rows = []
    # GSE149854 — authors' DESeq2 table (KO vs WT)
    # gene symbols are in the first (unnamed) column; "Gene_name" holds the description
    x = pd.read_excel(RAW / "GSE149854" / "GSE149854_ILC2_Gene_Result_new.xlsx")
    x = x.set_index(x.columns[0])
    x = x[~x.index.duplicated()]
    sign_ok = ((x["log2FoldChange"] > 0) == (x["mean_fpkm_KO"] > x["mean_fpkm_WT"])).mean()
    f = fisher_context(x, up, lfc_sign=-1 if sign_ok > 0.9 else 1)
    rows.append(
        {
            "dataset": "GSE149854",
            "question": "recall_up genes reduced in Batf KO (BATF-dependent)",
            **f,
            "note": f"authors' DESeq2 table; log2FC = KO/WT (sign check {sign_ok:.0%})",
        }
    )
    # GSE262210 — iron deprivation (DEF vs vehicle)
    c = pd.read_csv(
        RAW / "GSE262210" / "GSE262210_ILC2iron_gene_counts.txt.gz", sep="\t", index_col=0
    )
    sm, _ = _read_series_matrix(RAW / "GSE262210" / "GSE262210_series_matrix.txt.gz")
    titles = dict(
        zip(sm["!Sample_description"][0], sm["!Sample_title"][0], strict=True)
    )  # column -> title
    cond = {"vehicle": "vehicle", "DEF": "DEF", "DEF+FeCl3": "DEF_Fe"}
    meta = pd.DataFrame(
        {"cond": [cond[titles[col].rsplit("_rep", 1)[0]] for col in c.columns]}, index=c.columns
    )
    sel = meta["cond"].isin(["vehicle", "DEF"])
    r = _deseq(c.loc[:, sel].T, meta[sel], "~cond", ["cond", "DEF", "vehicle"])
    for name, genes in (("recall_up", up), ("recall_down", down)):
        f = fisher_context(r, genes)
        rows.append(
            {
                "dataset": "GSE262210",
                "question": f"{name} genes changed by iron deprivation (DEF vs vehicle)",
                **f,
                "note": "DESeq2, 3 v 3",
            }
        )
    # GSE291231 — stimulation (paired by mouse) and Blimp-1 KO
    s = pd.read_csv(RAW / "GSE291231" / "GSE291231_Stimulation_Raw_reads.csv.gz", index_col=0).drop(
        columns="Length"
    )
    s = s.groupby(level=0).sum()
    meta = pd.DataFrame(
        {
            "mouse": [x.split("_")[0] for x in s.columns],
            "cond": [x.split("_")[1] for x in s.columns],
        },
        index=s.columns,
    )
    r = _deseq(s.T, meta, "~mouse + cond", ["cond", "stim", "unstim"])
    g = gene_agreement(r["log2FoldChange"], up, down)
    rows.append(
        {
            "dataset": "GSE291231",
            "question": "signature direction under in-vitro stimulation (activation-shared?)",
            "in_set": f"{g['agree']}/{g['n_genes']} agree",
            "background": "-",
            "odds_ratio": float("nan"),
            "p": g["p"],
            "note": "DESeq2 paired by mouse, 4 v 4",
        }
    )
    k = pd.read_csv(RAW / "GSE291231" / "GSE291231_KO_Raw_reads.csv.gz", index_col=0).drop(
        columns="Length"
    )
    k = k.groupby(level=0).sum()
    meta = pd.DataFrame({"geno": [x.split("_")[0] for x in k.columns]}, index=k.columns)
    r = _deseq(k.T, meta, "~geno", ["geno", "KO", "WT"])
    f = fisher_context(r, up + down)
    rows.append(
        {
            "dataset": "GSE291231",
            "question": "signature genes changed in Blimp-1 KO",
            **f,
            "note": "DESeq2, 3 v 3",
        }
    )
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- report


def figure(results, out_dir: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels = [r for r in results if "scores" in r]
    seen, uniq = set(), []
    for r in panels:
        if r["dataset"] not in seen:
            uniq.append(r)
            seen.add(r["dataset"])
    fig, axes = plt.subplots(1, len(uniq), figsize=(6 * len(uniq), 4), squeeze=False)
    for ax, r in zip(axes[0], uniq, strict=True):
        order = {
            "GSE172258 RNA": ["2A", "1A", "2B", "1B"],
            "GSE81700": ["naive", "day 4", "week 2", "4 months"],
        }
        groups = order.get(r["dataset"], sorted(r["groups"].unique()))
        data = [r["scores"][(r["groups"] == g).to_numpy()].to_numpy() for g in groups]
        ax.boxplot(data, showfliers=False)
        for i, d in enumerate(data, start=1):
            ax.scatter([i] * len(d), d, s=6, alpha=0.4, c="#d97706")
        ax.set_xticks(range(1, len(groups) + 1), groups)
        ax.axhline(0, c="k", lw=0.5)
        ax.set_ylabel("frozen signature score (up − down)")
        ax.set_title(r["dataset"])
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "validation_scores.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _fmt(v):
    if v is None:
        return "n/a"
    if isinstance(v, float):
        return f"{v:.3g}"
    return str(v)


def render(results, acute, n_ilc2, ctx, fig, overall) -> str:
    lines = [
        "# Phase 10 — cross-dataset validation of the frozen signature",
        "",
        "Pre-registered rules: [`docs/validation_plan.md`](../../docs/validation_plan.md) (D-014). Signature: "
        f"`{SIGNATURE.name}`.",
        "",
        f"## Overall: **{overall}**",
        "",
        "## Core memory / recall datasets",
        "| dataset | comparison | units (test v ref) | Δ score | score p | gene agreement | up / down | gene p | verdict |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        s, g = r["score"], r["genes"]
        lines.append(
            f"| {r['dataset']} | {r['comparison']} | {s['n_test']} v {s['n_ref']} | {_fmt(s['delta'])} | "
            f"{_fmt(s['p'])} | {g['agree']}/{g['n_genes']} ({_fmt(g['fraction'])}) | {g['up_agree']} / "
            f"{g['down_agree']} | {_fmt(g['p'])} | **{r['verdict']}** |"
        )
    lines += ["", "Notes: " + "; ".join(f"{r['dataset']} — {r['note']}" for r in results), ""]
    n_sh = int(acute["acute_shared"].sum())
    by = acute.groupby("set")["acute_shared"].agg(["sum", "size"])
    lines += [
        "## Acute-activation reference — GSE149857 (influenza, lung ILC2s)",
        f"ILC2s analysed: day 5 {n_ilc2['day5']:,}, day 10 {n_ilc2['day10']:,}. Signature genes correlated with an "
        f"acute activation axis (without signature genes) in the signature's direction on both days: **{n_sh}** "
        f"(recall_up {int(by.loc['recall_up', 'sum'])}/{int(by.loc['recall_up', 'size'])}, recall_down "
        f"{int(by.loc['recall_down', 'sum'])}/{int(by.loc['recall_down', 'size'])} measurable).",
        "",
        "- **Acute-shared:** " + (", ".join(acute.index[acute["acute_shared"]]) or "none"),
        "- **Recall-specific candidates (measurable, not acute-shared):** "
        + ", ".join(
            acute.index[~acute["acute_shared"] & acute["rho_d5"].notna() & acute["rho_d10"].notna()]
        ),
        "",
        "- **Opposite to acute activation** (signature direction opposite to the acute correlation, FDR < 0.05 on "
        "both days; these separate recall from acute activation most clearly): "
        + (", ".join(acute.index[acute["acute_opposite"]]) or "none"),
        "",
        "## Mechanistic and context datasets",
        "| dataset | question | in signature | background | odds ratio | p | note |",
        "|---|---|---|---|---|---|---|",
        *[
            f"| {r.dataset} | {r.question} | {r.in_set} | {r.background} | {_fmt(r.odds_ratio)} | {_fmt(r.p)} | {r.note} |"
            for r in ctx.itertuples()
        ],
        "",
        "## Figures",
        f"![{fig.stem}](figures/{fig.name})",
        "",
    ]
    return "\n".join(lines)


def overall_category(results) -> str:
    core = [r["verdict"] for r in results]
    if "Replicated" in core:
        return "Replicated in ≥ 1 independent memory/recall comparison (see table for the rest)"
    if "Partially supported" in core:
        return "Partially supported"
    if all(v == "Not testable with available data" for v in core):
        return "Not testable with available data"
    return "Not supported"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--platform", type=Path, default=Path("references/GPL16570_symbols.tsv"))
    p.add_argument("--out-dir", type=Path, default=Path("results/validation"))
    p.add_argument("--report", type=Path, default=Path("reports/phase10/10_validation_report.md"))
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    import pandas as pd

    from .annotate import load_markers

    if not SIGNATURE.is_file():
        print(f"Error: {SIGNATURE} missing — run Phase 7 first.", file=sys.stderr)
        return 2
    up, down = signature_sets()
    results = gse172258_rna(up, down) + gse81700(up, down, args.platform)
    for r in results:
        print(f"[validate] {r['dataset']} · {r['comparison']}: {r['verdict']}")
    acute, n_ilc2 = gse149857_acute(up, down, load_markers())
    ctx = context(up, down)
    overall = overall_category(results)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {
                "dataset": r["dataset"],
                "comparison": r["comparison"],
                **{f"score_{k}": v for k, v in r["score"].items()},
                **{f"genes_{k}": v for k, v in r["genes"].items()},
                "verdict": r["verdict"],
                "note": r["note"],
            }
            for r in results
        ]
    ).to_csv(args.out_dir / "cross_dataset_comparison.csv", index=False)
    acute.to_csv(args.out_dir / "acute_shared_annotation.csv")
    ctx.to_csv(args.out_dir / "context_datasets.csv", index=False)
    fig = figure(results, args.out_dir / "figures")
    (args.report.parent / "figures").mkdir(parents=True, exist_ok=True)
    shutil.copy2(fig, args.report.parent / "figures" / fig.name)
    args.report.write_text(render(results, acute, n_ilc2, ctx, fig, overall), encoding="utf-8")
    print(f"[validate] overall: {overall} · report → {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
