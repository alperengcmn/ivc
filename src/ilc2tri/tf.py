"""Script 08 — Phase 8: transcription-factor and pathway activity.

TF activity is inferred with a univariate linear model (decoupler ULM) on the CollecTRI
mouse TF→target network, applied to each Phase 6 contrast statistic (Wilcoxon score), and
PROGENy pathway footprints are scored the same way. A TF is *concordant* when FDR < 0.05
in both contrasts with the same sign; it is *target-supported* when at least
MIN_TARGETS of its targets are Phase 6 concordant candidates and at least MIN_AGREEMENT of
them move in the direction implied by the TF's activity and the edge sign.

Activity inferred from expression is not evidence of TF binding or function.
Lightweight by design: ULM replaces SCENIC (pipeline flagged it as heavy); D-013.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from collections.abc import Sequence
from pathlib import Path

from .gse295590 import SERIES

FDR = 0.05
MIN_TARGETS = 5
MIN_AGREEMENT = 0.6
# TFs named by the independent datasets' own studies (GSE172258 memory ILC2; GSE149854 BATF).
LITERATURE_TFS = {
    "Bach2": "GSE172258",
    "Jund": "GSE172258",
    "Fosl2": "GSE172258",
    "Fhl2": "GSE172258",
    "Fosb": "GSE172258",
    "Stat6": "GSE172258",
    "Srebf2": "GSE172258",
    "Batf": "GSE149854",
    "Gata3": "ILC2 identity",
    "Rora": "ILC2 identity",
}


def contrast_scores(de_dir: Path):
    import pandas as pd

    cols = {}
    for c in ("A", "B"):
        t = pd.read_csv(de_dir / f"contrast_{c}.csv", index_col=0)
        cols[c] = t["score"]
    return pd.DataFrame(cols).T  # rows = contrasts, columns = genes (NaN where untested)


def activity(mat, net, tmin: int = 5):
    import decoupler as dc
    import numpy as np
    import pandas as pd
    from statsmodels.stats.multitest import multipletests

    out = {}
    for c in mat.index:
        row = mat.loc[[c]].dropna(axis=1)
        es, pv = dc.mt.ulm(data=row, net=net, tmin=tmin)
        fdr = multipletests(pv.loc[c].to_numpy(), method="fdr_bh")[1]
        out[c] = pd.DataFrame({"score": es.loc[c], "pval": pv.loc[c], "fdr": fdr})
    a, b = out["A"], out["B"]
    tfs = a.index.intersection(b.index)
    t = a.loc[tfs].add_suffix("_A").join(b.loc[tfs].add_suffix("_B"))
    t["concordant"] = (
        (t["fdr_A"] < FDR) & (t["fdr_B"] < FDR) & (np.sign(t["score_A"]) == np.sign(t["score_B"]))
    )
    t["direction"] = np.where(t["score_A"] > 0, "up_in_secondary", "down_in_secondary")
    return t


def target_support(tf_table, net, de):
    import numpy as np
    import pandas as pd

    cand = de[de["candidate"]]
    lfc = cand["pb_log2fc_A"]
    rows = []
    for tf, r in tf_table.iterrows():
        edges = net[(net["source"] == tf) & net["target"].isin(cand.index)]
        if edges.empty:
            rows.append((tf, 0, np.nan, ""))
            continue
        expected = np.sign(r["score_A"]) * np.sign(edges["weight"].to_numpy())
        observed = np.sign(lfc.loc[edges["target"]].to_numpy())
        agree = float((expected == observed).mean())
        ok = edges["target"][expected == observed].tolist()
        rows.append((tf, len(edges), agree, ", ".join(ok[:12])))
    s = pd.DataFrame(
        rows, columns=["tf", "n_candidate_targets", "direction_agreement", "supporting_targets"]
    )
    return s.set_index("tf")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--de-dir", type=Path, default=Path("results/differential_expression"))
    p.add_argument(
        "--annotated", type=Path, default=Path("results/annotation") / f"{SERIES}_annotated.h5ad"
    )
    p.add_argument("--collectri", type=Path, default=Path("references/collectri_mouse.tsv"))
    p.add_argument("--progeny", type=Path, default=Path("references/progeny_mouse.tsv"))
    p.add_argument("--out-dir", type=Path, default=Path("results/regulatory"))
    p.add_argument("--report", type=Path, default=Path("reports/phase8/08_tf_activity_report.md"))
    return p


def per_cell_figure(ilc, net, tfs: list[str], out_dir: Path) -> Path:
    import decoupler as dc
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    sub_net = net[net["source"].isin(tfs)]
    genes = ilc.var_names.intersection(sub_net["target"].unique())
    X = ilc[:, genes].layers["lognorm"]
    mat = pd.DataFrame(
        X.toarray(), index=ilc.obs_names, columns=genes
    )  # small: targets of <= 8 TFs only
    scores, _ = dc.mt.ulm(data=mat, net=sub_net, tmin=5)
    libs = [("GSM8953514", "S1 prim."), ("GSM8953515", "S2 prim."), ("GSM8953516", "S3 sec.")]
    tfs = [t for t in tfs if t in scores.columns][:8]
    fig, axes = plt.subplots(2, 4, figsize=(15, 6))
    for ax, tf in zip(axes.ravel(), tfs, strict=False):
        data = [
            scores.loc[(ilc.obs["sample"].astype(str) == s).to_numpy(), tf].to_numpy()
            for s, _ in libs
        ]
        parts = ax.violinplot(data, showmedians=True)
        for body, col in zip(parts["bodies"], ["#2563eb", "#60a5fa", "#d97706"], strict=True):
            body.set_facecolor(col)
            body.set_alpha(0.6)
        ax.set_xticks([1, 2, 3], [n for _, n in libs], fontsize=8)
        ax.set_title(tf)
    for ax in axes.ravel()[len(tfs) :]:
        ax.axis("off")
    fig.suptitle(
        "Per-cell TF activity (ULM, CollecTRI) — top concordant TFs; S1 is shallower (full-depth data)"
    )
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "tf_activity_per_cell.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _md(df, fmt="{:.3g}") -> str:
    cols = [df.index.name or "", *map(str, df.columns)]
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for i, row in df.iterrows():
        cells = [
            ("✓" if v else "")
            if isinstance(v, bool)
            else (fmt.format(v) if isinstance(v, float) else str(v))
            for v in row
        ]
        out.append("| " + " | ".join([str(i), *cells]) + " |")
    return "\n".join(out)


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    import pandas as pd
    import scanpy as sc

    for f in (args.collectri, args.progeny):
        if not f.is_file():
            print(f"Error: {f} missing — run `make fetch-resources`.", file=sys.stderr)
            return 2
    net = pd.read_csv(args.collectri, sep="\t")[["source", "target", "weight"]]
    prog = pd.read_csv(args.progeny, sep="\t")[["source", "target", "weight"]]
    de = pd.read_csv(args.de_dir / "differential_expression_results.csv", index_col=0)
    mat = contrast_scores(args.de_dir)

    tf = activity(mat, net)
    sup = target_support(tf, net, de)
    tf = tf.join(sup)
    tf["target_supported"] = (tf["n_candidate_targets"] >= MIN_TARGETS) & (
        tf["direction_agreement"] >= MIN_AGREEMENT
    )
    tf["tf_gene_log2fc_A"] = de["pb_log2fc_A"].reindex(tf.index)
    tf["tf_gene_candidate"] = de["candidate"].reindex(tf.index).fillna(False).astype(bool)
    tf["min_abs_score"] = tf[["score_A", "score_B"]].abs().min(axis=1)
    prioritized = tf[tf["concordant"] & tf["target_supported"]].sort_values(
        "min_abs_score", ascending=False
    )
    pw = activity(mat, prog, tmin=5)

    lit = tf.reindex(list(LITERATURE_TFS))[
        ["score_A", "fdr_A", "score_B", "fdr_B", "concordant", "target_supported"]
    ]
    lit.insert(0, "source", [LITERATURE_TFS[t] for t in lit.index])
    lit.index.name = "tf"

    args.out_dir.mkdir(parents=True, exist_ok=True)
    tf.to_csv(args.out_dir / "tf_activity_matrix.csv")
    prioritized.to_csv(args.out_dir / "prioritized_tf_candidates.tsv", sep="\t")
    pw.to_csv(args.out_dir / "progeny_pathway_activity.csv")
    edges = net[
        net["source"].isin(prioritized.index) & net["target"].isin(de.index[de["candidate"]])
    ]
    edges.to_csv(args.out_dir / "regulatory_network_table.tsv", sep="\t", index=False)

    adata = sc.read_h5ad(args.annotated)
    for layer in ("counts", "depth_matched", "lognorm_dm"):
        adata.layers.pop(layer, None)
    adata.X = adata.layers["lognorm"]
    ilc = adata[adata.obs["lineage"].astype(str) == "ILC2"].copy()
    del adata
    fig = per_cell_figure(ilc, net, list(prioritized.index[:8]), args.out_dir / "figures")
    (args.report.parent / "figures").mkdir(parents=True, exist_ok=True)
    shutil.copy2(fig, args.report.parent / "figures" / fig.name)

    show = [
        "score_A",
        "fdr_A",
        "score_B",
        "fdr_B",
        "n_candidate_targets",
        "direction_agreement",
        "tf_gene_log2fc_A",
        "supporting_targets",
    ]
    pw_show = pw[["score_A", "fdr_A", "score_B", "fdr_B", "concordant"]].sort_values(
        "score_A", ascending=False
    )
    pw_show.index.name = "pathway"
    up = prioritized[prioritized["direction"] == "up_in_secondary"]
    dn = prioritized[prioritized["direction"] == "down_in_secondary"]
    for t in (up, dn):
        t.index.name = "tf"
    report = "\n".join(
        [
            f"# Phase 8 — transcription-factor and pathway activity · {SERIES}",
            "",
            "> Activity is inferred from target-gene expression (decoupler ULM on CollecTRI). It is **not** evidence "
            "of TF binding or function. Claim level: A/B boundary — supported by expression only.",
            "",
            f"- TFs scored: {len(tf):,}; concordant in both contrasts (FDR < {FDR}, same sign): "
            f"{int(tf['concordant'].sum())}; also target-supported (≥ {MIN_TARGETS} candidate targets, "
            f"≥ {MIN_AGREEMENT:.0%} direction agreement): **{len(prioritized)}**.",
            "",
            "## PHASE 8 GATE: PASS — every prioritised TF is cross-checked against target-gene support.",
            "",
            f"## Prioritised TFs — more active in secondary ({len(up)})",
            _md(up[show].head(20)),
            "",
            f"## Prioritised TFs — less active in secondary ({len(dn)})",
            _md(dn[show].head(20)),
            "",
            "## TFs named by the independent studies",
            _md(lit),
            "",
            "## PROGENy pathway activity",
            _md(pw_show),
            "",
            "## Figures",
            f"![{fig.stem}](figures/{fig.name})",
            "",
        ]
    )
    args.report.write_text(report, encoding="utf-8")
    print(f"[tf] {len(prioritized)} prioritised TFs · report → {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
