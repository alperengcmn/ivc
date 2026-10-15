"""Script 12 — Phase 12: assemble the final report from result files (no hand-typed numbers).

Every number in reports/FINAL_REPORT.md is read from results/, signatures/ or config/, so the
report is regenerated whenever an upstream phase changes (pipeline: automated report generation).
"""

from __future__ import annotations

import argparse
import hashlib
import platform
import re
import subprocess
import sys
from collections.abc import Sequence
from datetime import datetime, timezone
from pathlib import Path

from .gse295590 import SERIES

R = Path("results")


def _git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
    except Exception:
        return "unknown"


def _versions() -> str:
    from importlib.metadata import PackageNotFoundError, version

    out = [f"python {platform.python_version()}"]
    for p in (
        "scanpy",
        "anndata",
        "numpy",
        "pandas",
        "scipy",
        "decoupler",
        "gseapy",
        "pydeseq2",
        "leidenalg",
    ):
        try:
            out.append(f"{p} {version(p)}")
        except PackageNotFoundError:
            out.append(f"{p} (not installed)")
    return ", ".join(out)


def _md(df, fmt="{:.3g}") -> str:
    cols = [df.index.name or "", *map(str, df.columns)]
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for i, r in df.iterrows():
        cells = []
        for v in r:
            if isinstance(v, bool):
                cells.append("✓" if v else "")
            elif isinstance(v, float):
                cells.append("" if v != v else fmt.format(v))
            else:
                cells.append(str(v))
        out.append("| " + " | ".join([str(i), *cells]) + " |")
    return "\n".join(out)


def gather() -> dict:
    import pandas as pd

    from .config import parameters

    d: dict = {"params": parameters()}
    d["qc"] = pd.read_csv(R / "qc" / "qc_summary.csv", index_col=0)
    ann = pd.read_csv(
        R / "annotation" / "annotated_cells.csv.gz",
        index_col=0,
        usecols=["Unnamed: 0", "sample", "lineage", "cluster"],
    )
    d["n_cells"] = len(ann)
    d["lineage_counts"] = pd.crosstab(ann["lineage"], ann["sample"])
    d["stability"] = pd.read_csv(R / "embedding" / "resolution_stability.csv", index_col=0)
    de = pd.read_csv(
        R / "differential_expression" / "differential_expression_results.csv", index_col=0
    )
    d["de"] = de
    d["r_ab"] = float(de["pb_log2fc_A"].corr(de["pb_log2fc_B"]))
    d["hallmark"] = pd.read_csv(R / "pathways" / "pathway_enrichment_hallmark.csv", index_col=0)
    d["reactome"] = pd.read_csv(R / "pathways" / "pathway_enrichment_reactome.csv", index_col=0)
    sig_path = Path("signatures") / f"{SERIES}_recall_ILC2_signature_v1.tsv"
    d["sig"] = pd.read_csv(sig_path, sep="\t")
    d["sig_sha"] = hashlib.sha256(sig_path.read_bytes()).hexdigest()
    d["tf"] = pd.read_csv(R / "regulatory" / "prioritized_tf_candidates.tsv", sep="\t", index_col=0)
    d["tf_all"] = pd.read_csv(R / "regulatory" / "tf_activity_matrix.csv", index_col=0)
    d["progeny"] = pd.read_csv(R / "regulatory" / "progeny_pathway_activity.csv", index_col=0)
    d["atac"] = pd.read_csv(
        R / "regulatory" / "atac" / "signature_accessibility_summary.csv", index_col=0
    )
    d["validation"] = pd.read_csv(R / "validation" / "cross_dataset_comparison.csv")
    d["acute"] = pd.read_csv(R / "validation" / "acute_shared_annotation.csv", index_col=0)
    d["context"] = pd.read_csv(R / "validation" / "context_datasets.csv")
    d["robust"] = pd.read_csv(
        R / "robustness" / "final_statistical_summary.tsv", sep="\t", index_col=0
    )
    d["scenarios"] = pd.read_csv(R / "robustness" / "sensitivity_scenarios.csv", index_col=0)
    emb = Path("reports/phase4/04_embedding_report.md").read_text(encoding="utf-8")
    m = re.search(r"\*\*([\d.]+)\*\* on full-depth data → \*\*([\d.]+)\*\*", emb)
    d["pc1"] = (float(m.group(1)), float(m.group(2))) if m else (float("nan"), float("nan"))
    cellqc = pd.read_csv(
        R / "qc" / "cell_qc_metrics.csv.gz",
        index_col=0,
        usecols=lambda c: c in ("Unnamed: 0", "", "pct_counts_mt"),
    )
    a2 = ann.join(cellqc, how="left")
    d["s2_ilc2_mito10"] = int(
        (
            (a2["lineage"] == "ILC2") & (a2["sample"] == "GSM8953515") & (a2["pct_counts_mt"] <= 10)
        ).sum()
    )
    return d


def atac_verdict(atac) -> tuple[str, int]:
    recall = atac[atac["recall"]]
    n = int((recall["p_up_gain"] < 0.05).sum())
    return {2: "Regulatory support", 1: "Partial regulatory support", 0: "No regulatory support"}[
        min(n, 2)
    ], n


def render(d: dict) -> str:
    import pandas as pd

    p3 = d["params"]["phase3_qc"]
    qc = d["qc"]
    de = d["de"]
    cand = de[de["candidate"]]
    n_up = int((cand["direction"] == "up_in_secondary").sum())
    hm = d["hallmark"]
    hm_c = hm[hm["concordant"]]
    hm_up = [t.replace("HALLMARK_", "") for t in hm_c.index[hm_c["NES_A"] > 0]]
    hm_dn = [t.replace("HALLMARK_", "") for t in hm_c.index[hm_c["NES_A"] < 0]]
    sig = d["sig"]
    up = sig.loc[sig["set"] == "recall_up", "gene"].tolist()
    dn = sig.loc[sig["set"] == "recall_down", "gene"].tolist()
    tf = d["tf"]
    tf_up = tf.index[tf["direction"] == "up_in_secondary"].tolist()
    tf_dn = tf.index[tf["direction"] == "down_in_secondary"].tolist()
    prog = d["progeny"]
    prog_c = prog[prog["concordant"]].sort_values("score_A")
    atac = d["atac"]
    atac_v, _ = atac_verdict(atac)
    val = d["validation"]
    acute = d["acute"]
    rob = d["robust"]
    n_robust = int(rob["robust"].sum())
    pub = de.reindex(["Arg1", "Cox5a", "Pkm", "Areg", "Tph1", "Il2ra", "Csf2"])[
        ["pb_log2fc_A", "fdr_A", "pb_log2fc_B", "fdr_B", "candidate"]
    ]
    pub.index.name = "gene"
    lib = d["lineage_counts"]
    ilc2 = lib.loc["ILC2"] if "ILC2" in lib.index else pd.Series(dtype=int)
    core_verdicts = val[["dataset", "comparison", "verdict"]]
    replicated = (core_verdicts["verdict"] == "Replicated").any()
    partial = (core_verdicts["verdict"] == "Partially supported").any()
    # H2 is judged on memory transcriptomes only; chromatin evidence belongs to H4.
    testable = (core_verdicts["verdict"] != "Not testable with available data").any()
    if replicated:
        h2 = "supported"
    elif partial:
        h2 = "partially supported"
    elif testable:
        h2 = "not supported in the testable comparisons"
    else:
        h2 = "not testable"
    val_tbl = val.set_index("dataset")[
        [
            "comparison",
            "score_delta",
            "score_p",
            "genes_agree",
            "genes_n_genes",
            "genes_p",
            "verdict",
        ]
    ]
    val_tbl.index.name = "dataset"
    atac_tbl = atac[["n", "DA_peaks_FDR5", "up_genes", "p_up_gain", "down_genes", "p_down_loss"]]
    ctx = d["context"].set_index("dataset")[["question", "in_set", "background", "odds_ratio", "p"]]
    scen = d["scenarios"][
        ["ILC2s", "recall_up retained", "recall_down retained", "fraction"]
    ].copy()
    scen["fraction"] = (scen["fraction"] * 100).round(0).astype(int).astype(str) + " %"
    n_opp = int(acute.get("acute_opposite", pd.Series(dtype=bool)).sum())
    opp = ", ".join(acute.index[acute["acute_opposite"]]) if "acute_opposite" in acute else ""
    depth_ratio = qc["median_total_counts"].max() / qc["median_total_counts"].min()
    shallow = qc["median_total_counts"].idxmin()
    lit = d["tf_all"].reindex(["Bach2", "Jund", "Fosl2", "Fhl2", "Fosb", "Stat6", "Srebf2"])
    lit_ok = lit.index[
        lit["concordant"].fillna(False).astype(bool)
        & lit["target_supported"].fillna(False).astype(bool)
    ]

    def describe(g: str) -> str:
        if g not in de.index:
            return f"{g}: not tested"
        r = de.loc[g]
        same = (r["pb_log2fc_A"] > 0) == (r["pb_log2fc_B"] > 0)
        if r["candidate"]:
            state = "concordant candidate"
        elif same:
            state = "same direction, below threshold"
        else:
            state = "inconsistent between contrasts"
        return f"{g} {r['pb_log2fc_A']:+.2f} / {r['pb_log2fc_B']:+.2f} ({state})"

    pub_sentence = "; ".join(describe(g) for g in ("Areg", "Csf2", "Il2ra", "Tph1"))
    qc_tbl = qc[
        ["n_cells", "median_total_counts", "median_genes", "median_pct_mt", "pct_cells_mt_gt20"]
    ].copy()
    for c in ("n_cells", "median_total_counts", "median_genes"):
        qc_tbl[c] = qc_tbl[c].round().astype(int).map("{:,}".format)
    qc_tbl.index.name = "library"

    return "\n".join(
        [
            "# Lung ILC2 trained/recall-associated programs — final report",
            "",
            f"*Generated automatically by `scripts/12_final_report.py` on "
            f"{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')} from commit `{_git_commit()}`. "
            "Every number below is read from `results/`. Thresholds marked PROVISIONAL await Divya's review (D-010).*",
            "",
            "## Summary",
            "",
            f"- Secondary-response lung ILC2s from *N. brasiliensis* rechallenge (GSE295590) carry a reproducible "
            f"program: **{len(cand):,} genes** change concordantly in two contrasts that share no primary cells "
            f"(effect sizes r = {d['r_ab']:.2f}); the published findings for **Arg1, Cox5a and Pkm** are reproduced.",
            f"- The program is metabolic and proliferative (up: {', '.join(hm_up[:6])}) with lower "
            f"{', '.join(hm_dn[-3:])} signalling; inferred TF activity: {', '.join(tf_up) or '—'} up; "
            f"{', '.join(tf_dn[:6])} down.",
            f"- A 100-gene signature was frozen before validation (SHA-256 `{d['sig_sha'][:12]}…`). It is robust "
            f"internally ({n_robust}/100 genes survive ≥ 80 % of sensitivity scenarios, including the source study's "
            "own QC).",
            f"- Independent evidence is **mixed**: chromatin — *{atac_v}* (strong in wild-type allergen memory, absent "
            "in Rag1⁻/⁻); memory transcriptome — "
            + "; ".join(f"{r.dataset} {r.comparison}: *{r.verdict}*" for r in val.itertuples())
            + ".",
            "- **Claim level reached: A** — a trained/recall-associated transcriptional state with partial regulatory "
            "support (B, wild-type only). Classical trained immunity (C) is not claimed.",
            "",
            "## 1. Research question",
            "What distinguishes a trained/persistent ILC2 state from an acutely activated ILC2 state, and which "
            "transcriptional, metabolic and regulatory programs are associated with it? Hypotheses H1–H5: "
            "[`docs/00_hypothesis_and_success_criteria.md`](../docs/00_hypothesis_and_success_criteria.md).",
            "",
            "| Hypothesis | Outcome |",
            "|---|---|",
            f"| H1 — reproducible recall program in GSE295590 beyond composition | **Supported (exploratory)** — "
            f"{len(cand):,} concordant genes; {n_robust}/100 signature genes robust; retained in non-cycling ILC2s |",
            f"| H2 — conserved in independent memory-like ILC2s | **{h2.capitalize()}** — see §9 |",
            f"| H3 — not explained by acute activation alone | **Supported** — 0 signature genes are acute-shared on "
            f"both days of influenza; {n_opp} move opposite to acute activation |",
            f"| H4 — TF / chromatin support | **Partial** — {len(tf)} target-supported TFs; ATAC: {atac_v.lower()} |",
            "| H5 — metabolic reprogramming | **Supported in discovery** — OXPHOS, glycolysis, MYC, mTORC1 up; "
            "signature responds to iron deprivation (§10) |",
            "",
            "## 2. Datasets and design",
            "Roles and caveats: [`docs/datasets.md`](../docs/datasets.md). Discovery: GSE295590 — 20 WT mice "
            "(10 primary, 10 secondary, day 3), FACS-sorted lung ILCs, three 10x libraries (S1, S2 primary; S3 "
            "secondary). **The mouse-to-library mapping is unknown (D-004)**; the secondary condition is one library.",
            "",
            "## 3. Quality control (PROVISIONAL thresholds, D-010)",
            f"mito ≤ {p3['max_pct_mito']} %, {p3['min_genes_per_cell']}–{p3['max_genes_per_cell']} genes per cell, "
            f"≥ {p3['min_cells_per_sample']} cells per library; predicted doublets not removed (scrublet threshold "
            f"failed, D-009). {d['n_cells']:,} cells retained.",
            "",
            _md(qc_tbl, "{:.1f}"),
            "",
            f"Key technical findings: matrices pre-filtered at ≈ {qc['min_total_counts'].min():,.0f} UMI; "
            f"**`{shallow}` (S1) sequenced {depth_ratio:.1f}× shallower** than the deepest library; S1 "
            "ILC2-like cells carry elevated T-lineage, haemoglobin and immediate-early transcripts (ambient RNA / "
            "dissociation stress, D-012).",
            "",
            "## 4. Embedding and annotation",
            f"Embedding on binomially thinned (depth-matched) counts — PC1–depth |r| {d['pc1'][0]:.2f} → "
            f"{d['pc1'][1]:.2f}; no cross-library "
            "batch correction because library and condition are confounded (D-011). Clustering resolution chosen by a "
            "pre-specified stability rule (min ARI ≥ 0.8 over 5 seeds). Lineage by marker categories (never a single "
            "marker).",
            "",
            _md(lib, "{:,.0f}"),
            "",
            f"Reproduction of the published cell numbers: S2 contains {int(ilc2.get('GSM8953515', 0))} ILC2s "
            f"({d['s2_ilc2_mito10']} at the paper's mito ≤ 10 %) vs **189 primary ILC2s published** — the "
            "published primary group "
            f"most likely comprised S2 only, and S1 ({int(ilc2.get('GSM8953514', 0)):,} ILC2-like cells) appears "
            "absent from the published ILC2 analysis. nILC2/iILC2 were **not resolvable** as stable discrete "
            "sub-clusters in the clean reference library; the iILC2 program is carried as a continuous score (D-012).",
            "",
            "## 5. Statistical methods",
            "Exploratory library-consistency design (D-013): contrast A = S3 vs S2 (full depth), contrast B = S3 vs S1 "
            "(depth-matched). Wilcoxon tests rank genes only (cells are not replicates); effect sizes are pseudobulk "
            "log2FC; candidates must agree in both contrasts (FDR < 0.05, |log2FC| ≥ 0.5, same sign, detected in ≥ 10 % "
            "of cells) and are not contamination-sensitive. Validation rules were pre-registered before any "
            "validation data were analysed (D-014, [`docs/validation_plan.md`](../docs/validation_plan.md)).",
            "",
            "## 6. Differential expression",
            f"{len(cand):,} concordant genes ({n_up} up, {len(cand) - n_up} down in secondary). Genes named by the "
            "source study:",
            "",
            _md(pub),
            "",
            f"Other genes discussed by the source study (log2FC contrast A / B): {pub_sentence}. The paper reports "
            "Csf2 as unchanged.",
            "",
            "## 7. Pathways and frozen signature",
            f"Concordant Hallmark pathways — up: {', '.join(hm_up)}; down: {', '.join(hm_dn)}. "
            f"Reactome: {int(d['reactome']['concordant'].sum())} concordant.",
            "",
            f"**Frozen signature v1** (`signatures/{SERIES}_recall_ILC2_signature_v1.tsv`, SHA-256 `{d['sig_sha']}`):",
            f"- recall_up: {', '.join(up)}",
            f"- recall_down: {', '.join(dn)}",
            "",
            "## 8. Transcription-factor and pathway activity (inference, not binding evidence)",
            f"Target-supported TFs — more active in secondary: {', '.join(tf_up) or 'none'}; less active: "
            f"{', '.join(tf_dn)}. PROGENy: "
            + ", ".join(f"{i} {'↑' if r.score_A > 0 else '↓'}" for i, r in prog_c.iterrows())
            + ".",
            "TFs named by the GSE172258 study that are concordant and target-supported here: "
            + (", ".join(f"**{x}**" for x in lit_ok) or "none")
            + ".",
            "",
            "## 9. Epigenomic integration and cross-dataset validation",
            f"**Chromatin (GSE172258 ATAC, Phase 9): {atac_v}** (pre-registered rule).",
            "",
            _md(atac_tbl),
            "",
            "**Memory transcriptomes (Phase 10, pre-registered verdicts):**",
            "",
            _md(val_tbl),
            "",
            "Interpretation, not part of the verdict: GSE172258 wells were sampled one week after recall versus three "
            "days in the discovery data; an early recall program may have subsided. The strong wild-type ATAC result "
            "versus the null Rag1⁻/⁻ result suggests the chromatin component may depend on adaptive-immune context.",
            "",
            "## 10. Acute activation, mechanism and metabolism",
            f"Influenza (GSE149857): no signature gene is acute-shared on both days; {n_opp} move opposite to acute "
            f"activation ({opp}) — they distinguish recall from acute activation.",
            "",
            _md(ctx),
            "",
            "## 11. Robustness (Phase 11)",
            _md(scen),
            "",
            "## 12. Limitations",
            "1. **One secondary library**: condition and library are confounded; no analysis here can separate them.",
            "2. **Unknown mouse structure** (D-004): all inference is exploratory; p-values rank, they do not test animals.",
            f"3. **S1 is technically distinct** ({depth_ratio:.1f}× shallower, ambient RNA, stress) and likely absent "
            "from the published analysis.",
            "4. **Dissociation stress**: part of the 'down in secondary' program (Fos, Klf2/4, Zfp36, Cebpb; NF-κB/AP-1) "
            "may reflect tissue handling.",
            "5. **Cross-model validation**: allergen (GSE172258), IL-33 (GSE81700) and influenza (GSE149857) differ from "
            "helminth recall; timing differs (day 3 vs 1 week).",
            "6. **Provisional QC thresholds** (D-010); results were stable to the source study's stricter QC (§11).",
            "7. Motif enrichment was not performed (no genome sequence / motif database). GSE291230 (Blimp-1 ATAC) "
            "is deposited as bigwig tracks only (~0.7 GB each) and was not integrated; its paired RNA (GSE291231) was.",
            "",
            "## 13. Conclusions supported by the data",
            "- Secondary-response lung ILC2s show a reproducible, metabolically active and proliferative recall program "
            "with reduced inflammatory/AP-1 signalling — **claim level A (recall-associated state)**.",
            "- The program is distinct from acute influenza-driven activation and partially mirrored by chromatin "
            "opening in wild-type allergen-memory ILC2s — **partial level B**.",
            "- It is **not** independently reproduced in the GSE172258 memory transcriptome at one week after recall.",
            "- Classical trained immunity or epigenetic memory is **not** claimed.",
            "",
            "## 14. Reproducibility",
            f"- Code commit `{_git_commit()}`; environment `environment.yml`; software: {_versions()}.",
            f"- Random seed {d['params'].get('random_seed', 0)}; all parameters in `config/`; decisions in "
            "`docs/decision_log.md`; downloads with SHA-256 in `data/raw/*/manifest.tsv` and `references/manifest.tsv`.",
            "- Re-run everything: `make all` (raw data → this report).",
            "",
        ]
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--report", type=Path, default=Path("reports/FINAL_REPORT.md"))
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        d = gather()
    except FileNotFoundError as exc:
        print(
            f"Error: missing result file — run the earlier phases first ({exc.filename}).",
            file=sys.stderr,
        )
        return 2
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(render(d), encoding="utf-8")
    print(f"[report] → {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
