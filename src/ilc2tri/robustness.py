"""Script 11 — Phase 11: statistical robustness of the frozen signature.

Each scenario re-runs both Phase 6 contrasts on a subset of ILC2s (or on depth-matched data)
and asks how many frozen signature genes are *retained*: same sign as frozen, FDR < 0.05 and
|pseudobulk log2FC| >= 0.5 in both contrasts. Leave-one-cluster-out uses pseudobulk effects
only (sign and |log2FC| >= 0.5 in both contrasts).

A gene is *robust* when retained in >= ROBUST_FRACTION of all scenarios. Genes from the
dissociation-stress program (van den Brink et al., 2017) are flagged: their difference
between libraries may reflect tissue handling rather than biology.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from .gse295590 import SERIES

ROBUST_FRACTION = 0.8
SIGNATURE = Path("signatures/GSE295590_recall_ILC2_signature_v1.tsv")
# Core dissociation-induced genes (van den Brink et al., Nat Methods 2017) — immediate-early / heat-shock.
DISSOCIATION = {
    "Fos",
    "Fosb",
    "Jun",
    "Junb",
    "Jund",
    "Egr1",
    "Atf3",
    "Ier2",
    "Ier3",
    "Hspa1a",
    "Hspa1b",
    "Dusp1",
    "Zfp36",
    "Nr4a1",
    "Klf2",
    "Klf4",
    "Klf6",
    "Socs3",
    "Btg2",
    "Cebpb",
    "Cebpd",
}


def retained(table, sig) -> pd.Series:  # noqa: F821
    import numpy as np

    from .de import FDR, MIN_LFC

    t = table.reindex(sig["gene"])
    exp = np.where(sig["set"].to_numpy() == "recall_up", 1, -1)
    ok = (
        (np.sign(t["pb_log2fc_A"]).to_numpy() == exp)
        & (np.sign(t["pb_log2fc_B"]).to_numpy() == exp)
        & (t["fdr_A"].to_numpy() < FDR)
        & (t["fdr_B"].to_numpy() < FDR)
        & (t["pb_log2fc_A"].abs().to_numpy() >= MIN_LFC)
        & (t["pb_log2fc_B"].abs().to_numpy() >= MIN_LFC)
    )
    import pandas as pd

    return pd.Series(np.nan_to_num(ok, nan=0).astype(bool), index=sig["gene"].to_numpy())


def scenarios(ilc):
    """name -> (description, boolean cell mask or None, contrast override or None)."""
    import numpy as np

    obs = ilc.obs
    lib = obs["sample"].astype(str)
    stress_cut = obs.groupby(lib)["state_stress_response"].transform(lambda s: s.quantile(0.75))
    contam_cut = obs.groupby(lib)["contamination_score"].transform("median")
    return {
        "mito_10": ("mito ≤ 10 % (source-study QC)", (obs["pct_counts_mt"] <= 10).to_numpy(), None),
        "genes_500": ("≥ 500 genes per cell", (obs["n_genes_by_counts"] >= 500).to_numpy(), None),
        "noncycling": (
            "non-cycling ILC2s (proliferation score ≤ 0)",
            (obs["state_proliferation"] <= 0).to_numpy(),
            None,
        ),
        "low_stress": (
            "dissociation-stress score ≤ 75th percentile per library",
            (obs["state_stress_response"] <= stress_cut).to_numpy(),
            None,
        ),
        "low_contamination": (
            "contamination score ≤ library median",
            (obs["contamination_score"] <= contam_cut).to_numpy(),
            None,
        ),
        "A_depth_matched": (
            "contrast A on depth-matched data",
            np.ones(len(obs), dtype=bool),
            "A_dm",
        ),
    }


def leave_one_cluster_out(ilc, sig):
    import numpy as np
    import pandas as pd

    from .de import CONTRASTS, MIN_LFC, pseudobulk_lfc

    lib = ilc.obs["sample"].astype(str).to_numpy()
    cl = ilc.obs["cluster"].astype(str).to_numpy()
    genes = [g for g in sig["gene"] if g in ilc.var_names]
    gi = [ilc.var_names.get_loc(g) for g in genes]
    exp = np.where(sig.set_index("gene").loc[genes, "set"].to_numpy() == "recall_up", 1, -1)
    rows = {}
    for c in sorted(set(cl), key=lambda x: int(x) if x.isdigit() else x):
        keep = cl != c
        ok = np.ones(len(genes), dtype=bool)
        for con in CONTRASTS:
            X = ilc.layers[con.counts]
            s = keep & (lib == con.secondary)
            p = keep & (lib == con.primary)
            if s.sum() < 20 or p.sum() < 20:
                ok &= False
                continue
            lfc_all = pseudobulk_lfc(X, s, p)
            lfc = lfc_all[gi]
            ok &= (np.sign(lfc) == exp) & (np.abs(lfc) >= MIN_LFC)
        rows[f"drop_cluster_{c}"] = pd.Series(ok, index=genes)
    return pd.DataFrame(rows).reindex(sig["gene"]).fillna(False)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--input", type=Path, default=Path("results/annotation") / f"{SERIES}_annotated.h5ad"
    )
    p.add_argument("--signature", type=Path, default=SIGNATURE)
    p.add_argument("--markers", type=Path, default=None)
    p.add_argument("--out-dir", type=Path, default=Path("results/robustness"))
    p.add_argument("--report", type=Path, default=Path("reports/phase11/11_robustness_report.md"))
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    import pandas as pd
    import scanpy as sc

    from .annotate import load_markers
    from .de import CONTRASTS, Contrast, concordant, excluded_genes, run_contrast

    if not args.input.is_file() or not args.signature.is_file():
        print("Error: run Phases 5–7 first.", file=sys.stderr)
        return 2
    sig = pd.read_csv(args.signature, sep="\t")
    adata = sc.read_h5ad(args.input)
    ilc = adata[adata.obs["lineage"].astype(str) == "ILC2"].copy()
    del adata
    excl = excluded_genes(
        ilc.var_names, load_markers(args.markers).get("contamination_sensitive", [])
    )
    a_dm = Contrast(
        "A", "S3 vs S2 (depth-matched)", "GSM8953516", "GSM8953515", "lognorm_dm", "depth_matched"
    )

    cols, summary = {}, []
    for key, (label, mask, override) in scenarios(ilc).items():
        cons = [a_dm, CONTRASTS[1]] if override == "A_dm" else list(CONTRASTS)
        ta, tb = (run_contrast(ilc, c, mask) for c in cons)
        conc = concordant(ta, tb, excl)
        r = retained(conc, sig)
        cols[key] = r
        n_cells = int(mask.sum())
        summary.append(
            {
                "scenario": label,
                "ILC2s": n_cells,
                "candidates": int(conc["candidate"].sum()),
                "recall_up retained": f"{int(r[sig['set'].eq('recall_up').to_numpy()].sum())}/50",
                "recall_down retained": f"{int(r[sig['set'].eq('recall_down').to_numpy()].sum())}/50",
                "fraction": float(r.mean()),
            }
        )
        print(f"[robust] {key}: {r.mean():.0%} of signature retained")
    loco = leave_one_cluster_out(ilc, sig)
    for c in loco.columns:
        summary.append(
            {
                "scenario": f"leave out {c.replace('drop_cluster_', 'cluster ')} (pseudobulk only)",
                "ILC2s": int((ilc.obs["cluster"].astype(str) != c.split("_")[-1]).sum()),
                "candidates": float("nan"),
                "recall_up retained": f"{int(loco[c][sig['set'].eq('recall_up').to_numpy()].sum())}/50",
                "recall_down retained": f"{int(loco[c][sig['set'].eq('recall_down').to_numpy()].sum())}/50",
                "fraction": float(loco[c].mean()),
            }
        )
    allm = pd.concat([pd.DataFrame(cols), loco.set_axis(pd.DataFrame(cols).index)], axis=1)
    per_gene = sig.set_index("gene")[["set"]].copy()
    per_gene["fraction_retained"] = allm.mean(axis=1).to_numpy()
    per_gene["robust"] = per_gene["fraction_retained"] >= ROBUST_FRACTION
    per_gene["dissociation_gene"] = per_gene.index.isin(DISSOCIATION)
    summ = pd.DataFrame(summary).set_index("scenario")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    allm.to_csv(args.out_dir / "robustness_analysis.csv")
    per_gene.to_csv(args.out_dir / "final_statistical_summary.tsv", sep="\t")
    summ.to_csv(args.out_dir / "sensitivity_scenarios.csv")

    robust_up = per_gene[(per_gene["set"] == "recall_up") & per_gene["robust"]].index.tolist()
    robust_dn = per_gene[(per_gene["set"] == "recall_down") & per_gene["robust"]].index.tolist()
    diss = per_gene[per_gene["dissociation_gene"]].index.tolist()
    n_dn = int((per_gene.loc[diss, "set"] == "recall_down").sum())
    diss_dir = f"{n_dn} of {len(diss)} in recall_down, i.e. higher in the primary libraries"

    def md(df):
        cols_ = [df.index.name or "", *map(str, df.columns)]
        out = ["| " + " | ".join(cols_) + " |", "|" + "---|" * len(cols_)]
        for i, r in df.iterrows():
            out.append(
                "| "
                + " | ".join(
                    [str(i)]
                    + [
                        ""
                        if (isinstance(v, float) and v != v)
                        else (
                            f"{v:.0%}"
                            if k == "fraction"
                            else (f"{v:,}" if isinstance(v, int) else str(v))
                        )
                        for k, v in r.items()
                    ]
                )
                + " |"
            )
        return "\n".join(out)

    report = "\n".join(
        [
            f"# Phase 11 — statistical robustness · {SERIES}",
            "",
            "## PHASE 11 GATE: PASS — all pre-specified sensitivity analyses completed.",
            "",
            "A signature gene is *retained* in a scenario when it keeps its frozen sign, FDR < 0.05 and |log2FC| ≥ 0.5 in "
            f"both contrasts (leave-one-cluster-out: sign and |log2FC| only). *Robust* = retained in ≥ "
            f"{ROBUST_FRACTION:.0%} of {allm.shape[1]} scenarios.",
            "",
            md(summ),
            "",
            "The non-cycling scenario is partly circular (retention in non-cycling ILC2s was a selection criterion "
            "for the frozen signature, Phase 7) and is shown for completeness only.",
            "",
            f"## Robust genes: {len(robust_up) + len(robust_dn)} of 100",
            f"- **recall_up ({len(robust_up)})**: {', '.join(robust_up)}",
            f"- **recall_down ({len(robust_dn)})**: {', '.join(robust_dn)}",
            "",
            "## Dissociation-stress flag",
            f"{len(diss)} signature genes belong to the dissociation-induced program (van den Brink et al., 2017): "
            f"{', '.join(diss)} ({diss_dir}). The primary libraries carry higher stress scores (Phase 5), so these "
            "differences may reflect tissue handling between libraries; interpret the NF-κB/AP-1/TNFα 'down in "
            "secondary' results (Phases 7–8) with this caveat.",
            "",
            "## What these analyses cannot address",
            "- Library and condition are confounded (one secondary library); no scenario can separate them.",
            "- The mouse-to-library structure remains unknown (D-004), so no animal-level inference is possible.",
            "",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(report, encoding="utf-8")
    print(
        f"[robust] robust genes: {len(robust_up)} up, {len(robust_dn)} down · report → {args.report}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
