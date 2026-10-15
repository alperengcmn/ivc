"""Script 07 — Phase 7: pathway enrichment and the frozen discovery signature.

Pathways: pre-ranked GSEA (gseapy) on each Phase 6 contrast with MSigDB mouse Hallmark and
Reactome. A pathway is *concordant* when FDR < 0.05 in both contrasts with the same NES sign.

Signature (frozen here, before any validation-dataset expression values are analysed):
  core genes = Phase 6 concordant candidates that
    * keep bootstrap sign stability >= 0.95, and
    * remain concordant (same sign, FDR < 0.05 in both contrasts) in non-cycling ILC2s,
  so the signature is not simply a cell-cycle program that any activated dataset would match.
  The top SIGNATURE_SIZE genes per direction (ranked by the smaller |log2FC| of the two
  contrasts) form `recall_up` / `recall_down`.

The frozen file is written to signatures/ with a SHA-256 in signatures/MANIFEST.tsv and is
never edited afterwards; a changed signature gets a new version number.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
from collections.abc import Sequence
from datetime import datetime, timezone
from pathlib import Path

from .gse295590 import SERIES

FDR = 0.05
MIN_STABILITY = 0.95
SIGNATURE_SIZE = 50
VERSION = "v1"
GSEA_PERMUTATIONS = 1000
GENE_SETS = {
    "hallmark": "references/msigdb_hallmark_mouse.gmt",
    "reactome": "references/msigdb_reactome_mouse.gmt",
}


def read_gmt(path: Path) -> dict[str, list[str]]:
    sets = {}
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) > 2:
            sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def prerank(rnk, gene_sets: dict, seed: int):
    import gseapy as gp

    res = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        permutation_num=GSEA_PERMUTATIONS,
        min_size=15,
        max_size=500,
        seed=seed,
        threads=4,
        outdir=None,
        verbose=False,
    )
    df = res.res2d.set_index("Term")
    return (
        df[["NES", "FDR q-val", "Lead_genes"]]
        .rename(columns={"FDR q-val": "fdr"})
        .astype({"NES": float, "fdr": float})
    )


def concordant_pathways(a, b):
    import numpy as np

    terms = a.index.intersection(b.index)
    t = (
        a.loc[terms, ["NES", "fdr"]]
        .add_suffix("_A")
        .join(b.loc[terms, ["NES", "fdr"]].add_suffix("_B"))
    )
    t["concordant"] = (
        (t["fdr_A"] < FDR) & (t["fdr_B"] < FDR) & (np.sign(t["NES_A"]) == np.sign(t["NES_B"]))
    )
    t["direction"] = np.where(t["NES_A"] > 0, "up_in_secondary", "down_in_secondary")
    return t.sort_values("NES_A", ascending=False)


def freeze_signature(de, strat, out_dir: Path):
    import numpy as np
    import pandas as pd

    cand = de[de["candidate"]].copy()
    s = strat.reindex(cand.index)
    keep_strat = (
        (s["fdr_A"] < FDR)
        & (s["fdr_B"] < FDR)
        & (np.sign(s["pb_log2fc_A"]) == np.sign(cand["pb_log2fc_A"]))
        & (np.sign(s["pb_log2fc_B"]) == np.sign(cand["pb_log2fc_B"]))
    )
    cand["retained_noncycling"] = keep_strat.fillna(False).to_numpy()
    core = cand[(cand["bootstrap_stability"] >= MIN_STABILITY) & cand["retained_noncycling"]]
    rows = []
    for direction, name in (("up_in_secondary", "recall_up"), ("down_in_secondary", "recall_down")):
        top = (
            core[core["direction"] == direction]
            .sort_values("min_abs_lfc", ascending=False)
            .head(SIGNATURE_SIZE)
        )
        for rank, (gene, r) in enumerate(top.iterrows(), start=1):
            rows.append(
                {
                    "set": name,
                    "rank": rank,
                    "gene": gene,
                    "log2fc_A": round(r["pb_log2fc_A"], 3),
                    "log2fc_B": round(r["pb_log2fc_B"], 3),
                    "fdr_A": r["fdr_A"],
                    "fdr_B": r["fdr_B"],
                    "bootstrap_stability": r["bootstrap_stability"],
                }
            )
    sig = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{SERIES}_recall_ILC2_signature_{VERSION}.tsv"
    if path.exists():
        existing = pd.read_csv(path, sep="\t")
        if list(existing["gene"]) != list(sig["gene"]):
            raise RuntimeError(
                f"{path} is frozen and differs from the recomputed signature; bump VERSION."
            )
    else:
        sig.to_csv(path, sep="\t", index=False, float_format="%.4g")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = out_dir / "MANIFEST.tsv"
    if not manifest.exists():
        manifest.write_text("file\tversion\tsha256\tfrozen_utc\trule\n")
    if path.name not in manifest.read_text():
        stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        rule = (
            f"Phase 6 concordant (A & B), bootstrap >= {MIN_STABILITY}, retained in non-cycling ILC2s, "
            f"top {SIGNATURE_SIZE} per direction by min |log2FC|"
        )
        with manifest.open("a") as fh:
            fh.write(f"{path.name}\t{VERSION}\t{digest}\t{stamp}\t{rule}\n")
    return sig, cand, path, digest


def score_figure(ilc, sig, out_dir: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc

    ilc.X = ilc.layers["lognorm_dm"].copy()
    up = [g for g in sig.loc[sig["set"] == "recall_up", "gene"] if g in ilc.var_names]
    dn = [g for g in sig.loc[sig["set"] == "recall_down", "gene"] if g in ilc.var_names]
    sc.tl.score_genes(ilc, up, score_name="sig_up", random_state=0, use_raw=False)
    sc.tl.score_genes(ilc, dn, score_name="sig_down", random_state=0, use_raw=False)
    ilc.obs["recall_score"] = ilc.obs["sig_up"] - ilc.obs["sig_down"]
    libs = ["GSM8953514", "GSM8953515", "GSM8953516"]
    labels = ["S1 primary\n(depth-matched)", "S2 primary", "S3 secondary"]
    data = [
        ilc.obs.loc[ilc.obs["sample"].astype(str) == s, "recall_score"].to_numpy() for s in libs
    ]
    fig, ax = plt.subplots(figsize=(6, 4))
    parts = ax.violinplot(data, showmedians=True)
    for body, col in zip(parts["bodies"], ["#2563eb", "#60a5fa", "#d97706"], strict=True):
        body.set_facecolor(col)
        body.set_alpha(0.6)
    ax.set_xticks([1, 2, 3], labels)
    ax.set_ylabel("recall score (up − down)")
    ax.set_title(
        f"Frozen signature {VERSION} in discovery data\n(in-sample: separation is expected, not evidence)"
    )
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "signature_score_discovery.png"
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


def render(paths: dict, sig, cand, sig_path, digest, figs) -> str:
    lines = [f"# Phase 7 — pathways and frozen discovery signature · {SERIES}", ""]
    for name, t in paths.items():
        conc = t[t["concordant"]]
        lines += [
            f"## {name.capitalize()} — concordant pathways ({len(conc)} of {len(t)} tested)",
            _md(conc.drop(columns=["concordant"]).head(25)),
            "",
        ]
    n_core = int(
        ((cand["bootstrap_stability"] >= MIN_STABILITY) & cand["retained_noncycling"]).sum()
    )
    up = sig.loc[sig["set"] == "recall_up", "gene"].tolist()
    dn = sig.loc[sig["set"] == "recall_down", "gene"].tolist()
    lines += [
        f"## Frozen signature `{sig_path.name}`",
        f"- SHA-256: `{digest}`",
        f"- Phase 6 candidates: {len(cand):,} → retained in non-cycling ILC2s and bootstrap ≥ {MIN_STABILITY}: "
        f"**{n_core:,}** → top {SIGNATURE_SIZE} per direction.",
        f"- **recall_up ({len(up)})**: {', '.join(up)}",
        f"- **recall_down ({len(dn)})**: {', '.join(dn)}",
        "",
        "Frozen before any validation-dataset expression values were analysed. Phase 10 tests it in independent "
        "datasets; the discovery-data score below is in-sample and is shown only as a sanity check.",
        "",
        "## Phase 7 gate",
        "The pipeline gate (signature stable on ≥ 1 compatible dataset) is evaluated in Phase 10. Internal "
        "stability is shown by the bootstrap and non-cycling filters above.",
        "",
        "## Figures",
        *[f"![{p.stem}](figures/{p.name})" for p in figs],
        "",
    ]
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--de-dir", type=Path, default=Path("results/differential_expression"))
    p.add_argument(
        "--annotated", type=Path, default=Path("results/annotation") / f"{SERIES}_annotated.h5ad"
    )
    p.add_argument("--out-dir", type=Path, default=Path("results/pathways"))
    p.add_argument("--signature-dir", type=Path, default=Path("signatures"))
    p.add_argument(
        "--report", type=Path, default=Path("reports/phase7/07_pathways_signature_report.md")
    )
    p.add_argument("--params", type=Path, default=None)
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    import pandas as pd
    import scanpy as sc

    from .config import parameters

    for f in GENE_SETS.values():
        if not Path(f).is_file():
            print(f"Error: {f} missing — run `make fetch-resources`.", file=sys.stderr)
            return 2
    seed = int(parameters(args.params).get("random_seed", 0))
    de = pd.read_csv(args.de_dir / "differential_expression_results.csv", index_col=0)
    strat = pd.read_csv(args.de_dir / "noncycling_stratum.csv", index_col=0)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    pathways = {}
    for name, gmt in GENE_SETS.items():
        sets = read_gmt(Path(gmt))
        res = {}
        for c in ("A", "B"):
            rnk = pd.read_csv(
                args.de_dir / f"ranked_{c}.rnk", sep="\t", header=None, index_col=0
            ).iloc[:, 0]
            rnk = rnk[~rnk.index.duplicated()].dropna()
            res[c] = prerank(rnk, sets, seed)
            res[c].to_csv(args.out_dir / f"gsea_{name}_{c}.csv")
        pathways[name] = concordant_pathways(res["A"], res["B"])
        pathways[name].to_csv(args.out_dir / f"pathway_enrichment_{name}.csv")
        print(f"[pathways] {name}: {int(pathways[name]['concordant'].sum())} concordant")

    sig, cand, sig_path, digest = freeze_signature(de, strat, args.signature_dir)
    cand.to_csv(args.out_dir / "candidate_signature_genes.tsv", sep="\t")
    print(f"[pathways] signature frozen: {sig_path} ({len(sig)} genes)")

    adata = sc.read_h5ad(args.annotated)
    ilc = adata[adata.obs["lineage"].astype(str) == "ILC2"].copy()
    figs = [score_figure(ilc, sig, args.out_dir / "figures")]
    (args.report.parent / "figures").mkdir(parents=True, exist_ok=True)
    for f in figs:
        shutil.copy2(f, args.report.parent / "figures" / f.name)
    args.report.write_text(render(pathways, sig, cand, sig_path, digest, figs), encoding="utf-8")
    print(f"[pathways] report → {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
