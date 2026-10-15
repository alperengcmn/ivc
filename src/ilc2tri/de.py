"""Script 06 — Phase 6: exploratory, library-consistency differential expression in ILC2s.

Why not a standard replicate-based test: the mouse-to-library mapping of GSE295590 is
unknown (D-004), and the secondary condition is a single library. Cell-level p-values
therefore treat cells as replicates (pseudoreplication) and are used for *ranking only*.

What makes a gene a candidate instead — it must agree across two independent contrasts
that share no primary cells:

  A. secondary S3 vs primary S2 — full depth (similar depth: ~8.5k vs ~7.4k UMI/cell)
  B. secondary S3 vs primary S1 — depth-matched (S1 is ~5.6x shallower; D-009)

Candidate = FDR < 0.05 in both, |pseudobulk log2FC| >= 0.5 in both, same sign, detected
in >= 10 % of cells of the higher group, and not contamination-sensitive (D-012).

Robustness reported here: (i) bootstrap sign stability of the pseudobulk effect,
(ii) retention when proliferating cells are excluded (composition check).
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from .gse295590 import SERIES

FDR = 0.05
MIN_LFC = 0.5
MIN_DETECT = 0.10
MIN_DETECT_TEST = (
    0.05  # genes detected in < 5 % of both groups are not tested (as in the source study)
)
N_BOOT = 200
PROLIF_MAX = 0.0  # non-cycling stratum: proliferation score <= 0
PUBLISHED_GENES = [
    "Arg1",
    "Cox5a",
    "Pkm",
    "Areg",
    "Tph1",
    "Il2ra",
    "Il1rl1",
    "Mki67",
    "Il7r",
    "Csf2",
]


@dataclass(frozen=True)
class Contrast:
    name: str
    label: str
    secondary: str
    primary: str
    layer: str  # log-normalised layer used for the test
    counts: str  # count layer used for pseudobulk effect sizes


CONTRASTS = (
    Contrast("A", "S3 vs S2 (full depth)", "GSM8953516", "GSM8953515", "lognorm", "counts"),
    Contrast(
        "B", "S3 vs S1 (depth-matched)", "GSM8953516", "GSM8953514", "lognorm_dm", "depth_matched"
    ),
)


# Lineage/ambient-sensitive gene families (D-012): TCR segments, NK/cytotoxic and T-lineage genes,
# haemoglobins, mitochondrial and ribosomal transcripts.
EXCLUDE_PATTERN = re.compile(
    r"^(mt-|Rp[sl]\d|Hb[ab]-|Tr[abdg][vjc]|Tcr[bg]-|Klr|Gzm|Cd3[deg]$|Cd8[ab]|Ms4a4b$|Ms4a6b$|Ccl5$|Nkg7$)"
)


def excluded_genes(var_names, contamination: list[str]) -> set[str]:
    return set(contamination) | {g for g in var_names if EXCLUDE_PATTERN.match(g)}


def pseudobulk_lfc(counts, mask_sec, mask_pri, pseudocount: float = 1.0):
    """log2 ratio of library-level CPM (sum of counts over cells)."""
    import numpy as np

    s = np.asarray(counts[mask_sec].sum(axis=0)).ravel()
    p = np.asarray(counts[mask_pri].sum(axis=0)).ravel()
    cpm_s = s / max(s.sum(), 1) * 1e6
    cpm_p = p / max(p.sum(), 1) * 1e6
    return np.log2((cpm_s + pseudocount) / (cpm_p + pseudocount))


def run_contrast(ilc, c: Contrast, cell_mask=None):
    import numpy as np
    import pandas as pd
    import scanpy as sc

    lib = ilc.obs["sample"].astype(str).to_numpy()
    keep = np.isin(lib, [c.secondary, c.primary])
    if cell_mask is not None:
        keep &= cell_mask
    sub = ilc[keep].copy()
    sub.obs["group"] = np.where(
        sub.obs["sample"].astype(str) == c.secondary, "secondary", "primary"
    )
    counts = sub.layers[c.counts]
    sec = (sub.obs["group"] == "secondary").to_numpy()
    det_s = np.asarray((counts[sec] > 0).mean(axis=0)).ravel()
    det_p = np.asarray((counts[~sec] > 0).mean(axis=0)).ravel()
    tested = (det_s >= MIN_DETECT_TEST) | (det_p >= MIN_DETECT_TEST)
    sub = sub[:, tested].copy()
    sc.tl.rank_genes_groups(
        sub,
        "group",
        groups=["secondary"],
        reference="primary",
        method="wilcoxon",
        layer=c.layer,
        use_raw=False,
        tie_correct=True,
        pts=False,
    )
    r = sc.get.rank_genes_groups_df(sub, group="secondary").set_index("names")
    out = pd.DataFrame(index=sub.var_names)
    out["score"] = r["scores"]
    out["pval"] = r["pvals"]
    out["fdr"] = r["pvals_adj"]
    out["pb_log2fc"] = pseudobulk_lfc(sub.layers[c.counts], sec, ~sec)
    out["pct_secondary"] = det_s[tested]
    out["pct_primary"] = det_p[tested]
    out.attrs.update(n_secondary=int(sec.sum()), n_primary=int((~sec).sum()))
    return out


def concordant(a, b, excluded: set[str]):
    import numpy as np

    genes = a.index.intersection(b.index)
    a, b = a.loc[genes], b.loc[genes]
    same = np.sign(a["pb_log2fc"]) == np.sign(b["pb_log2fc"])
    sig = (a["fdr"] < FDR) & (b["fdr"] < FDR)
    big = (a["pb_log2fc"].abs() >= MIN_LFC) & (b["pb_log2fc"].abs() >= MIN_LFC)
    up = a["pb_log2fc"] > 0
    det = (
        np.where(
            up,
            np.minimum(a["pct_secondary"], b["pct_secondary"]),
            np.minimum(a["pct_primary"], b["pct_primary"]),
        )
        >= MIN_DETECT
    )
    clean = ~genes.isin(list(excluded))
    table = (
        a[["pb_log2fc", "fdr", "score"]]
        .add_suffix("_A")
        .join(b[["pb_log2fc", "fdr", "score"]].add_suffix("_B"))
    )
    table["direction"] = np.where(up, "up_in_secondary", "down_in_secondary")
    table["contamination_sensitive"] = ~clean
    table["candidate"] = same & sig & big & det & clean
    table["min_abs_lfc"] = np.minimum(a["pb_log2fc"].abs(), b["pb_log2fc"].abs())
    return table


def bootstrap_stability(ilc, genes, n_boot: int = N_BOOT, seed: int = 0):
    """Fraction of bootstrap resamples (cells within each library) in which both contrasts keep the
    observed sign and |log2FC| >= MIN_LFC, per gene."""
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)
    lib = ilc.obs["sample"].astype(str).to_numpy()
    gi = [ilc.var_names.get_loc(g) for g in genes]
    obs_sign = {}
    hits = np.zeros(len(gi))
    layers = {c.name: ilc.layers[c.counts][:, gi].tocsr() for c in CONTRASTS}
    for c in CONTRASTS:
        X = layers[c.name]
        obs_sign[c.name] = np.sign(pseudobulk_lfc(X, lib == c.secondary, lib == c.primary))
    idx = {s: np.flatnonzero(lib == s) for s in np.unique(lib)}
    for _ in range(n_boot):
        ok = np.ones(len(gi), dtype=bool)
        for c in CONTRASTS:
            X = layers[c.name]
            ss = rng.choice(idx[c.secondary], len(idx[c.secondary]))
            pp = rng.choice(idx[c.primary], len(idx[c.primary]))
            s = np.asarray(X[ss].sum(axis=0)).ravel()
            p = np.asarray(X[pp].sum(axis=0)).ravel()
            tot_s = ilc.layers[c.counts][ss].sum()
            tot_p = ilc.layers[c.counts][pp].sum()
            lfc = np.log2((s / tot_s * 1e6 + 1) / (p / tot_p * 1e6 + 1))
            ok &= (np.sign(lfc) == obs_sign[c.name]) & (np.abs(lfc) >= MIN_LFC)
        hits += ok
    return pd.Series(hits / n_boot, index=genes, name="bootstrap_stability")


def figures(tables: dict, conc, out_dir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    for ax, c in zip(axes[:2], CONTRASTS, strict=True):
        t = tables[c.name]
        y = -np.log10(t["fdr"].clip(lower=1e-300))
        cand = t.index.isin(conc.index[conc["candidate"]])
        ax.scatter(t["pb_log2fc"], y, s=3, c="#cbd5e1", rasterized=True)
        ax.scatter(
            t.loc[cand, "pb_log2fc"],
            y[cand],
            s=5,
            c="#dc2626",
            rasterized=True,
            label="concordant candidate",
        )
        for g in PUBLISHED_GENES:
            if g in t.index:
                ax.annotate(g, (t.loc[g, "pb_log2fc"], y[g]), fontsize=7)
        ax.axvline(MIN_LFC, ls=":", c="k", lw=0.7)
        ax.axvline(-MIN_LFC, ls=":", c="k", lw=0.7)
        ax.set_xlabel("pseudobulk log2FC (secondary / primary)")
        ax.set_ylabel("-log10 FDR (cell-level; ranking only)")
        ax.set_title(f"Contrast {c.name}: {c.label}")
        ax.legend(fontsize=8, frameon=False)
    ax = axes[2]
    ax.scatter(conc["pb_log2fc_A"], conc["pb_log2fc_B"], s=3, c="#cbd5e1", rasterized=True)
    m = conc["candidate"]
    ax.scatter(
        conc.loc[m, "pb_log2fc_A"], conc.loc[m, "pb_log2fc_B"], s=5, c="#dc2626", rasterized=True
    )
    r = np.corrcoef(conc["pb_log2fc_A"], conc["pb_log2fc_B"])[0, 1]
    ax.axhline(0, c="k", lw=0.5)
    ax.axvline(0, c="k", lw=0.5)
    ax.set_xlabel("log2FC  A: S3 vs S2")
    ax.set_ylabel("log2FC  B: S3 vs S1 (depth-matched)")
    ax.set_title(f"Agreement between contrasts (r = {r:.2f})")
    fig.tight_layout()
    paths.append(out_dir / "de_overview.png")
    fig.savefig(paths[-1], dpi=150)
    plt.close(fig)
    return paths


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


def render(tables, conc, stab, strat, published, figs, r_ab) -> str:
    cand = conc[conc["candidate"]].copy()
    cand["bootstrap_stability"] = stab
    n_up = int((cand["direction"] == "up_in_secondary").sum())
    n_dn = len(cand) - n_up
    med_stab = float(stab.median()) if len(stab) else float("nan")
    top_up = (
        cand[cand["direction"] == "up_in_secondary"]
        .sort_values("min_abs_lfc", ascending=False)
        .head(30)
    )
    top_dn = (
        cand[cand["direction"] == "down_in_secondary"]
        .sort_values("min_abs_lfc", ascending=False)
        .head(20)
    )
    cols = ["pb_log2fc_A", "pb_log2fc_B", "fdr_A", "fdr_B", "bootstrap_stability"]
    sizes = "; ".join(
        f"{c.name}: {tables[c.name].attrs['n_secondary']:,} vs {tables[c.name].attrs['n_primary']:,} "
        "ILC2s"
        for c in CONTRASTS
    )
    return "\n".join(
        [
            f"# Phase 6 — exploratory differential expression in ILC2s · {SERIES}",
            "",
            "> **Exploratory.** Each condition is represented by library-level data with an unknown mouse "
            "structure (D-004). Cell-level p-values are pseudoreplicated and used **only to rank genes**. "
            "Evidence comes from agreement between two contrasts that share no primary cells, bootstrap "
            "stability, and (Phase 10) independent datasets. Claim level: **A — recall-associated state**.",
            "",
            f"- Contrasts — {sizes}",
            f"- Candidate rule: FDR < {FDR} in both, |pseudobulk log2FC| ≥ {MIN_LFC} in both, same sign, detected in "
            f"≥ {MIN_DETECT:.0%} of the higher group, not contamination-sensitive.",
            f"- **Concordant candidates: {len(cand):,}** ({n_up} up, {n_dn} down in secondary).",
            f"- Genome-wide agreement of effect sizes between A and B: r = **{r_ab:.2f}**.",
            f"- Bootstrap sign stability of candidates ({N_BOOT} resamples): median **{med_stab:.2f}**.",
            f"- Non-cycling stratum (proliferation score ≤ {PROLIF_MAX}): **{strat:.0%}** of candidates remain "
            "concordant (same sign, FDR < 0.05 in both) — the signal is not explained by proliferating cells alone.",
            "",
            "## Genes reported in the source study",
            _md(published),
            "",
            "## Top concordant genes — up in secondary",
            _md(top_up[cols]),
            "",
            "## Top concordant genes — down in secondary",
            _md(top_dn[cols]),
            "",
            "Full tables: `results/differential_expression/`.",
            "",
            "## Figures",
            *[f"![{p.stem}](figures/{p.name})" for p in figs],
            "",
        ]
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--input", type=Path, default=Path("results/annotation") / f"{SERIES}_annotated.h5ad"
    )
    p.add_argument("--out-dir", type=Path, default=Path("results/differential_expression"))
    p.add_argument("--report", type=Path, default=Path("reports/phase6/06_de_report.md"))
    p.add_argument("--params", type=Path, default=None)
    p.add_argument("--markers", type=Path, default=None)
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    import numpy as np
    import pandas as pd
    import scanpy as sc

    from .annotate import load_markers
    from .config import parameters

    if not args.input.is_file():
        print(f"Error: {args.input} not found — run Phase 5 first.", file=sys.stderr)
        return 2
    seed = int(parameters(args.params).get("random_seed", 0))
    adata = sc.read_h5ad(args.input)
    ilc = adata[adata.obs["lineage"].astype(str) == "ILC2"].copy()
    excl = excluded_genes(
        ilc.var_names, load_markers(args.markers).get("contamination_sensitive", [])
    )
    print(f"[de] {ilc.n_obs:,} ILC2s")

    tables = {}
    for c in CONTRASTS:
        tables[c.name] = run_contrast(ilc, c)
        print(f"[de] contrast {c.name}: {len(tables[c.name]):,} genes tested")
    conc = concordant(tables["A"], tables["B"], excl)
    cand_genes = list(conc.index[conc["candidate"]])
    stab = bootstrap_stability(ilc, cand_genes, seed=seed) if cand_genes else pd.Series(dtype=float)

    noncyc = (ilc.obs["state_proliferation"] <= PROLIF_MAX).to_numpy()
    sa, sb = (run_contrast(ilc, c, noncyc) for c in CONTRASTS)
    sc_ = concordant(sa, sb, excl)
    kept = sc_.reindex(cand_genes)
    same_dir = (np.sign(kept["pb_log2fc_A"]) == np.sign(conc.loc[cand_genes, "pb_log2fc_A"])) & (
        np.sign(kept["pb_log2fc_B"]) == np.sign(conc.loc[cand_genes, "pb_log2fc_B"])
    )
    retained = (
        float(((kept["fdr_A"] < FDR) & (kept["fdr_B"] < FDR) & same_dir).mean())
        if cand_genes
        else 0.0
    )

    pub = [g for g in PUBLISHED_GENES if g in conc.index]
    published = conc.loc[
        pub,
        ["pb_log2fc_A", "fdr_A", "pb_log2fc_B", "fdr_B", "candidate", "contamination_sensitive"],
    ]
    published.index.name = "gene"
    r_ab = float(np.corrcoef(conc["pb_log2fc_A"], conc["pb_log2fc_B"])[0, 1])

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for c in CONTRASTS:
        tables[c.name].to_csv(args.out_dir / f"contrast_{c.name}.csv")
        tables[c.name]["score"].sort_values(ascending=False).to_csv(
            args.out_dir / f"ranked_{c.name}.rnk", sep="\t", header=False
        )
    out = conc.join(stab)
    out.to_csv(args.out_dir / "differential_expression_results.csv")
    sc_.to_csv(args.out_dir / "noncycling_stratum.csv")
    figs = figures(tables, conc, args.out_dir / "figures")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    (args.report.parent / "figures").mkdir(exist_ok=True)
    for f in figs:
        shutil.copy2(f, args.report.parent / "figures" / f.name)
    args.report.write_text(
        render(tables, conc, stab, retained, published, figs, r_ab), encoding="utf-8"
    )
    print(f"[de] {len(cand_genes)} concordant candidates · report → {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
