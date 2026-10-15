"""Script 05 — Phase 5: lineage annotation, ILC2 confirmation and nILC2/iILC2 subsets.

1. Every cell gets a score per marker category (config/markers.yaml), computed on the
   depth-matched log-normalised data so libraries are comparable.
2. Clusters are labelled by the lineage category with the highest cluster-level z-score;
   no population is defined by a single marker. Clusters with a strong second lineage are
   flagged as possible doublets/mixtures.
3. ILC2 cells are re-embedded and sub-clustered; sub-clusters are labelled nILC2 / iILC2
   following the source study's definition (IL-5 vs IL-13 + Ki-67 + Il2ra), or left
   'ILC2_unassigned' when neither program is clearly dominant.
4. The published cell numbers are compared with ours as a reproduction check.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from collections.abc import Sequence
from pathlib import Path

from .gse295590 import BY_GEO, SERIES

SUBSET_MARGIN = 0.5  # z-score margin needed to call nILC2 vs iILC2
MIN_LINEAGE_SCORE = 0.1  # best lineage score below this → 'unassigned'
MIXED_RATIO = 0.5  # second-best lineage ≥ 50 % of the best → cluster flagged as mixed
# Source study (Fig. 3): nILC2 3,579 · iILC2 1,412 · primary 189 · secondary 4,802 ILC2s.
PUBLISHED = {"nILC2": 3579, "iILC2": 1412, "primary": 189, "secondary": 4802}


def load_markers(path: Path | None = None) -> dict:
    from .config import CONFIG_DIR, load_yaml

    return load_yaml(path or CONFIG_DIR / "markers.yaml")


def score_categories(adata, categories: dict[str, list[str]], prefix: str, seed: int = 0) -> dict:
    """Add `<prefix><name>` score columns. Returns {category: missing genes}."""
    import scanpy as sc

    missing = {}
    for name, genes in categories.items():
        present = [g for g in genes if g in adata.var_names]
        missing[name] = [g for g in genes if g not in adata.var_names]
        if len(present) < 2:
            raise ValueError(f"category {name!r}: fewer than 2 marker genes present")
        sc.tl.score_genes(
            adata, present, score_name=f"{prefix}{name}", random_state=seed, use_raw=False
        )
    return missing


def zscore_by_group(obs, columns: list[str], group: str):
    means = obs.groupby(group, observed=True)[columns].mean()
    return (means - means.mean()) / means.std(ddof=0).replace(0, 1), means


def label_lineage(adata, lineages: list[str]):
    """Label each cluster with the lineage whose *raw* mean score is highest.

    score_genes values are already centred on matched control genes, so raw means are
    comparable across categories; z-scoring across clusters is not used here because a
    lineage that dominates many clusters would be penalised. A cluster whose best score is
    below MIN_LINEAGE_SCORE is 'unassigned'; one whose second-best lineage reaches
    MIXED_RATIO of the best is flagged as mixed (possible doublets / contamination).
    """
    import numpy as np

    cols = [f"lin_{x}" for x in lineages]
    means = adata.obs.groupby("cluster", observed=True)[cols].mean()
    vals = means.to_numpy()
    order = np.argsort(-vals, axis=1)
    best = vals[np.arange(len(vals)), order[:, 0]]
    second = vals[np.arange(len(vals)), order[:, 1]]
    top = [
        lineages[i] if b >= MIN_LINEAGE_SCORE else "unassigned"
        for i, b in zip(order[:, 0], best, strict=True)
    ]
    table = means.copy()
    table.columns = lineages
    table["label"] = top
    table["second"] = [lineages[i] for i in order[:, 1]]
    table["second_to_best"] = np.where(best > 0, second / best, np.nan).round(2)
    table["mixed_flag"] = (second >= MIN_LINEAGE_SCORE) & (second >= MIXED_RATIO * best)
    table["n_cells"] = adata.obs["cluster"].value_counts().reindex(table.index).to_numpy()
    adata.obs["lineage"] = (
        adata.obs["cluster"].map(dict(zip(table.index, top, strict=True))).astype("category")
    )
    adata.obs["cluster_mixed_flag"] = adata.obs["cluster"].map(table["mixed_flag"]).astype(bool)
    return table


def subcluster_ilc2(adata, reference: str, seed: int = 0):
    """Resolve nILC2 / iILC2 on a clean reference library, then transfer to the others.

    Sub-clusters are sought only in the `reference` library (deep, low contamination),
    because a joint embedding separates libraries rather than states (D-012). If no
    resolution meets the Phase 4 stability rule, subsets are declared *not resolvable*
    and every ILC2 keeps the label 'ILC2'; the iILC2 program is then carried as a
    continuous per-cell score (`sub_iILC2`) instead of a discrete label.
    """
    import numpy as np
    import scanpy as sc
    from sklearn.neighbors import NearestCentroid

    from .embed import MIN_ARI, select_resolution

    ilc = adata[adata.obs["lineage"] == "ILC2"].copy()
    ilc.X = ilc.layers["lognorm_dm"].copy()
    sc.pp.highly_variable_genes(ilc, n_top_genes=2000, batch_key="sample", flavor="seurat")
    sc.pp.pca(ilc, n_comps=30, mask_var="highly_variable", random_state=seed)
    sc.pp.neighbors(ilc, n_pcs=20, random_state=seed)
    sc.tl.umap(ilc, random_state=seed)

    ref = ilc[ilc.obs["sample"].astype(str) == reference].copy()
    sc.pp.highly_variable_genes(ref, n_top_genes=2000, flavor="seurat")
    sc.pp.pca(ref, n_comps=30, mask_var="highly_variable", random_state=seed)
    sc.pp.neighbors(ref, n_pcs=20, random_state=seed)
    res, stability = select_resolution(ref, key="ref_subcluster")
    resolvable = bool(
        stability.loc[res, "min_ari"] >= MIN_ARI and stability.loc[res, "n_clusters"] >= 2
    )
    info = {
        "reference": reference,
        "resolution": res,
        "stability": stability,
        "resolvable": resolvable,
        "table": None,
    }

    if resolvable:
        z, means = zscore_by_group(ref.obs, ["sub_iILC2"], "ref_subcluster")
        labels = np.where(z["sub_iILC2"] > SUBSET_MARGIN, "iILC2", "nILC2")
        table = means.rename(columns={"sub_iILC2": "iILC2_program"})
        table["z"] = z["sub_iILC2"].round(2)
        table["subset"] = labels
        info["table"] = table
        ref_labels = ref.obs["ref_subcluster"].map(dict(zip(table.index, labels, strict=True)))
        genes = ref.var_names[ref.var["highly_variable"]]
        clf = NearestCentroid().fit(ref[:, genes].X.toarray(), ref_labels.to_numpy())
        ilc.obs["ilc2_subset"] = clf.predict(ilc[:, genes].X.toarray())
        ilc.obs.loc[ref.obs_names, "ilc2_subset"] = ref_labels.to_numpy()
    else:
        ilc.obs["ilc2_subset"] = "ILC2"
    ilc.obs["ilc2_subset"] = ilc.obs["ilc2_subset"].astype("category")
    adata.obs["ilc2_subset"] = "non-ILC2"
    adata.obs.loc[ilc.obs_names, "ilc2_subset"] = ilc.obs["ilc2_subset"].astype(str)
    adata.obs["ilc2_subset"] = adata.obs["ilc2_subset"].astype("category")
    return ilc, info


def figures(adata, ilc, out_dir: Path, marker_genes: list[str]) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import scanpy as sc

    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    xy = adata.obsm["X_umap"]
    cmap = plt.get_cmap("tab10")
    for i, lab in enumerate(sorted(adata.obs["lineage"].astype(str).unique())):
        m = (adata.obs["lineage"].astype(str) == lab).to_numpy()
        axes[0].scatter(
            xy[m, 0],
            xy[m, 1],
            s=1.5,
            color=cmap(i % 10),
            label=f"{lab} ({m.sum():,})",
            rasterized=True,
        )
    axes[0].legend(markerscale=6, fontsize=7, frameon=False, loc="best")
    axes[0].set_title("Lineage (all cells)")
    ixy = ilc.obsm["X_umap"]
    if set(ilc.obs["ilc2_subset"].astype(str)) <= {"ILC2"}:
        sca = axes[1].scatter(
            ixy[:, 0], ixy[:, 1], s=1.5, c=ilc.obs["sub_iILC2"], cmap="magma", rasterized=True
        )
        fig.colorbar(sca, ax=axes[1], label="iILC2 program score")
        axes[1].set_title("iILC2 program (subsets not resolvable)")
    else:
        colors = {"nILC2": "#2563eb", "iILC2": "#dc2626"}
        for lab, col in colors.items():
            m = (ilc.obs["ilc2_subset"].astype(str) == lab).to_numpy()
            axes[1].scatter(
                ixy[m, 0],
                ixy[m, 1],
                s=1.5,
                color=col,
                label=f"{lab} ({m.sum():,})",
                rasterized=True,
            )
        axes[1].legend(markerscale=6, fontsize=8, frameon=False)
        axes[1].set_title("ILC2 subsets (reference-mapped)")
    group = {"GSM8953514": "#2563eb", "GSM8953515": "#60a5fa", "GSM8953516": "#d97706"}
    for s, col in group.items():
        m = (ilc.obs["sample"].astype(str) == s).to_numpy()
        tag = BY_GEO[s].group if s in BY_GEO else ""
        axes[2].scatter(
            ixy[m, 0],
            ixy[m, 1],
            s=1.5,
            color=col,
            label=f"{s} {tag} ({m.sum():,})",
            rasterized=True,
        )
    axes[2].legend(markerscale=6, fontsize=8, frameon=False)
    axes[2].set_title("ILC2s by library")
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    paths.append(out_dir / "umap_annotation.png")
    fig.savefig(paths[-1], dpi=150)
    plt.close(fig)

    genes = [g for g in marker_genes if g in adata.var_names]
    adata.obs["lineage_subset"] = np.where(
        adata.obs["lineage"].astype(str) == "ILC2",
        adata.obs["ilc2_subset"].astype(str),
        adata.obs["lineage"].astype(str),
    )
    dp = sc.pl.dotplot(
        adata,
        genes,
        groupby="lineage_subset",
        layer="lognorm_dm",
        standard_scale="var",
        show=False,
        return_fig=True,
    )
    dp.savefig(out_dir / "marker_dotplot.png", dpi=150, bbox_inches="tight")
    plt.close("all")
    paths.append(out_dir / "marker_dotplot.png")
    return paths


def _md(df, floatfmt="{:.2f}") -> str:
    cols = [df.index.name or "", *map(str, df.columns)]
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for idx, row in df.iterrows():
        cells = []
        for v in row:
            if isinstance(v, bool):
                cells.append("⚠" if v else "")
            elif isinstance(v, float):
                cells.append(floatfmt.format(v))
            elif isinstance(v, int):
                cells.append(f"{v:,}")
            else:
                cells.append(str(v))
        out.append("| " + " | ".join([str(idx), *cells]) + " |")
    return "\n".join(out)


def _stab(df) -> str:
    lines = ["| resolution | clusters | min ARI | mean ARI |", "|---|---|---|---|"]
    lines += [
        f"| {r} | {int(x.n_clusters)} | {x.min_ari:.3f} | {x.mean_ari:.3f} |"
        for r, x in df.iterrows()
    ]
    return "\n".join(lines)


def render(lineage_table, info, counts, repro, missing, contam, gate_ok, min_ilc2, figs) -> str:
    miss = "; ".join(f"{k}: {', '.join(v)}" for k, v in missing.items() if v) or "none"
    if info["resolvable"]:
        subset_text = [
            f"Stable sub-clusters found in the reference library `{info['reference']}` at resolution "
            f"{info['resolution']}; labelled by the iILC2 program (z > {SUBSET_MARGIN} → iILC2) and transferred "
            "to the other libraries by nearest centroid.",
            "",
            _md(info["table"]),
        ]
    else:
        subset_text = [
            f"**nILC2 / iILC2 are not resolvable as discrete, stable sub-clusters** in the reference library "
            f"`{info['reference']}`: no resolution reaches min ARI ≥ 0.8 across 5 Leiden seeds (table below). "
            "ILC2 states form a continuum here, so no discrete subset labels are assigned (D-012). The iILC2 "
            "program (Il13, Il2ra, Il1rl1, Mki67, Top2a — the source study's description) is carried as a "
            "continuous per-cell score `sub_iILC2`, and Phase 6 stratifies by proliferation instead.",
        ]
    return "\n".join(
        [
            f"# Phase 5 — ILC2 identification and subset annotation · {SERIES}",
            "",
            f"- **PHASE 5 GATE: {'PASS' if gate_ok else 'FAIL'}** — ILC2 confirmed by ILC2-identity *and* "
            f"type-2 effector categories; ≥ {min_ilc2} ILC2 cells required per library.",
            f"- Marker genes not present in the data: {miss}",
            "",
            "## Lineage per cluster (mean category scores; label = highest raw score)",
            _md(lineage_table),
            "",
            f"⚠ = second-best lineage score ≥ {MIXED_RATIO:.0%} of the best: mixed cluster (possible doublets or "
            f"ambient contamination). Best score < {MIN_LINEAGE_SCORE} → 'unassigned'.",
            "",
            "## ILC2 subsets",
            *subset_text,
            "",
            _stab(info["stability"]),
            "",
            "## Cells per library",
            _md(counts, "{:,.0f}"),
            "",
            "## Reproduction check against the published numbers",
            "The source study used mito ≤ 10 % (we use ≤ 20 %, D-010), so its rule is applied here for comparison.",
            "",
            _md(repro, "{:,.0f}"),
            "",
            "## Lineage-contamination score in ILC2s, by library",
            "Mean of the highest non-ILC2 lineage score per ILC2 cell (depth-matched data). Higher values mean more "
            "T/NK/erythroid transcripts inside ILC2-classified cells.",
            "",
            _md(contam, "{:.3f}"),
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
        "--input", type=Path, default=Path("results/embedding") / f"{SERIES}_embedded.h5ad"
    )
    p.add_argument(
        "--out", type=Path, default=Path("results/annotation") / f"{SERIES}_annotated.h5ad"
    )
    p.add_argument("--report", type=Path, default=Path("reports/phase5/05_annotation_report.md"))
    p.add_argument("--params", type=Path, default=None)
    p.add_argument("--markers", type=Path, default=None)
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    import pandas as pd
    import scanpy as sc

    from .config import parameters

    if not args.input.is_file():
        print(f"Error: {args.input} not found — run Phase 4 first.", file=sys.stderr)
        return 2
    params = parameters(args.params)
    seed = int(params.get("random_seed", 0))
    min_ilc2 = int((params.get("phase5_annotation") or {}).get("min_confirmed_ilc2_cells") or 0)
    markers = load_markers(args.markers)

    adata = sc.read_h5ad(args.input)
    adata.X = adata.layers["lognorm_dm"].copy()
    missing = {}
    missing |= score_categories(adata, markers["lineage"], "lin_", seed)
    missing |= score_categories(adata, markers["state"], "state_", seed)
    missing |= score_categories(adata, markers["subset"], "sub_", seed)

    lineages = list(markers["lineage"])
    lineage_table = label_lineage(adata, lineages)
    ilc_clusters = lineage_table.index[lineage_table["label"] == "ILC2"]
    effector = adata.obs.groupby("cluster", observed=True)["state_type2_effector"].mean()
    confirmed = bool(len(ilc_clusters)) and bool(
        (effector.loc[ilc_clusters] > effector.mean()).any()
    )

    reference = (params.get("phase5_annotation") or {}).get(
        "subset_reference_library", "GSM8953516"
    )
    ilc, info = subcluster_ilc2(adata, reference, seed)
    counts = pd.crosstab(adata.obs["ilc2_subset"].astype(str), adata.obs["sample"].astype(str))
    counts.index.name = "class"
    ilc_per_lib = counts.drop(index="non-ILC2", errors="ignore").sum()
    gate_ok = confirmed and bool((ilc_per_lib >= min_ilc2).all())

    is_ilc = adata.obs["lineage"].astype(str) == "ILC2"
    lib = adata.obs["sample"].astype(str)
    strict = adata.obs["pct_counts_mt"] <= 10
    rows = {
        "ILC2 · GSM8953515 (primary S2), mito ≤ 10 %": (
            int((is_ilc & strict & (lib == "GSM8953515")).sum()),
            PUBLISHED["primary"],
        ),
        "ILC2 · GSM8953514 + GSM8953515 (all primary), mito ≤ 10 %": (
            int((is_ilc & strict & lib.isin(["GSM8953514", "GSM8953515"])).sum()),
            PUBLISHED["primary"],
        ),
        "ILC2 · GSM8953516 (secondary), mito ≤ 10 %": (
            int((is_ilc & strict & (lib == "GSM8953516")).sum()),
            PUBLISHED["secondary"],
        ),
    }
    repro = pd.DataFrame(rows, index=["this_analysis", "published"]).T.astype(float)
    repro.index.name = "population"

    non_ilc = [f"lin_{x}" for x in lineages if x != "ILC2"]
    adata.obs["contamination_score"] = adata.obs[non_ilc].max(axis=1)
    contam = (
        adata.obs.loc[is_ilc]
        .groupby("sample", observed=True)["contamination_score"]
        .agg(["mean", "median"])
        .rename_axis("library")
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    keep = [
        "sample",
        "infection_group",
        "cluster",
        "lineage",
        "cluster_mixed_flag",
        "ilc2_subset",
        "contamination_score",
        "doublet_score",
        *[c for c in adata.obs if c.startswith(("lin_", "state_", "sub_"))],
    ]
    adata.obs[keep].to_csv(args.out.parent / "annotated_cells.csv.gz")
    lineage_table.to_csv(args.out.parent / "lineage_cluster_table.csv")
    info["stability"].to_csv(args.out.parent / "ilc2_subset_stability.csv")
    if info["table"] is not None:
        info["table"].to_csv(args.out.parent / "ilc2_subcluster_table.csv")
    adata.write_h5ad(args.out, compression="gzip")
    ilc.write_h5ad(args.out.parent / f"{SERIES}_ilc2.h5ad", compression="gzip")

    marker_genes = [
        "Gata3",
        "Il1rl1",
        "Il7r",
        "Thy1",
        "Klrg1",
        "Arg1",
        "Il5",
        "Calca",
        "Il13",
        "Il2ra",
        "Mki67",
        "Areg",
        "Cd3e",
        "Trac",
        "Rorc",
        "Il17a",
        "Nkg7",
        "Ncr1",
        "Cd19",
        "Cd79a",
        "Lyz2",
        "Cpa3",
        "Hbb-bs",
        "Epcam",
    ]
    figs = figures(adata, ilc, args.out.parent / "figures", marker_genes)
    (args.report.parent / "figures").mkdir(parents=True, exist_ok=True)
    for f in figs:
        shutil.copy2(f, args.report.parent / "figures" / f.name)
    args.report.write_text(
        render(lineage_table, info, counts, repro, missing, contam, gate_ok, min_ilc2, figs),
        encoding="utf-8",
    )
    print(f"[annotate] gate {'PASS' if gate_ok else 'FAIL'} · report → {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
