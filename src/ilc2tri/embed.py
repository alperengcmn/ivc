"""Script 04 — Phase 4: depth-aware normalisation, dimensionality reduction and clustering.

Why depth-matched: GSM8953514 is sequenced ~5.6× shallower than the other libraries
(D-009). Embedding raw-depth data makes PC1 a depth axis. We therefore build the
embedding from *binomially thinned* counts: every library is thinned to the median
depth of the shallowest one (seeded, reproducible). Full-depth counts are kept in
`layers["counts"]` for differential expression.

No cross-library batch correction is applied: the only secondary library is also a
separate batch, so any correction across libraries would remove the biology under study
(D-011). Batch effects are *evaluated* and reported instead.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from .gse295590 import SERIES

RESOLUTION_GRID = (0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2)
N_HVG = 2000
N_PCS = 30
SEEDS = (0, 1, 2, 3, 4)  # clustering stability check
MIN_ARI = 0.8  # pre-specified: highest resolution whose min ARI across seeds >= MIN_ARI


def thin_to_shallowest(counts, samples, seed: int = 0):
    """Binomially thin each library so its median depth equals the shallowest library's."""
    import numpy as np
    import scipy.sparse as sp

    X = sp.csr_matrix(counts, dtype=np.float32, copy=True)
    totals = np.asarray(X.sum(axis=1)).ravel()
    samples = np.asarray(samples).astype(str)
    medians = {s: float(np.median(totals[samples == s])) for s in np.unique(samples)}
    target = min(medians.values())
    rng = np.random.default_rng(seed)
    factors = {}
    for s, med in sorted(medians.items()):
        p = min(1.0, target / med) if med > 0 else 1.0
        factors[s] = p
        if p >= 1.0:
            continue
        for r in np.flatnonzero(samples == s):
            lo, hi = X.indptr[r], X.indptr[r + 1]
            X.data[lo:hi] = rng.binomial(X.data[lo:hi].astype(np.int64), p)
    X.eliminate_zeros()
    return X, medians, factors


def lognorm(X):
    import numpy as np
    import scipy.sparse as sp

    X = sp.csr_matrix(X, dtype=np.float32, copy=True)
    totals = np.asarray(X.sum(axis=1)).ravel()
    totals[totals == 0] = 1
    X = sp.diags(1e4 / totals) @ X
    X.data = np.log1p(X.data)
    return sp.csr_matrix(X, dtype=np.float32)


def select_resolution(adata, grid=RESOLUTION_GRID, seeds=SEEDS, min_ari=MIN_ARI, key="cluster"):
    """Pre-specified rule: highest resolution whose min ARI (seed 0 vs others) >= min_ari.

    Writes the chosen clustering to adata.obs[key]; returns (resolution, stability table).
    """
    import pandas as pd
    import scanpy as sc
    from sklearn.metrics import adjusted_rand_score

    rows, labels = [], {}
    for res in grid:
        runs = []
        for s in seeds:
            sc.tl.leiden(
                adata,
                resolution=res,
                random_state=s,
                flavor="igraph",
                n_iterations=2,
                directed=False,
                key_added="_tmp",
            )
            runs.append(adata.obs["_tmp"].copy())
        aris = [adjusted_rand_score(runs[0], r) for r in runs[1:]]
        labels[res] = runs[0]
        rows.append(
            {
                "resolution": res,
                "n_clusters": runs[0].nunique(),
                "min_ari": min(aris),
                "mean_ari": sum(aris) / len(aris),
            }
        )
    del adata.obs["_tmp"]
    table = pd.DataFrame(rows).set_index("resolution").round(3)
    ok = table.index[table["min_ari"] >= min_ari]
    chosen = float(max(ok)) if len(ok) else float(min(grid))
    adata.obs[key] = labels[chosen].astype("category")
    return chosen, table


def embed(adata, seed: int = 0) -> dict:
    """Populate layers, PCA, neighbours, UMAP and Leiden clusterings. Returns diagnostics."""
    import numpy as np
    import pandas as pd
    import scanpy as sc

    samples = adata.obs["sample"].astype(str).to_numpy()
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    thinned, medians, factors = thin_to_shallowest(adata.layers["counts"], samples, seed)
    adata.layers["depth_matched"] = thinned
    adata.layers["lognorm"] = lognorm(adata.layers["counts"])
    adata.layers["lognorm_dm"] = lognorm(thinned)
    adata.X = adata.layers["lognorm_dm"].copy()

    diag: dict = {"medians": medians, "thinning": factors}
    # PC1–depth correlation before thinning (diagnostic only)
    full = adata.copy()
    full.X = adata.layers["lognorm"].copy()
    sc.pp.highly_variable_genes(full, n_top_genes=N_HVG, batch_key="sample", flavor="seurat")
    sc.pp.pca(full, n_comps=N_PCS, mask_var="highly_variable", random_state=seed)
    depth_full = np.log1p(np.asarray(adata.layers["counts"].sum(axis=1)).ravel())
    diag["pc1_depth_r_full"] = float(abs(np.corrcoef(full.obsm["X_pca"][:, 0], depth_full)[0, 1]))
    del full

    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, batch_key="sample", flavor="seurat")
    sc.pp.pca(adata, n_comps=50, mask_var="highly_variable", random_state=seed)
    depth_dm = np.log1p(np.asarray(thinned.sum(axis=1)).ravel())
    diag["pc1_depth_r_matched"] = float(abs(np.corrcoef(adata.obsm["X_pca"][:, 0], depth_dm)[0, 1]))
    sc.pp.neighbors(adata, n_pcs=N_PCS, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    diag["resolution"], diag["stability"] = select_resolution(adata)
    diag["n_clusters"] = int(adata.obs["cluster"].nunique())
    diag["composition"] = pd.crosstab(adata.obs["cluster"], adata.obs["sample"])
    return diag


def figures(adata, out_dir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    out_dir.mkdir(parents=True, exist_ok=True)
    xy = adata.obsm["X_umap"]
    paths = []
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    samples = sorted(adata.obs["sample"].astype(str).unique())
    palette = {"GSM8953514": "#2563eb", "GSM8953515": "#60a5fa", "GSM8953516": "#d97706"}
    for s in samples:
        m = (adata.obs["sample"].astype(str) == s).to_numpy()
        axes[0].scatter(
            xy[m, 0], xy[m, 1], s=1.5, c=palette.get(s, "#888"), label=s, rasterized=True
        )
    axes[0].legend(markerscale=6, fontsize=8, frameon=False)
    axes[0].set_title("Library (depth-matched embedding)")
    cl = adata.obs["cluster"].astype(str).to_numpy()
    cmap = plt.get_cmap("tab20")
    for i, c in enumerate(sorted(set(cl), key=int)):
        m = cl == c
        axes[1].scatter(xy[m, 0], xy[m, 1], s=1.5, color=cmap(i % 20), rasterized=True)
        axes[1].text(*np.median(xy[m], axis=0), c, fontsize=9, weight="bold", ha="center")
    axes[1].set_title("Leiden clusters")
    sca = axes[2].scatter(
        xy[:, 0],
        xy[:, 1],
        s=1.5,
        c=np.log10(adata.obs["total_counts"]),
        cmap="viridis",
        rasterized=True,
    )
    fig.colorbar(sca, ax=axes[2], label="log10 UMI (full depth)")
    axes[2].set_title("Original sequencing depth")
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    paths.append(out_dir / "umap_overview.png")
    fig.savefig(paths[-1], dpi=150)
    plt.close(fig)
    return paths


def _md_table(df) -> str:
    cols = [df.index.name or "", *map(str, df.columns)]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += [
        "| " + " | ".join([str(i), *[f"{v:,}" for v in row]]) + " |" for i, row in df.iterrows()
    ]
    return "\n".join(lines)


def _stab_table(df) -> str:
    lines = ["| resolution | clusters | min ARI | mean ARI |", "|---|---|---|---|"]
    lines += [
        f"| {r} | {int(row.n_clusters)} | {row.min_ari:.3f} | {row.mean_ari:.3f} |"
        for r, row in df.iterrows()
    ]
    return "\n".join(lines)


def render(diag: dict, n_cells: int, figs: list[Path]) -> str:
    stable = float(diag["stability"].loc[diag["resolution"], "min_ari"]) >= MIN_ARI
    reduced = diag["pc1_depth_r_matched"] < 0.3
    gate = "PASS" if stable and reduced else "REVIEW"
    th = ", ".join(f"`{s}` × {p:.3f}" for s, p in diag["thinning"].items())
    med = ", ".join(f"`{s}` {m:,.0f}" for s, m in diag["medians"].items())
    return "\n".join(
        [
            f"# Phase 4 — normalisation, embedding and clustering · {SERIES}",
            "",
            f"- Cells: **{n_cells:,}** (after provisional Phase 3 filtering, D-010)",
            f"- **PHASE 4 GATE: {gate}** — clustering stable across {len(SEEDS)} seeds (min ARI ≥ {MIN_ARI}: {stable}); "
            f"depth no longer drives PC1 (|r| < 0.3: {reduced}).",
            "",
            "## Depth handling",
            f"- Median UMI per cell (full depth): {med}",
            f"- Binomial thinning factors to the shallowest library: {th}",
            f"- |r| between PC1 and log depth: **{diag['pc1_depth_r_full']:.2f}** on full-depth data → "
            f"**{diag['pc1_depth_r_matched']:.2f}** on depth-matched data.",
            "- Full-depth counts are retained (`layers['counts']`) for differential expression.",
            "",
            "## Batch evaluation (no correction applied — D-011)",
            "The secondary condition has a single library, so library and condition are fully confounded; "
            "any cross-library correction would remove the biology under study. Composition per cluster:",
            "",
            _md_table(diag["composition"]),
            "",
            "## Resolution selection (pre-specified rule)",
            f"Highest resolution whose minimum ARI across {len(SEEDS)} Leiden seeds is ≥ {MIN_ARI} → "
            f"**{diag['resolution']}** ({diag['n_clusters']} clusters).",
            "",
            _stab_table(diag["stability"]),
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
    p.add_argument("--input", type=Path, default=Path("results/qc") / f"{SERIES}_qc_filtered.h5ad")
    p.add_argument(
        "--out", type=Path, default=Path("results/embedding") / f"{SERIES}_embedded.h5ad"
    )
    p.add_argument("--report", type=Path, default=Path("reports/phase4/04_embedding_report.md"))
    p.add_argument("--params", type=Path, default=None)
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    import shutil

    import scanpy as sc

    from .config import parameters

    if not args.input.is_file():
        print(
            f"Error: {args.input} not found — run Phase 3 with agreed thresholds first.",
            file=sys.stderr,
        )
        return 2
    seed = int(parameters(args.params).get("random_seed", 0))
    adata = sc.read_h5ad(args.input)
    print(f"[embed] {adata.n_obs:,} cells · seed {seed}")
    diag = embed(adata, seed=seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.obs[["sample", "cluster"]].to_csv(args.out.parent / "clustering_table.csv")
    diag["stability"].to_csv(args.out.parent / "resolution_stability.csv")
    adata.write_h5ad(args.out, compression="gzip")
    figs = figures(adata, args.out.parent / "figures")
    (args.report.parent / "figures").mkdir(parents=True, exist_ok=True)
    for f in figs:
        shutil.copy2(f, args.report.parent / "figures" / f.name)
    args.report.write_text(render(diag, adata.n_obs, figs), encoding="utf-8")
    print(f"[embed] report → {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
