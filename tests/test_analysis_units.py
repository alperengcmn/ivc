"""Unit tests for the decision rules of Phases 4–11 (fast, synthetic inputs only)."""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("scanpy")


def test_thinning_matches_shallowest_median():
    import scipy.sparse as sp

    from ilc2tri.embed import thin_to_shallowest

    rng = np.random.default_rng(0)
    deep = rng.poisson(8, (300, 200))
    shallow = rng.poisson(1, (300, 200))
    X = sp.csr_matrix(np.vstack([deep, shallow]))
    samples = ["deep"] * 300 + ["shallow"] * 300
    Y, medians, factors = thin_to_shallowest(X, samples, seed=0)
    tot = np.asarray(Y.sum(axis=1)).ravel()
    assert factors["shallow"] == 1.0 and factors["deep"] < 0.2
    assert abs(np.median(tot[:300]) - medians["shallow"]) / medians["shallow"] < 0.05
    assert (Y.toarray() <= X.toarray()).all()  # thinning never adds counts


def test_concordance_rule():
    from ilc2tri.de import concordant

    idx = ["up_ok", "flip", "small", "ns", "contam"]
    a = pd.DataFrame(
        {
            "pb_log2fc": [1.0, 1.0, 0.2, 1.0, 2.0],
            "fdr": [1e-5, 1e-5, 1e-5, 0.2, 1e-9],
            "score": 1.0,
            "pct_secondary": 0.5,
            "pct_primary": 0.5,
        },
        index=idx,
    )
    b = a.copy()
    b.loc["flip", "pb_log2fc"] = -1.0
    t = concordant(a, b, excluded={"contam"})
    assert t.loc["up_ok", "candidate"]
    assert not t.loc[["flip", "small", "ns", "contam"], "candidate"].any()


def test_contamination_patterns():
    from ilc2tri.de import excluded_genes

    ex = excluded_genes(
        [
            "Trdv4",
            "Tcrg-C1",
            "Klrb1c",
            "Gzmb",
            "Cd3e",
            "Cd38",
            "Hbb-bs",
            "Hbegf",
            "Rpl13",
            "mt-Co1",
            "Il13",
        ],
        [],
    )
    assert {"Trdv4", "Tcrg-C1", "Klrb1c", "Gzmb", "Cd3e", "Hbb-bs", "Rpl13", "mt-Co1"} <= ex
    assert not {"Cd38", "Hbegf", "Il13"} & ex


@pytest.mark.parametrize(
    ("score_p", "score_opp", "gene_p", "gene_opp", "n", "expected"),
    [
        (0.01, 0.99, 0.01, 0.99, 80, "Replicated"),
        (0.01, 0.99, 0.30, 0.70, 80, "Partially supported"),
        (None, None, 0.01, 0.99, 80, "Partially supported"),
        (0.40, 0.60, 0.30, 0.70, 80, "Not supported"),
        (0.99, 0.01, 0.01, 0.99, 80, "Not supported"),  # significant in the opposite direction
        (0.01, 0.99, 0.01, 0.99, 10, "Not testable with available data"),
    ],
)
def test_preregistered_verdicts(score_p, score_opp, gene_p, gene_opp, n, expected):
    from ilc2tri.validation import verdict

    assert (
        verdict(
            {"p": score_p, "p_opposite": score_opp},
            {"n_genes": n, "p": gene_p, "p_opposite": gene_opp},
        )
        == expected
    )


def test_gene_agreement_counts_directions():
    from ilc2tri.validation import gene_agreement

    lfc = pd.Series({"u1": 1.0, "u2": -0.5, "d1": -1.0, "d2": -2.0})
    g = gene_agreement(lfc, ["u1", "u2"], ["d1", "d2", "missing"])
    assert g["n_genes"] == 4 and g["agree"] == 3
    assert g["up_agree"] == "1/2" and g["down_agree"] == "2/2"


def test_nearest_tss_within_window():
    from ilc2tri.chromatin import nearest_gene

    tss = {"chr1": ([1_000, 200_000], ["GeneA", "GeneB"])}
    assert nearest_gene("chr1:1100-1300", tss) == ("GeneA", 200)
    assert nearest_gene("chr1:120000-120200", tss) == (None, None)  # > 50 kb from both
    assert nearest_gene("chr2:1-100", tss) == (None, None)


def test_signature_file_matches_manifest():
    import hashlib
    from pathlib import Path

    root = Path(__file__).resolve().parents[1] / "signatures"
    manifest = pd.read_csv(root / "MANIFEST.tsv", sep="\t")
    for row in manifest.itertuples():
        assert hashlib.sha256((root / row.file).read_bytes()).hexdigest() == row.sha256
