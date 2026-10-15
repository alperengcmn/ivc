# Pipeline (v2) — phases, gates and status

Each phase ends in a **gate**: the next phase does not start until the gate condition holds.
Status legend: ✅ done · 🔄 in progress · ⏸ blocked · ⬜ not started.

| # | Phase | Key outputs | Gate to next phase | Status |
|---|-------|-------------|--------------------|--------|
| 0 | Hypothesis & success criteria | [`00_hypothesis_and_success_criteria.md`](00_hypothesis_and_success_criteria.md) | thresholds agreed | 🟡 provisional thresholds (D-010) — Divya to review |
| 1 | Dataset discovery & validation | `validated_metadata`, compatibility report, `selected_datasets` | metadata validated **and** ATAC availability checked | ✅ structure validated; mouse mapping open (D-004) |
| 2 | Data acquisition | raw files + `download_log.tsv` | all accessions downloaded & organised | ✅ all 8 series + resources, SHA-256 logged |
| 3 | Single-cell QC | `qc_summary.csv`, filtered object, QC figures, `qc_report.md` | QC meets Phase 0 thresholds | ✅ PASS with provisional thresholds |
| 4 | Normalisation & dimensionality reduction | PCA/UMAP, clustering table | clustering stable, batch effects evaluated | ✅ PASS (depth-matched, D-011) |
| 5 | ILC2 identification & subsets | annotated cells, marker tables, `ilc2_subset_report.md` | ≥ 1 marker category confirms ILC2; ≥ N cells | ✅ PASS; n/iILC2 not resolvable (D-012) |
| 6 | Differential expression | DE table, volcano plots, ranked lists | ≥ K biological replicates, FDR < 0.05 | ✅ exploratory library-consistency design (D-013) |
| 7 | Pathways & signatures | enrichment, gene-set scores, candidate signature | signature stable on ≥ 1 compatible dataset | ✅ signature v1 frozen; external stability → Phase 10 |
| 8 | TF activity | TF activity matrix, network, prioritised TFs | candidates cross-checked against target genes | ✅ PASS (decoupler ULM) |
| 9 | Epigenomic integration (conditional) | peaks, motifs, region–gene links | status recorded: performed / not performed / not testable | ✅ performed — partial support; motifs not performed |
| 10 | Cross-dataset validation | comparison table, validation report | every signature has a validation category | ✅ pre-registered verdicts (D-014) |
| 11 | Statistical robustness | sensitivity analyses | sensitivity analyses complete | ✅ PASS |
| 12 | Final report & reproducibility | report, environment, seeds, configs | — | ✅ `reports/FINAL_REPORT.md` (generated) |

**Validation categories (Phase 10):** Replicated · Partially supported · Dataset-specific ·
Not supported · Not testable with available data.

**Phase 9 fallback:** if no compatible chromatin data exist, Phase 9 is not run but stays in the
pipeline, and the final report states "Epigenomic integration: not performed — no compatible
dataset identified."

## Script map

| Script | Phase | Purpose |
|--------|-------|---------|
| `scripts/00_fetch_gse295590.sh`, `fetch_geo.py`, `fetch_resources.py` | 2 | downloads with checksums |
| `scripts/01_validate_sample_metadata.py` | 1 | structural validation of sample metadata |
| `scripts/02_validate_geo_sra_crosswalk.py` | 1 | GEO → BioSample → SRX → SRR crosswalk |
| `scripts/03_qc.py` | 3 | QC metrics, doublets, flags, filtering |
| `scripts/04_embed.py` | 4 | depth-matched embedding, stability-selected clustering |
| `scripts/05_annotate.py` | 5 | lineage annotation, ILC2 subsets, reproduction check |
| `scripts/06_de.py` | 6 | exploratory library-consistency DE |
| `scripts/07_pathways.py` | 7 | GSEA, frozen signature |
| `scripts/08_tf_activity.py` | 8 | CollecTRI TF and PROGENy pathway activity |
| `scripts/09_chromatin.py` | 9 | GSE172258 ATAC vs signature |
| `scripts/10_validate.py` | 10 | pre-registered cross-dataset validation |
| `scripts/11_robustness.py` | 11 | sensitivity analyses |
| `scripts/12_final_report.py` | 12 | final report from results |
