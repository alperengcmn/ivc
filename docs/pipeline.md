# Pipeline (v2) — phases, gates and status

Each phase ends in a **gate**: the next phase does not start until the gate condition holds.
Status legend: ✅ done · 🔄 in progress · ⏸ blocked · ⬜ not started.

| # | Phase | Key outputs | Gate to next phase | Status |
|---|-------|-------------|--------------------|--------|
| 0 | Hypothesis & success criteria | [`00_hypothesis_and_success_criteria.md`](00_hypothesis_and_success_criteria.md) | thresholds agreed | 🔄 thresholds TBD |
| 1 | Dataset discovery & validation | `validated_metadata`, compatibility report, `selected_datasets` | metadata validated **and** ATAC availability checked | 🔄 GSE295590 structure ✅; replicate mapping ❓ |
| 2 | Data acquisition | raw files + `download_log.tsv` | all accessions downloaded & organised | ⬜ `make fetch` ready |
| 3 | Single-cell QC | `qc_summary.csv`, filtered object, QC figures, `qc_report.md` | QC meets Phase 0 thresholds | ⬜ |
| 4 | Normalisation & dimensionality reduction | PCA/UMAP, clustering table | clustering stable, batch effects evaluated | ⬜ |
| 5 | ILC2 identification & subsets | annotated cells, marker tables, `ilc2_subset_report.md` | ≥ 1 marker category confirms ILC2; ≥ N cells | ⬜ |
| 6 | Differential expression | DE table, volcano plots, ranked lists | ≥ K biological replicates, FDR < 0.05 | ⏸ blocked (D-004) |
| 7 | Pathways & signatures | enrichment, gene-set scores, candidate signature | signature stable on ≥ 1 compatible dataset | ⬜ |
| 8 | TF activity | TF activity matrix, network, prioritised TFs | candidates cross-checked against target genes | ⬜ heavy compute |
| 9 | Epigenomic integration (conditional) | peaks, motifs, region–gene links | status recorded: performed / not performed / not testable | ⬜ |
| 10 | Cross-dataset validation | comparison table, validation report | every signature has a validation category | ⬜ |
| 11 | Statistical robustness | sensitivity analyses | sensitivity analyses complete | ⬜ |
| 12 | Final report & reproducibility | report, environment, seeds, configs | — | ⬜ |

**Validation categories (Phase 10):** Replicated · Partially supported · Dataset-specific ·
Not supported · Not testable with available data.

**Phase 9 fallback:** if no compatible chromatin data exist, Phase 9 is not run but stays in the
pipeline, and the final report states "Epigenomic integration: not performed — no compatible
dataset identified."

## Script map

| Script | Phase | Purpose |
|--------|-------|---------|
| `scripts/00_fetch_gse295590.sh` | 2 | download + checksum + log GSE295590 processed matrices |
| `scripts/01_validate_sample_metadata.py` | 1 | structural validation of sample metadata |
| `scripts/02_validate_geo_sra_crosswalk.py` | 1 | GEO → BioSample → SRX → SRR crosswalk validation |
| `scripts/03_…` | 3 | QC — next |
