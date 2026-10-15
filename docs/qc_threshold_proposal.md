# Phase 3 QC thresholds — proposal for GSE295590

**Status:** PROPOSED by Alperen · awaiting Divya's review. `config/parameters.yaml` keeps `null`
until both sign off (CONTRIBUTING: thresholds are agreed, not tuned). Once agreed, copy the
block below into the config, log it as **D-010**, and re-run `make qc`.

Evidence: [`reports/qc/03_qc_report.md`](../reports/qc/03_qc_report.md) (report-only run on the
deposited matrices, 2026-10-15) and `results/qc/qc_figures/`.

## What the data show

| Observation | Evidence | Consequence |
|---|---|---|
| Matrices are already filtered at ≈ 500 UMI/cell | minimum UMI 505 / 533 / 500 in S1 / S2 / S3 | no UMI floor needed on our side; sub-500 barcodes are unrecoverable |
| **S1 is sequenced ~5.6× shallower** than S2/S3 | median UMI 1,283 vs 7,160 / 7,215; median genes 732 vs ~2,100–2,400 | **No single gene-count cut-off can be "fair"** to all libraries; depth is a confounder for Phase 4 onwards |
| S3 contains a distinct damaged-cell population | 10.1 % of S3 cells > 20 % mito; 9.0 % have < 500 genes **and** > 20 % mito; visible as a separate cloud in `qc_counts_vs_genes.png` | a mito ceiling removes it cleanly |
| S2 has a small version of the same tail | 1.4 % of cells > 20 % mito | same rule applies |
| Scrublet's automatic threshold failed | 0 % predicted doublets in S1 (9,012 cells) and S3 (8,677 cells) | do **not** filter on `predicted_doublet`; assess doublets per cluster in Phase 5 |
| MAD-based mito ceilings are very tight | suggested 3.9–6.6 % → would remove 8–15 % of cells | too stringent here; see biological caveat below |

### Biological caveat on the mito ceiling

The source study reports **increased oxidative phosphorylation** in trained (secondary) ILC2s.
Mitochondrial transcript fraction partly tracks that biology, so an aggressive ceiling would
preferentially remove secondary-response cells and bias the very comparison this project is about.
The ceiling should therefore remove the clearly damaged cloud (> 20 % mito, low complexity) and no
more; a stricter ceiling is tested as a sensitivity analysis, not used as the default.

## Proposed values

```yaml
phase3_qc:
  min_cells_per_sample: 500      # library-level floor; S2 (811 after filtering) passes, flagged as low recovery
  max_pct_mito: 20               # removes the damaged high-mito cloud; spares OXPHOS-high cells
  min_genes_per_cell: 200        # standard floor; S1 minimum is 247, so S1 is untouched by it
  max_genes_per_cell: 8000       # safeguard only (~S3 99.9th percentile); doublets handled in Phase 5
  doublet_method: scrublet
  remove_predicted_doublets: false
```

## Effect (dry run on the real data)

| Sample | Before | After | Removed |
|---|---|---|---|
| S1 · GSM8953514 | 9,012 | 9,011 | 1 (0.0 %) |
| S2 · GSM8953515 | 838 | 811 | 27 (3.2 %) |
| S3 · GSM8953516 | 8,677 | 7,753 | 924 (10.6 %) |
| **Total** | **18,527** | **17,575** | **952 (5.1 %)** |

Phase 3 gate under these values: **PASS** (every sample ≥ 500 cells after filtering).

## Sensitivity analyses to carry into Phase 11

| Variant | S1 | S2 | S3 |
|---|---|---|---|
| Proposed (mito ≤ 20 %) | 9,011 | 811 | 7,753 |
| Strict mito ≤ 10 % | 8,995 | 782 | 7,545 |
| Strict genes ≥ 500 | 8,514 | 772 | 7,671 |

Findings that change direction between these variants are reported as threshold-sensitive.

## Questions for Divya

1. Is 20 % an acceptable mito ceiling given the OXPHOS argument, or would you prefer 15 %?
2. Should S2 be kept as a separate library, or flagged as "low recovery — exploratory only"?
3. How should the 5.6× depth difference of S1 be described in the limitations section?
