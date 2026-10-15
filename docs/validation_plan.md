# Pre-registered validation plan (Phases 9–10)

**Written and committed before any validation-dataset expression or accessibility values were
analysed** (only file headers and GEO sample metadata had been inspected). The frozen signature
is `signatures/GSE295590_recall_ILC2_signature_v1.tsv` (SHA-256 in `signatures/MANIFEST.tsv`).
Any change to these rules after results are seen must be logged as a deviation in
[`decision_log.md`](decision_log.md) and reported next to the original result.

## Two signature metrics, used for every dataset

1. **Score** — per sample unit: mean z-score (log expression, z across the dataset's units) of
   `recall_up` genes minus mean z of `recall_down` genes. Compared between groups with a
   one-sided Mann–Whitney test in the expected direction (memory/recall > comparator).
2. **Gene direction agreement** — for each signature gene measured and expressed in the
   dataset, the sign of its log2 fold change (memory/recall vs comparator). Agreement = fraction
   with the expected sign; one-sided binomial test against 0.5.

A dataset is *not testable* when fewer than 20 of the 100 signature genes are measurable.

## Core memory/recall datasets — verdict categories

| Dataset | Comparison (expected higher score) | Unit |
|---|---|---|
| GSE172258 RNA | **1B** (Alternaria-sensitised, after recall) vs **2B** (saline, after recall) — primary test | 20/40-cell well |
| GSE172258 RNA | 1A (sensitised, before recall) vs 2A (saline, before recall) — resting memory | well |
| GSE81700 | 4 months vs naive; week 2 vs naive; day 4 vs naive (lung only) | mouse sample |

- **Replicated** — score test p < 0.05 **and** gene agreement p < 0.05, both in the expected direction.
- **Partially supported** — exactly one of the two is significant in the expected direction, or the
  score test is not computable (< 3 units per group) and gene agreement is significant.
- **Not supported** — neither is significant, or either is significant in the opposite direction.
- **Not testable with available data** — < 20 signature genes measurable.

Caveat fixed in advance: GSE172258 wells are pooled from the same animals per group, so well-level
p-values are exploratory; GSE81700 has 2–3 mice per group.

## Acute-activation reference (GSE149857, influenza, scRNA-seq day 5 / day 10)

Not a verdict: each signature gene is annotated **acute-shared** if, within lung ILC2s, it
correlates with an activation axis (type-2 effector + activation + proliferation scores, computed
*without* signature genes) in the signature's direction (Spearman |ρ| ≥ 0.1, FDR < 0.05) on both days.
Genes that are not acute-shared are reported as *recall-specific candidates*.

## Mechanistic / context datasets (no verdict; enrichment only)

| Dataset | Question | Test |
|---|---|---|
| GSE149854 | Are `recall_up` genes BATF-dependent? | Fisher: fraction down in *Batf* KO (padj < 0.05) vs all genes |
| GSE262210 | Do metabolic `recall_up` genes respond to iron deprivation? | Fisher: fraction changed by DEF (padj < 0.05) vs all genes; direction reported |
| GSE291231 | Activation-shared (stim vs unstim) and Blimp-1 dependence | Gene direction agreement as above; Fisher for KO |

## Chromatin (Phase 9) — GSE172258 ATAC-seq

DESeq2 (pydeseq2) per comparison: Rag1⁻/⁻ Alt/Alt vs Sal/Alt (3 v 3), WT Alt/Alt vs Sal/Alt (3 v 3),
Rag1⁻/⁻ Alt/– vs Sal/– (2 v 3). Peaks are assigned to the nearest RefSeq TSS (mm10) within 50 kb; a
gene's accessibility change is the mean Wald statistic of its peaks.

- **Regulatory support** — `recall_up` genes have higher gene-level accessibility statistics than all
  other genes (one-sided Mann–Whitney p < 0.05) in **both** recall comparisons (Rag1⁻/⁻ and WT).
- **Partial** — in one of them. **None** — in neither.
- Motif enrichment requires a genome sequence and motif database that this project does not ship;
  it is *not performed* and stated as such in the final report (pipeline Phase 9 fallback).
