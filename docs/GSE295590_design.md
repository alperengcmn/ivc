# GSE295590 — experimental design and replicate structure

Source basis: the original publication (Methods, Figure 3 —
[PMC13382995](https://pmc.ncbi.nlm.nih.gov/articles/PMC13382995/)), GEO records, the SRA run table,
and the initial raw-data inspection.

## Working model

| Level | What we know | Status |
|-------|--------------|--------|
| Experimental animals | 20 WT mice — 10 primary + 10 secondary (Figure 3 scRNA-seq experiment) | ✅ Confirmed |
| Infection | 650 *N. brasiliensis* L3; pyrantel pamoate after primary infection; sampling day 3 after primary / day 3 after secondary | ✅ Confirmed |
| Cell isolation | MACS negative enrichment → FACS Live CD45⁺ Lin⁻ CD90.2⁺ CD127⁺ | ✅ Confirmed |
| nILC2 / iILC2 | Identified **computationally** after sequencing, not sorted separately | ✅ Confirmed |
| GEO samples | 3: GSM8953514 (primary, S1), GSM8953515 (primary, S2), GSM8953516 (secondary, S3) | ✅ Confirmed |
| SRA | 3 BioSamples / 3 Experiments / 6 runs (2 per Experiment) | ✅ Confirmed — validated by Script 02 |
| Library sequencing | Uniquely indexed libraries pooled for sequencing (NovaSeq 6000) | ✅ Confirmed |
| Pre-library mouse pooling | Whether several mice contributed to one 10x library | ❓ **Unresolved** |
| Mouse-to-library mapping | Which mice contributed to each BioSample | ❓ **Unresolved** |
| Biological replicate structure | Independent animal-level replicates | ❓ **Unresolved** |

## Accession crosswalk

| GEO | BioSample | SRA Experiment | SRR runs | Condition |
|-----|-----------|----------------|----------|-----------|
| GSM8953514 | SAMN48137559 | SRX28549631 | SRR33304601, SRR33304602 | primary |
| GSM8953515 | SAMN48137558 | SRX28549632 | SRR33304599, SRR33304600 | primary |
| GSM8953516 | SAMN48137557 | SRX28549633 | SRR33304597, SRR33304598 | secondary |

## Raw-data observations (reproduced in Phase 3 — [QC report](../reports/qc/03_qc_report.md))

| Sample | Cells (first inspection) | Note |
|--------|--------------------------|------|
| S1 · GSM8953514 | ~9,012 | — |
| S2 · GSM8953515 | ~838 | ~10× fewer cells — investigate (recovery, input, library quality); **do not** use it to infer mouse numbers |
| S3 · GSM8953516 | ~8,677 | higher mitochondrial fraction — investigate, do not exclude automatically |

Phase 3 confirmed all three counts exactly and added: matrices pre-filtered at ≈ 500 UMI; **S1 sequenced ~5.6× shallower** than S2/S3 (median 1,283 vs ~7,200 UMI); S3's mito signal is a distinct ~10 % damaged-cell cloud rather than a sample-wide shift.

## What we do **not** assume

- that the three GEO samples are three independent biological replicates;
- that the six SRR runs are six replicates (they are two runs per Experiment);
- any split such as 5 + 5 primary mice across S1 and S2;
- that n = 20 animals means n = 20 libraries;
- mouse pooling, merely because sequencing libraries were pooled.

## Consequence for the analysis

The biological unit is the animal; cells are nested within libraries. Until the
mouse-to-library structure is known, **formal animal-level differential expression is blocked**.
Allowed meanwhile: technical QC, cell-state characterisation, annotation, composition analysis and
explicitly labelled *exploratory* signature discovery.

The question to the authors is drafted in [`author_query.md`](author_query.md).
