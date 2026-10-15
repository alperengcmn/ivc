# Phase 0 — Hypothesis and success criteria

> Pipeline output: `hypothesis_and_success_criteria.md`. This file is the reference point
> for every downstream comparison. Changes require an entry in [`decision_log.md`](decision_log.md).

## Central question

**What distinguishes a trained/persistent ILC2 state from an acutely activated ILC2 state, and
which transcriptional, metabolic and regulatory programs are associated with that state?**

## Working hypotheses

| ID | Hypothesis | Primary data | Tested in |
|----|------------|--------------|-----------|
| H1 | Secondary-response lung ILC2s (GSE295590) carry a reproducible transcriptional program that differs from primary-response ILC2s beyond differences in nILC2/iILC2 composition. | GSE295590 | Phases 3–7 |
| H2 | Part of that program is conserved in independently generated memory-like ILC2s, at pathway level if not gene level. | GSE172258, GSE81700 | Phase 10 |
| H3 | Part of that program is **not** explained by acute activation / tissue repair alone. | GSE149857 | Phase 10 |
| H4 | Conserved features are associated with candidate TF activity and, where data allow, altered chromatin accessibility. | GSE149854, GSE172258 ATAC, GSE291230/1 | Phases 8–9 |
| H5 | Metabolic reprogramming is a consistent component of the trained/memory-like state. | GSE295590, GSE262210 | Phase 7 + context |

The template hypothesis from the pipeline — *"after inflammatory stimulus X, does program Y in lung
ILC2s remain detectable Z days after the stimulus is withdrawn?"* — cannot be tested directly in
GSE295590, because both groups are sampled 3 days after an active exposure. H1 is therefore framed
as a **recall-associated** difference, and persistence is addressed through GSE172258 / GSE81700.

## Success thresholds (to be agreed)

Values are deliberately `null` in [`config/parameters.yaml`](../config/parameters.yaml) until both
collaborators sign off. Gated steps must not run while their threshold is `null`.

| Phase | Criterion | Threshold | Status |
|-------|-----------|-----------|--------|
| 3 · QC | minimum cells per sample | X = **TBD** | set after inspecting S2 (838 cells in first inspection) |
| 3 · QC | mitochondrial fraction | Y = **TBD** | set per dataset after inspecting distributions (S3 elevated) |
| 5 · Annotation | confirmed ILC2 cells | N = **TBD** | — |
| 5 · Annotation | contamination score | M = **TBD** | — |
| 6 · DE | biological replicates | K = **TBD** | **blocked** — replicate structure of GSE295590 unresolved |
| 6 · DE | FDR | < 0.05 | agreed (pipeline v2) |

## Claim ladder

Every result in this project is labelled with the highest rung its evidence supports.

| Rung | Claim | Evidence needed |
|------|-------|-----------------|
| A | *Trained/recall-associated transcriptional state* | Computational identification (this project) |
| B | *Regulatory / epigenetic support* | TF activity + compatible chromatin data |
| C | *Classical trained immunity / causal memory* | Functional, longitudinal, mechanistic evidence — **outside the reach of this reanalysis** |

## Signature rules

1. Candidate genes need: statistical support in GSE295590, meaningful effect size, biological
   coherence, no obvious technical/QC explanation, and consideration of nILC2 vs iILC2 composition.
2. Reproducibility in at least one independent memory-like dataset; pathway-level coherence preferred.
3. "Core" (reproducible) and "supporting" (dataset-specific) features are reported separately.
4. **The discovery signature is frozen (versioned, committed) before any independent validation.**
