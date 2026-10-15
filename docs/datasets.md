# Datasets

Machine-readable registry: [`config/datasets.yaml`](../config/datasets.yaml).

## Hierarchy

| Dataset | Role | Priority | Model / assay |
|---------|------|----------|---------------|
| **GSE295590** | Discovery — trained/recall reference | 🔴 Core | *N. brasiliensis* primary vs secondary, lung ILC2, scRNA-seq |
| **GSE172258** | Independent memory validation (+ chromatin) | 🔴 Core | Repeated *Alternaria* exposure, scRNA-seq + ATAC-seq |
| **GSE149857** | Acute inflammatory reference | 🔴 Core | Influenza, lung ILCs, scRNA-seq |
| GSE81700 | Memory-like validation, signature level only | 🟠 High | Naive / IL-33-activated / long-lived ILC2, microarray |
| GSE149854 | BATF-dependent regulatory context | 🟠 High | Influenza, WT vs *Batf* KO, bulk RNA-seq (n = 3/group) |
| GSE262210 | Metabolic context | 🟡 Supporting | Iron availability and ILC2 metabolism |
| GSE291230 / GSE291231 | Additional regulatory reference | 🟡 Supporting | IL-9 / Blimp-1, paired ATAC-seq / RNA-seq |

## How each dataset answers a question

```
Q1  What does the trained state look like?        GSE295590
Q2  Which features are memory/training-associated? GSE172258 + GSE81700
Q3  Which features are just acute activation?      GSE149857
Q4  What might regulate the state?                 GSE149854 + GSE172258 ATAC + GSE291230/1
Q5  Is metabolism part of the phenotype?           GSE295590 + GSE262210
```

## Caveats that travel with each dataset

- **GSE295590** — primary vs secondary are both sampled at day 3 after an active exposure, so the
  comparison mixes prior exposure with the rechallenge response. Replicate structure unresolved;
  see [`GSE295590_design.md`](GSE295590_design.md).
- **GSE172258** — allergen vs helminth model: look for conserved programs, not identical gene lists.
- **GSE81700** — microarray: gene/signature-level validation only, no expression-value integration.
- **GSE149854, GSE262210, GSE291230/1** — context and mechanism, not evidence of trained immunity.

Datasets are never merged directly when species, tissue, condition, technology or isolation
strategy are incompatible (pipeline Phase 1, rule 3).
