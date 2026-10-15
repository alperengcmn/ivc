# Lung ILC2 trained/recall-associated programs — final report

*Generated automatically by `scripts/12_final_report.py` on 2026-10-15 10:18 UTC from commit `94dce71`. Every number below is read from `results/`. Thresholds marked PROVISIONAL await Divya's review (D-010).*

## Summary

- Secondary-response lung ILC2s from *N. brasiliensis* rechallenge (GSE295590) carry a reproducible program: **883 genes** change concordantly in two contrasts that share no primary cells (effect sizes r = 0.91); the published findings for **Arg1, Cox5a and Pkm** are reproduced.
- The program is metabolic and proliferative (up: OXIDATIVE_PHOSPHORYLATION, MYC_TARGETS_V1, E2F_TARGETS, DNA_REPAIR, MTORC1_SIGNALING, G2M_CHECKPOINT) with lower HYPOXIA, ESTROGEN_RESPONSE_EARLY, TNFA_SIGNALING_VIA_NFKB signalling; inferred TF activity: Myc up; Rela, Stat5a, Runx1, Stat1, Ets1, Stat2 down.
- A 100-gene signature was frozen before validation (SHA-256 `537d38d95778…`). It is robust internally (97/100 genes survive ≥ 80 % of sensitivity scenarios, including the source study's own QC).
- Independent evidence is **mixed**: chromatin — *Partial regulatory support* (strong in wild-type allergen memory, absent in Rag1⁻/⁻); memory transcriptome — GSE172258 RNA after recall: Alt-sensitised (1B) vs saline (2B): *Not supported*; GSE172258 RNA before recall: Alt-sensitised (1A) vs saline (2A): *Not supported*; GSE81700 lung ILC2 4 months after IL-33 vs naive: *Not supported*; GSE81700 lung ILC2 week 2 after IL-33 vs naive: *Partially supported*; GSE81700 lung ILC2 day 4 after IL-33 vs naive: *Partially supported*.
- **Interpretation:** the signature tracks *recent* type-2 activation — it is strongly shared with IL-33-driven activation and fades over weeks — so the secondary-response program looks like an intensified activation state rather than a persistent memory signature.
- **Claim level reached: A** — a recall-associated transcriptional state with partial regulatory support (B, wild-type only). Classical trained immunity (C) is not claimed.

## 1. Research question
What distinguishes a trained/persistent ILC2 state from an acutely activated ILC2 state, and which transcriptional, metabolic and regulatory programs are associated with it? Hypotheses H1–H5: [`docs/00_hypothesis_and_success_criteria.md`](../docs/00_hypothesis_and_success_criteria.md).

| Hypothesis | Outcome |
|---|---|
| H1 — reproducible recall program in GSE295590 beyond composition | **Supported (exploratory)** — 883 concordant genes; 97/100 signature genes robust; retained in non-cycling ILC2s |
| H2 — conserved in independent memory-like ILC2s | **Partially supported** (pre-registered rule) — support comes from recently activated ILC2s (IL-33 day 4: 69/93 genes agree, *Partially supported*; week 2: 55/93 genes agree, *Partially supported*), not from long-lived memory (4 months: 46/93 genes agree, *Not supported*; Alternaria 1 week after recall: 36/80 genes agree, *Not supported*) |
| H3 — not explained by acute activation alone | **Partly** — not shared with influenza-driven activation (0 acute-shared genes; 3 opposite), but strongly shared with IL-33-driven activation (GSE81700 day 4: 69/93 genes agree, *Partially supported*) |
| H4 — TF / chromatin support | **Partial** — 14 target-supported TFs; ATAC: partial regulatory support |
| H5 — metabolic reprogramming | **Supported in discovery** — OXPHOS, glycolysis, MYC, mTORC1 up; signature responds to iron deprivation (§10) |

## 2. Datasets and design
Roles and caveats: [`docs/datasets.md`](../docs/datasets.md). Discovery: GSE295590 — 20 WT mice (10 primary, 10 secondary, day 3), FACS-sorted lung ILCs, three 10x libraries (S1, S2 primary; S3 secondary). **The mouse-to-library mapping is unknown (D-004)**; the secondary condition is one library.

## 3. Quality control (PROVISIONAL thresholds, D-010)
mito ≤ 20 %, 200–8000 genes per cell, ≥ 500 cells per library; predicted doublets not removed (scrublet threshold failed, D-009). 17,575 cells retained.

| library | n_cells | median_total_counts | median_genes | median_pct_mt | pct_cells_mt_gt20 |
|---|---|---|---|---|---|
| GSM8953514 | 9,012 | 1,283 | 732 | 2.0 | 0.0 |
| GSM8953515 | 838 | 7,160 | 2,386 | 2.3 | 1.4 |
| GSM8953516 | 8,677 | 7,215 | 2,095 | 3.7 | 10.1 |

Key technical findings: matrices pre-filtered at ≈ 500 UMI; **`GSM8953514` (S1) sequenced 5.6× shallower** than the deepest library; S1 ILC2-like cells carry elevated T-lineage, haemoglobin and immediate-early transcripts (ambient RNA / dissociation stress, D-012).

## 4. Embedding and annotation
Embedding on binomially thinned (depth-matched) counts — PC1–depth |r| 0.77 → 0.17; no cross-library batch correction because library and condition are confounded (D-011). Clustering resolution chosen by a pre-specified stability rule (min ARI ≥ 0.8 over 5 seeds). Lineage by marker categories (never a single marker).

| lineage | GSM8953514 | GSM8953515 | GSM8953516 |
|---|---|---|---|
| B_cell | 4 | 29 | 31 |
| ILC2 | 8245 | 221 | 6622 |
| NK | 59 | 355 | 498 |
| T_cell | 689 | 149 | 262 |
| mast | 12 | 3 | 64 |
| unassigned | 2 | 54 | 276 |

Reproduction of the published cell numbers: S2 contains 221 ILC2s (197 at the paper's mito ≤ 10 %) vs **189 primary ILC2s published** — the published primary group most likely comprised S2 only, and S1 (8,245 ILC2-like cells) appears absent from the published ILC2 analysis. nILC2/iILC2 were **not resolvable** as stable discrete sub-clusters in the clean reference library; the iILC2 program is carried as a continuous score (D-012).

## 5. Statistical methods
Exploratory library-consistency design (D-013): contrast A = S3 vs S2 (full depth), contrast B = S3 vs S1 (depth-matched). Wilcoxon tests rank genes only (cells are not replicates); effect sizes are pseudobulk log2FC; candidates must agree in both contrasts (FDR < 0.05, |log2FC| ≥ 0.5, same sign, detected in ≥ 10 % of cells) and are not contamination-sensitive. Validation rules were pre-registered before any validation data were analysed (D-014, [`docs/validation_plan.md`](../docs/validation_plan.md)).

## 6. Differential expression
883 concordant genes (498 up, 385 down in secondary). Genes named by the source study:

| gene | pb_log2fc_A | fdr_A | pb_log2fc_B | fdr_B | candidate |
|---|---|---|---|---|---|
| Arg1 | 1.48 | 1.27e-47 | 1.92 | 0 | ✓ |
| Cox5a | 1.64 | 2.45e-66 | 1.23 | 5.82e-200 | ✓ |
| Pkm | 1.14 | 4.8e-34 | 1.11 | 5.08e-132 | ✓ |
| Areg | 0.338 | 1.9e-13 | 0.228 | 3.01e-42 |  |
| Tph1 |  |  |  |  |  |
| Il2ra | 0.03 | 0.626 | 1.21 | 2.62e-117 |  |
| Csf2 | -4.88 | 7.15e-30 | -5.13 | 0 | ✓ |

Other genes discussed by the source study (log2FC contrast A / B): Areg +0.34 / +0.23 (same direction, below threshold); Csf2 -4.88 / -5.13 (concordant candidate); Il2ra +0.03 / +1.21 (same direction, below threshold); Tph1: not tested. The paper reports Csf2 as unchanged.

## 7. Pathways and frozen signature
Concordant Hallmark pathways — up: OXIDATIVE_PHOSPHORYLATION, MYC_TARGETS_V1, E2F_TARGETS, DNA_REPAIR, MTORC1_SIGNALING, G2M_CHECKPOINT, FATTY_ACID_METABOLISM, GLYCOLYSIS, PEROXISOME; down: COMPLEMENT, IL6_JAK_STAT3_SIGNALING, KRAS_SIGNALING_UP, TGF_BETA_SIGNALING, KRAS_SIGNALING_DN, INFLAMMATORY_RESPONSE, HYPOXIA, ESTROGEN_RESPONSE_EARLY, TNFA_SIGNALING_VIA_NFKB. Reactome: 45 concordant.

**Frozen signature v1** (`signatures/GSE295590_recall_ILC2_signature_v1.tsv`, SHA-256 `537d38d957780a607e1632b74b4c7060801b02cf7609ea4260862d70d3bcec1a`):
- recall_up: Il1r2, Tff1, Fabp5, Cxcl3, Ctla4, Mthfd2, Cst7, Dut, H2-Eb1, Tnfrsf4, H2-DMa, Selenoh, H2-Ab1, Dctpp1, Trerf1, H2-Aa, Ttc39c, Batf, Lgals1, Il13, Siva1, Serpinb6a, Cd74, Prdx4, AW112010, Cisd1, Tomm5, Nudcd2, Pgk1, Tpi1, Arg1, Impa2, Glrx, Mrpl12, Mpp6, Vim, Aprt, Ecm1, Cpd, Hint1, Nhp2, Sh3bgrl, Alad, Nudt21, Erg28, Nme1, Cox5a, Gapdh, Npm1, Mrpl51
- recall_down: Csf2, Klf2, Plaur, Spry2, Lpcat2, Tnfrsf9, S1pr1, Klf3, Isg15, Tle4, Camk2n1, Uba52, Bcl2l11, Bloc1s1, Gm11808, Snhg9, Rnaset2a, Rnaset2b, Fcer1g, Calca, Zfp36, Gem, Cebpb, Pdcd4, Kmt2d, Klf4, Maff, Susd6, Kit, Fos, Malat1, Sdc4, Rab27b, Ets2, B930036N10Rik, Tcf7, Klf10, Cxcr4, Jag1, H2-Q4, Serpinb9, Rasgrp2, AC149090.1, Socs2, Ifi203, Trim12a, Ccl4, Hexim1, Hes1, Chmp1b

## 8. Transcription-factor and pathway activity (inference, not binding evidence)
Target-supported TFs — more active in secondary: Myc; less active: Rela, Stat5a, Runx1, Stat1, Ets1, Stat2, Nfkb1, Creb1, Irf1, Jund, Foxo3, Nfkb2, Fos. PROGENy: Hypoxia ↓, NFkB ↓, p53 ↓, TNFa ↓, TGFb ↓, JAK-STAT ↓, PI3K ↑, EGFR ↑, MAPK ↑, Estrogen ↑.
TFs named by the GSE172258 study that are concordant and target-supported here: **Jund**.

## 9. Epigenomic integration and cross-dataset validation
**Chromatin (GSE172258 ATAC, Phase 9): Partial regulatory support** (pre-registered rule).

| comparison | n | DA_peaks_FDR5 | up_genes | p_up_gain | down_genes | p_down_loss |
|---|---|---|---|---|---|---|
| Rag1⁻/⁻ Alt/Alt vs Sal/Alt | 3 v 3 | 4131 | 49 | 0.353 | 47 | 0.119 |
| WT Alt/Alt vs Sal/Alt | 3 v 3 | 24735 | 50 | 1.99e-09 | 47 | 2.91e-06 |
| Rag1⁻/⁻ Alt/– vs Sal/– (resting memory) | 2 v 3 | 1829 | 49 | 0.35 | 47 | 0.657 |

**Memory transcriptomes (Phase 10, pre-registered verdicts):**

| dataset | comparison | score_delta | score_p | genes_agree | genes_n_genes | genes_p | verdict |
|---|---|---|---|---|---|---|---|
| GSE172258 RNA | after recall: Alt-sensitised (1B) vs saline (2B) | -0.0287 | 0.872 | 36 | 80 | 0.843 | Not supported |
| GSE172258 RNA | before recall: Alt-sensitised (1A) vs saline (2A) | 0.0166 | 0.312 | 44 | 80 | 0.217 | Not supported |
| GSE81700 | lung ILC2 4 months after IL-33 vs naive | 0.203 |  | 46 | 93 | 0.582 | Not supported |
| GSE81700 | lung ILC2 week 2 after IL-33 vs naive | 0.311 |  | 55 | 93 | 0.0483 | Partially supported |
| GSE81700 | lung ILC2 day 4 after IL-33 vs naive | 1.86 |  | 69 | 93 | 1.65e-06 | Partially supported |

Interpretation, not part of the verdict: GSE172258 wells were sampled one week after recall versus three days in the discovery data; an early recall program may have subsided. The strong wild-type ATAC result versus the null Rag1⁻/⁻ result suggests the chromatin component may depend on adaptive-immune context.

## 10. Acute activation, mechanism and metabolism
Influenza (GSE149857): no signature gene is acute-shared on both days; 3 move opposite to acute activation (Calca, Csf2, Plaur) — they distinguish recall from acute activation.

| dataset | question | in_set | background | odds_ratio | p |
|---|---|---|---|---|---|
| GSE149854 | recall_up genes reduced in Batf KO (BATF-dependent) | 3/31 | 3752/6369 | 0.0747 | 1 |
| GSE262210 | recall_up genes changed by iron deprivation (DEF vs vehicle) | 47/50 | 9078/13532 | 7.69 | 5.83e-06 |
| GSE262210 | recall_down genes changed by iron deprivation (DEF vs vehicle) | 32/46 | 9093/13536 | 1.12 | 0.433 |
| GSE291231 | signature direction under in-vitro stimulation (activation-shared?) | 59/95 agree | - |  | 0.0117 |
| GSE291231 | signature genes changed in Blimp-1 KO | 24/98 | 1076/16277 | 4.58 | 2.09e-08 |

## 11. Robustness (Phase 11)
| scenario | ILC2s | recall_up retained | recall_down retained | fraction |
|---|---|---|---|---|
| mito ≤ 10 % (source-study QC) | 14855 | 50/50 | 50/50 | 100 % |
| ≥ 500 genes per cell | 14579 | 50/50 | 50/50 | 100 % |
| non-cycling ILC2s (proliferation score ≤ 0) | 9102 | 50/50 | 50/50 | 100 % |
| dissociation-stress score ≤ 75th percentile per library | 11316 | 50/50 | 46/50 | 96 % |
| contamination score ≤ library median | 9143 | 50/50 | 45/50 | 95 % |
| contrast A on depth-matched data | 15088 | 50/50 | 49/50 | 99 % |
| leave out cluster 0 (pseudobulk only) | 6762 | 45/50 | 50/50 | 95 % |
| leave out cluster 1 (pseudobulk only) | 10899 | 50/50 | 50/50 | 100 % |
| leave out cluster 5 (pseudobulk only) | 12515 | 50/50 | 50/50 | 100 % |

## 12. Limitations
1. **One secondary library**: condition and library are confounded; no analysis here can separate them.
2. **Unknown mouse structure** (D-004): all inference is exploratory; p-values rank, they do not test animals.
3. **S1 is technically distinct** (5.6× shallower, ambient RNA, stress) and likely absent from the published analysis.
4. **Dissociation stress**: part of the 'down in secondary' program (Fos, Klf2/4, Zfp36, Cebpb; NF-κB/AP-1) may reflect tissue handling.
5. **Cross-model validation**: allergen (GSE172258), IL-33 (GSE81700) and influenza (GSE149857) differ from helminth recall; timing differs (day 3 vs 1 week).
6. **Provisional QC thresholds** (D-010); results were stable to the source study's stricter QC (§11).
7. Motif enrichment was not performed (no genome sequence / motif database). GSE291230 (Blimp-1 ATAC) is deposited as bigwig tracks only (~0.7 GB each) and was not integrated; its paired RNA (GSE291231) was.

## 13. Conclusions supported by the data
- Secondary-response lung ILC2s show a reproducible, metabolically active and proliferative recall program with reduced inflammatory/AP-1 signalling — **claim level A (recall-associated state)**.
- It differs from influenza-driven activation but closely resembles IL-33-driven activation (day 4), weakening by week 2 and absent at 4 months; it is not reproduced in the Alternaria memory transcriptome one week after recall.
- Chromatin at signature genes opens in wild-type allergen-memory ILC2s (not Rag1⁻/⁻) — **partial level B**.
- Together, the data support an **intensified type-2 activation / recall state** in secondary ILC2s; they do not show that this state persists as a resting memory program.
- Classical trained immunity or epigenetic memory is **not** claimed.

## 14. Reproducibility
- Code commit `94dce71`; environment `environment.yml`; software: python 3.13.16, scanpy 1.12.4, anndata 0.13.4, numpy 2.5.3, pandas 3.0.5, scipy 1.18.1, decoupler 2.2.0, gseapy 1.3.1, pydeseq2 0.5.4, leidenalg 0.12.0.
- Random seed 0; all parameters in `config/`; decisions in `docs/decision_log.md`; downloads with SHA-256 in `data/raw/*/manifest.tsv` and `references/manifest.tsv`.
- Re-run everything: `make all` (raw data → this report).
