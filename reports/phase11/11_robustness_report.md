# Phase 11 — statistical robustness · GSE295590

## PHASE 11 GATE: PASS — all pre-specified sensitivity analyses completed.

A signature gene is *retained* in a scenario when it keeps its frozen sign, FDR < 0.05 and |log2FC| ≥ 0.5 in both contrasts (leave-one-cluster-out: sign and |log2FC| only). *Robust* = retained in ≥ 80% of 9 scenarios.

| scenario | ILC2s | candidates | recall_up retained | recall_down retained | fraction |
|---|---|---|---|---|---|
| mito ≤ 10 % (source-study QC) | 14,855 | 937.0 | 50/50 | 50/50 | 100% |
| ≥ 500 genes per cell | 14,579 | 945.0 | 50/50 | 50/50 | 100% |
| non-cycling ILC2s (proliferation score ≤ 0) | 9,102 | 549.0 | 50/50 | 50/50 | 100% |
| dissociation-stress score ≤ 75th percentile per library | 11,316 | 818.0 | 50/50 | 46/50 | 96% |
| contamination score ≤ library median | 9,143 | 808.0 | 50/50 | 45/50 | 95% |
| contrast A on depth-matched data | 15,088 | 715.0 | 50/50 | 49/50 | 99% |
| leave out cluster 0 (pseudobulk only) | 6,762 |  | 45/50 | 50/50 | 95% |
| leave out cluster 1 (pseudobulk only) | 10,899 |  | 50/50 | 50/50 | 100% |
| leave out cluster 5 (pseudobulk only) | 12,515 |  | 50/50 | 50/50 | 100% |

The non-cycling scenario is partly circular (retention in non-cycling ILC2s was a selection criterion for the frozen signature, Phase 7) and is shown for completeness only.

## Robust genes: 97 of 100
- **recall_up (50)**: Il1r2, Tff1, Fabp5, Cxcl3, Ctla4, Mthfd2, Cst7, Dut, H2-Eb1, Tnfrsf4, H2-DMa, Selenoh, H2-Ab1, Dctpp1, Trerf1, H2-Aa, Ttc39c, Batf, Lgals1, Il13, Siva1, Serpinb6a, Cd74, Prdx4, AW112010, Cisd1, Tomm5, Nudcd2, Pgk1, Tpi1, Arg1, Impa2, Glrx, Mrpl12, Mpp6, Vim, Aprt, Ecm1, Cpd, Hint1, Nhp2, Sh3bgrl, Alad, Nudt21, Erg28, Nme1, Cox5a, Gapdh, Npm1, Mrpl51
- **recall_down (47)**: Csf2, Klf2, Plaur, Spry2, Lpcat2, Tnfrsf9, S1pr1, Klf3, Tle4, Camk2n1, Uba52, Bcl2l11, Bloc1s1, Gm11808, Snhg9, Rnaset2a, Rnaset2b, Fcer1g, Calca, Zfp36, Gem, Cebpb, Pdcd4, Kmt2d, Klf4, Maff, Susd6, Kit, Fos, Malat1, Sdc4, Rab27b, Ets2, B930036N10Rik, Tcf7, Klf10, Cxcr4, H2-Q4, Serpinb9, Rasgrp2, AC149090.1, Socs2, Ifi203, Trim12a, Hexim1, Hes1, Chmp1b

## Dissociation-stress flag
5 signature genes belong to the dissociation-induced program (van den Brink et al., 2017): Klf2, Zfp36, Cebpb, Klf4, Fos (5 of 5 in recall_down, i.e. higher in the primary libraries). The primary libraries carry higher stress scores (Phase 5), so these differences may reflect tissue handling between libraries; interpret the NF-κB/AP-1/TNFα 'down in secondary' results (Phases 7–8) with this caveat.

## What these analyses cannot address
- Library and condition are confounded (one secondary library); no scenario can separate them.
- The mouse-to-library structure remains unknown (D-004), so no animal-level inference is possible.
