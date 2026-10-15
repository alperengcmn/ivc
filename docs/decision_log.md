# Decision log

Every decision that changes scope, a dataset's role, a threshold or a statistical choice gets an
entry. Newest first. Format: **ID · date · decision** — context — consequence.

---

### D-008 · 2026-10-15 · Phase 3 QC runs report-only until thresholds are agreed
- **Context:** Phase 0 thresholds X/Y are not yet agreed; S2 is small and S3 mito-rich.
- **Decision:** Script 03 removes no cells while any `phase3_qc` threshold is `null`. It flags
  per-sample MAD outliers (5 MAD counts/genes on log scale, 3 MAD upper for mito %) and reports
  `suggested_max_pct_mt` as a discussion starting point only. Thresholds are set per dataset after
  both collaborators review the report; once set, the same script filters and evaluates the gate.
- **Consequence:** no threshold is tuned to a result; S2/S3 are investigated, never silently dropped.

### D-007 · 2026-10-15 · Validators report scoped status labels
- **Context:** a bare "PASS" could be read as "biological design validated" (Divya's review).
- **Decision:** Script 01 reports `STRUCTURAL METADATA VALIDATION: PASS/FAIL`; Script 02 reports
  `GEO-SRA CROSSWALK VALIDATION: PASS/FAIL`. Each report opens with a scope note.
- **Consequence:** a passing report can never be quoted as evidence of replicate structure.

### D-006 · 2026-10-15 · SRR run counts are observations, not expectations about replication
- **Context:** Script 02 originally checked "expected 2" runs per Experiment.
- **Decision:** Script 02 now checks the exact SRR IDs listed in the current SRA run table and
  words the count as "an observation about sequencing runs, not a replicate count".
- **Consequence:** stronger check (exact IDs, not just counts) with unambiguous wording.

### D-005 · 2026-10-15 · Crosswalk validation split from sample validation
- **Decision:** Script 01 = GEO sample-level structure; Script 02 = GEO → BioSample → SRX → SRR.
- **Consequence:** each script has one responsibility; both share one reference module
  (`src/ilc2tri/gse295590.py`).

### D-004 · 2026-10-15 · Formal DE blocked until replicate structure is known
- **Context:** 20 mice vs 3 deposited samples; mapping and pooling unresolved.
- **Decision:** no animal-level DE inference on GSE295590 until resolved or documented as a
  limitation with an agreed alternative. Exploratory analyses are labelled as such.
- **Next action:** send [`author_query.md`](author_query.md).

### D-003 · 2026-10-15 · Signature frozen before validation
- **Decision:** the discovery signature is versioned and committed before touching validation datasets.

### D-002 · 2026-10-15 · Dataset hierarchy adopted
- **Decision:** GSE295590 discovery; GSE172258 + GSE81700 memory validation; GSE149857 acute
  reference; GSE149854 + ATAC sets regulatory; GSE262210 metabolic. See [`datasets.md`](datasets.md).

### D-001 · 2026-10-15 · Claim ladder adopted
- **Decision:** results are labelled A (trained/recall-associated state), B (regulatory support)
  or C (classical trained immunity). This project targets A and, where data allow, B.
