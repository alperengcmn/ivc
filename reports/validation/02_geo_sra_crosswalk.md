# GEO–SRA crosswalk validation — GSE295590

- Input: `metadata/GSE295590/sra_crosswalk.tsv`
- **GEO-SRA CROSSWALK VALIDATION: PASS**

> **Scope.** This status covers the accession relationships GEO → BioSample → SRA Experiment → SRR only. It does **not** validate the underlying biological replicate structure.

## Findings
- Run-level records: **6**
- Required columns: `geo_accession`, `biosample`, `sra_experiment`, `srr_accession`

### Per-sample crosswalk
- `GSM8953514` → `SAMN48137559` → `SRX28549631`: 2 distinct SRR run(s) observed; the current SRA run table lists 2 for this experiment (an observation about sequencing runs, not a replicate count).
- `GSM8953515` → `SAMN48137558` → `SRX28549632`: 2 distinct SRR run(s) observed; the current SRA run table lists 2 for this experiment (an observation about sequencing runs, not a replicate count).
- `GSM8953516` → `SAMN48137557` → `SRX28549633`: 2 distinct SRR run(s) observed; the current SRA run table lists 2 for this experiment (an observation about sequencing runs, not a replicate count).

### SRR runs by SRA Experiment
- `SRX28549631`: 2 distinct SRR run(s) — `SRR33304601`, `SRR33304602`.
- `SRX28549632`: 2 distinct SRR run(s) — `SRR33304599`, `SRR33304600`.
- `SRX28549633`: 2 distinct SRR run(s) — `SRR33304597`, `SRR33304598`.

### Interpretation safeguards
- SRR run counts describe sequencing runs, not biological replicates.
- Mouse-to-library mapping: **UNKNOWN / not established**.
- Pre-library pooling status: **UNKNOWN / not established**.
- Independent biological replicate structure: **UNKNOWN / not inferred**.
- This validation covers accession relationships only; it does not establish the experimental unit.

### Result
- Validation errors: **0**
