# Sample metadata validation — GSE295590

- Input: `metadata/GSE295590/samples.tsv`
- **STRUCTURAL METADATA VALIDATION: PASS**

> **Scope.** This status covers the *structure* of the sample metadata only (identifiers, group labels, completeness). It does **not** establish the biological replicate structure, the mouse-to-library mapping, or pre-library pooling.

## Findings
- Records read: **3**
- Columns: `geo_accession`, `sample_title`, `infection_group`, `mouse_id`, `biological_replicate`
- Required columns: `geo_accession`, `infection_group`
- GEO sample set matches the three deposited GSE295590 samples.
- Infection-group labels: populated, valid, and consistent with deposited sample titles.
- Blank cells by column: `mouse_id`: 3; `biological_replicate`: 3.

### Replicate metadata
- Replicate field `mouse_id`: 0/3 record(s) populated.
- Replicate field `biological_replicate`: 0/3 record(s) populated.
- **UNRESOLVED:** Replicate fields exist but are empty: no replicate information has been recorded yet.
- Mouse-to-library mapping: **UNKNOWN / not established**.
- Pre-library pooling status: **UNKNOWN / not established**.
- Independent biological replicate structure: **UNKNOWN / not inferred**.

### Result
- Validation errors: **0**
