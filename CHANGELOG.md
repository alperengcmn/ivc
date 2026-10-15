# Changelog

## [0.1.0] — 2026-10-15

### Added
- Repository structure following pipeline v2 (config, metadata, docs, scripts, results, reports).
- `ilc2tri` package: strict table reader, shared report helpers, GSE295590 reference module.
- Script 01 — structural sample-metadata validation.
- Script 02 — GEO → BioSample → SRA Experiment → SRR crosswalk validation.
- `scripts/00_fetch_gse295590.sh` — download with checksum and download log (Phase 2).
- Test suite (19 tests), ruff configuration, VS Code workspace settings and tasks.
- Docs: hypothesis & success criteria, dataset hierarchy, GSE295590 design note, pipeline status,
  decision log, author query draft.

### Changed (relative to the scripts in the shared Google Doc)
- Scoped status labels: `STRUCTURAL METADATA VALIDATION` and `GEO-SRA CROSSWALK VALIDATION` (D-007).
- Script 02 checks exact SRR IDs instead of a run count, and words counts as observations (D-006).
- Explicit `--delimiter` now always overrides the file extension (previously ignored for `.csv`/`.tsv` in Script 02).
- Line numbers in error messages use the reader's physical line (correct with quoted newlines).
- Empty replicate columns are reported as "no replicate information recorded yet" rather than
  as populated fields.
