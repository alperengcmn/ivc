# Changelog

## [0.3.1] — 2026-10-15

### Fixed
- QC report is now produced entirely by Script 03: it copies its figures to `reports/qc/figures/` and
  embeds them, so the committed report renders on GitHub and a rerun reproduces it exactly.
- `make qc` crashed with `AttributeError: Can only use .cat accessor with a 'category' dtype` on
  older stacks (scanpy 1.10 / anndata 0.10 / pandas 2.2, as resolved by conda-forge): scrublet drops
  categorical dtypes. QC no longer depends on the dtype, restores it after scrublet, and refuses to
  continue if scrublet ever changes cell order.

### Added
- Tests that run the real scrublet path (the earlier suite only used `--no-doublets`), 28 total.
- GitHub Actions CI on an *oldest* and a *latest* dependency stack.

## [0.3.0] — 2026-10-15

### Added
- First Phase 3 run on the deposited GSE295590 matrices; committed report `reports/qc/03_qc_report.md`.
- QC report now detects: pre-filtered matrices (shared UMI floor), between-library depth
  differences, silent scrublet failure (0 % doublets in large libraries), high-mito tails.
- `docs/qc_threshold_proposal.md` — proposed thresholds with dry-run effects and sensitivity table.
- `remove_predicted_doublets` parameter (default `false`, D-009).

### Fixed
- Filter-mode gate is evaluated on per-sample counts *after* filtering (was before).
- Cell counts in report notes print as integers.

## [0.2.0] — 2026-10-15

### Added
- Script 03 — Phase 3 QC (`make qc`): 10x triplet discovery by GSM prefix, per-cell metrics
  (UMI, genes, mito %, ribo %), per-sample scrublet doublets, MAD outlier flags, violin and
  scatter figures, Markdown report. Report-only until thresholds are agreed (D-008).
- Synthetic GSE295590-shaped test data: small S2, mito-rich S3 (5 new tests, 24 total).
- `analysis` optional dependency group; `ilc2-qc` command; VS Code task.

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
