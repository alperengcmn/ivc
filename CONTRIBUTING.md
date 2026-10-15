# Contributing

This is a two-person collaboration (Alperen Göçmen · Divya). These rules keep the analysis
reviewable and the claims honest.

## Workflow

1. **Branch per task** — `phase3/qc-gse295590`, `docs/author-reply`, `fix/crosswalk-wording`.
2. **Pull request for everything that reaches `main`.** The other collaborator reviews.
   Scientific review (is the interpretation right?) counts as much as code review.
3. **Commits** use a short prefix: `feat:`, `fix:`, `docs:`, `test:`, `data:`, `refactor:`.
4. Before opening a PR:
   ```bash
   make lint test validate
   ```

## Scientific rules (non-negotiable)

- **Cells are not replicates.** Inference uses the sample/animal as the unit.
- **Never infer** mouse identity, pooling or replicate assignments that no source states.
  Unknown stays `UNKNOWN` and is written down.
- **Label every claim** with the claim ladder (A / B / C) from
  [`docs/00_hypothesis_and_success_criteria.md`](docs/00_hypothesis_and_success_criteria.md).
  Use "trained/recall-associated" — never "proves trained immunity".
- **Thresholds are agreed, not tuned.** Change one only with a [decision log](docs/decision_log.md) entry.
- **Freeze before you validate.** The discovery signature is committed before validation data is touched.

## Code rules

- Logic in `src/ilc2tri/`, thin CLIs in `scripts/NN_name.py`, tests in `tests/`.
- Every script: `--help`, exit codes `0` pass / `1` fail / `2` I/O error, a Markdown report.
- Parameters in `config/`, never hard-coded in notebooks.
- Fixed random seed (`config/parameters.yaml → random_seed`).
- Notebooks are for exploration; anything a result depends on moves into a script.
  Clear outputs before committing.

## Data rules

- `data/` and `results/` are git-ignored — everything is re-creatable from accessions + scripts.
- Small, reviewable evidence (validation reports, metadata tables, frozen signatures) **is** committed.
- Every download is logged with URL, SHA-256 and UTC date.
