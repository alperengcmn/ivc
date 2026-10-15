# ILC2 project — common commands. Run `make help` to list them.
SHELL       := /bin/bash
PY          ?= python
META        := metadata/GSE295590
REPORTS     := reports/validation

.DEFAULT_GOAL := help
.PHONY: help env install validate validate-samples validate-crosswalk test lint format fetch clean

help: ## Show available commands
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-20s\033[0m %s\n",$$1,$$2}'

env: ## Create the conda environment (ilc2-tri)
	conda env create -f environment.yml

install: ## Install the package in editable mode with dev tools (no conda)
	$(PY) -m pip install -e ".[dev]"

validate: validate-samples validate-crosswalk ## Run all Phase 1 metadata validations

validate-samples: ## Script 01 — structural sample-metadata validation
	$(PY) scripts/01_validate_sample_metadata.py $(META)/samples.tsv --report $(REPORTS)/01_sample_metadata.md

validate-crosswalk: ## Script 02 — GEO → BioSample → SRX → SRR crosswalk validation
	$(PY) scripts/02_validate_geo_sra_crosswalk.py $(META)/sra_crosswalk.tsv --report $(REPORTS)/02_geo_sra_crosswalk.md

test: ## Run the test suite
	$(PY) -m pytest

lint: ## Lint with ruff
	ruff check src tests scripts

format: ## Auto-format with ruff
	ruff format src tests scripts

fetch: ## Phase 2 — download GSE295590 processed matrices into data/raw/
	bash scripts/00_fetch_gse295590.sh

clean: ## Remove caches (never touches data/ or results/)
	find . -type d \( -name __pycache__ -o -name .pytest_cache -o -name .ruff_cache \) -prune -exec rm -rf {} +
