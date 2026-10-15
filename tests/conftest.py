import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

META = ROOT / "metadata" / "GSE295590"


@pytest.fixture
def write_tsv(tmp_path):
    """Write rows (lists of strings) as a TSV and return its path."""

    def _write(rows, name="table.tsv"):
        path = tmp_path / name
        path.write_text("\n".join("\t".join(r) for r in rows) + "\n", encoding="utf-8")
        return path

    return _write


@pytest.fixture
def samples_rows():
    return [line.split("\t") for line in (META / "samples.tsv").read_text().splitlines()]


@pytest.fixture
def crosswalk_rows():
    return [line.split("\t") for line in (META / "sra_crosswalk.tsv").read_text().splitlines()]
