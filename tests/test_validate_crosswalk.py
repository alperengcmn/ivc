from conftest import META
from ilc2tri.tabular import read_table
from ilc2tri.validate_crosswalk import STATUS_LABEL, main, render, validate


def run(path):
    return validate(read_table(path))


def test_committed_crosswalk_passes():
    f = run(META / "sra_crosswalk.tsv")
    assert f.passed, f.lines


def test_status_line_and_wording():
    report = render(META / "sra_crosswalk.tsv", run(META / "sra_crosswalk.tsv"))
    assert f"{STATUS_LABEL}: PASS" in report
    assert "not a replicate count" in report
    assert "SRR run counts describe sequencing runs, not biological replicates." in report


def test_swapped_biosample_fails(write_tsv, crosswalk_rows):
    rows = [r[:] for r in crosswalk_rows]
    rows[1][1] = "SAMN48137558"  # belongs to GSM8953515
    f = run(write_tsv(rows))
    assert not f.passed
    assert any("linked to multiple GEO samples" in line for line in f.lines)


def test_missing_run_fails(write_tsv, crosswalk_rows):
    f = run(write_tsv(crosswalk_rows[:-1]))
    assert not f.passed
    assert any("SRR33304598" in line for line in f.lines if "ERROR" in line)


def test_unlisted_run_fails(write_tsv, crosswalk_rows):
    rows = [r[:] for r in crosswalk_rows]
    rows[1][3] = "SRR99999999"
    assert not run(write_tsv(rows)).passed


def test_duplicate_srr_fails(write_tsv, crosswalk_rows):
    assert not run(write_tsv(crosswalk_rows + [crosswalk_rows[1]])).passed


def test_condition_mismatch_fails(write_tsv, crosswalk_rows):
    rows = [r[:] for r in crosswalk_rows]
    rows[5][4] = "primary"  # GSM8953516 run labelled primary
    assert not run(write_tsv(rows)).passed


def test_missing_columns_short_circuits_with_safeguards(write_tsv):
    f = run(write_tsv([["geo_accession", "srr_accession"], ["GSM8953514", "SRR33304601"]]))
    assert not f.passed
    assert any("UNKNOWN" in line for line in f.lines)


def test_cli_exit_codes(tmp_path, write_tsv, crosswalk_rows):
    assert main([str(META / "sra_crosswalk.tsv"), "--report", str(tmp_path / "r.md")]) == 0
    assert main([str(write_tsv(crosswalk_rows[:-1])), "--report", str(tmp_path / "b.md")]) == 1
    assert main([str(tmp_path / "missing.tsv")]) == 2
