from conftest import META
from ilc2tri.tabular import read_table
from ilc2tri.validate_samples import STATUS_LABEL, main, render, validate


def run(path):
    return validate(read_table(path))


def test_committed_metadata_passes():
    f = run(META / "samples.tsv")
    assert f.passed, f.lines


def test_status_line_is_scoped():
    report = render(META / "samples.tsv", run(META / "samples.tsv"))
    assert f"{STATUS_LABEL}: PASS" in report
    assert "does **not** establish the biological replicate structure" in report


def test_replicate_structure_is_never_inferred():
    text = "\n".join(run(META / "samples.tsv").lines)
    assert "UNKNOWN / not inferred" in text
    assert "UNRESOLVED" in text


def test_duplicate_geo_fails(write_tsv, samples_rows):
    rows = samples_rows + [samples_rows[1]]
    assert not run(write_tsv(rows)).passed


def test_group_conflict_fails(write_tsv, samples_rows):
    rows = [r[:] for r in samples_rows]
    rows[3][2] = "primary"  # GSM8953516 is secondary
    f = run(write_tsv(rows))
    assert not f.passed
    assert any("conflicts" in line for line in f.lines)


def test_invalid_group_reported_once(write_tsv, samples_rows):
    rows = [r[:] for r in samples_rows]
    rows[1][2] = "tertiary"
    f = run(write_tsv(rows))
    assert f.errors == 1


def test_malformed_geo_id_fails(write_tsv, samples_rows):
    rows = [r[:] for r in samples_rows]
    rows[1][0] = "SRR33304601"
    assert not run(write_tsv(rows)).passed


def test_missing_required_column_fails(write_tsv):
    assert not run(write_tsv([["geo_accession"], ["GSM8953514"]])).passed


def test_ragged_row_fails(write_tsv, samples_rows):
    rows = samples_rows + [["GSM8953599", "x"]]
    assert not run(write_tsv(rows)).passed


def test_cli_exit_codes(tmp_path, write_tsv, samples_rows):
    report = tmp_path / "r.md"
    assert main([str(META / "samples.tsv"), "--report", str(report)]) == 0
    assert report.read_text().startswith("# Sample metadata validation")
    bad = write_tsv(samples_rows + [samples_rows[1]])
    assert main([str(bad), "--report", str(tmp_path / "bad.md")]) == 1
    assert main([str(tmp_path / "does_not_exist.tsv")]) == 2
