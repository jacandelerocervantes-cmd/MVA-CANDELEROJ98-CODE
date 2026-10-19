"""Tests for Track 1 submission CSV validation against variant_ranking schema."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import jsonschema
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
PROJECT_ROOT = REPO_ROOT.parent
SCHEMAS_DIR = REPO_ROOT / "shared-contracts" / "schemas"
GATK_CSV = PROJECT_ROOT / "submissions" / "track1" / "candeleroj98_gatk-manta-fusion.csv"
CSV_TARGETS = [GATK_CSV]

# The submission CSV is patient-derived and is not versioned (see .gitignore): these tests run only where it exists locally.
pytestmark = [
    pytest.mark.track1,
    pytest.mark.skipif(not GATK_CSV.exists(), reason="submission CSV is patient-derived and not published"),
]


@pytest.fixture
def schema_validator() -> jsonschema.Draft202012Validator:
    """Loads and compiles Draft202012Validator for variant_ranking.schema.json."""
    schema_path = SCHEMAS_DIR / "variant_ranking.schema.json"
    assert schema_path.exists(), f"Schema not found: {schema_path}"
    with open(schema_path, "r", encoding="utf-8") as f:
        schema = json.load(f)
    return jsonschema.Draft202012Validator(schema)


@pytest.mark.parametrize("csv_path", CSV_TARGETS)
def test_submission_csv_exists_and_row_count(csv_path: Path) -> None:
    """Verifies that the CSV exists and contains between 1 and 10 predictions."""
    assert csv_path.exists(), f"Submission CSV missing at {csv_path}"
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    assert 1 <= len(rows) <= 10, f"Expected 1-10 rows, found {len(rows)}"


@pytest.mark.parametrize("csv_path", CSV_TARGETS)
def test_submission_csv_columns(csv_path: Path) -> None:
    """Verifies that all required header columns are present."""
    expected_cols = [
        "proband_id",
        "chrom_1",
        "pos_1",
        "ref_1",
        "alt_1",
        "chrom_2",
        "pos_2",
        "ref_2",
        "alt_2",
        "epcr",
        "finding_type",
        "notes",
    ]
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
    assert header == expected_cols


@pytest.mark.parametrize("csv_path", CSV_TARGETS)
def test_submission_csv_epcr_monotonic_ranking(csv_path: Path) -> None:
    """Verifies that EPCR values are in (0, 1] and strictly descending."""
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        epcr_vals = [float(row["epcr"]) for row in reader]

    for val in epcr_vals:
        assert 0.0 < val <= 1.0, f"EPCR {val} is outside (0, 1]"

    # Verify strictly monotonic descending order
    for i in range(len(epcr_vals) - 1):
        assert epcr_vals[i] > epcr_vals[i + 1], (
            f"EPCR ranking not strictly descending: {epcr_vals[i]} <= {epcr_vals[i+1]}"
        )


@pytest.mark.parametrize("csv_path", CSV_TARGETS)
def test_submission_csv_schema_conformance(
    csv_path: Path,
    schema_validator: jsonschema.Draft202012Validator,
) -> None:
    """Parses and validates every row against the JSON Schema Draft 2020-12."""
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader):
            instance: dict[str, Any] = {
                "proband_id": row["proband_id"],
                "chrom_1": row["chrom_1"],
                "pos_1": int(row["pos_1"]),
                "ref_1": row["ref_1"],
                "alt_1": row["alt_1"],
                "epcr": float(row["epcr"]),
                "finding_type": row["finding_type"],
            }
            if row["notes"]:
                instance["notes"] = row["notes"]

            # Compound-het handling
            has_chrom_2 = bool(row["chrom_2"].strip())
            has_pos_2 = bool(row["pos_2"].strip())
            has_ref_2 = bool(row["ref_2"].strip())
            has_alt_2 = bool(row["alt_2"].strip())

            if has_chrom_2 or has_pos_2 or has_ref_2 or has_alt_2:
                assert has_chrom_2 and has_pos_2 and has_ref_2 and has_alt_2, (
                    f"Row {idx+1} has partial compound-het fields"
                )
                instance["chrom_2"] = row["chrom_2"]
                instance["pos_2"] = int(row["pos_2"])
                instance["ref_2"] = row["ref_2"]
                instance["alt_2"] = row["alt_2"]

            errors = list(schema_validator.iter_errors(instance))
            assert not errors, (
                f"Row {idx+1} schema error: {[e.message for e in errors]}"
            )


@pytest.mark.parametrize("csv_path", CSV_TARGETS)
def test_submission_csv_exact_physical_line_count(csv_path: Path) -> None:
    """Verifies the CSV has 1 header + 1-10 data rows (the challenge accepts up to 10) and no ghost empty lines."""
    with open(csv_path, "rb") as f:
        content = f.read()

    lines = content.splitlines()
    assert 2 <= len(lines) <= 11, f"Expected 2-11 physical lines (header + 1-10 rows) in {csv_path.name}, got {len(lines)}"
    assert not content.endswith(b"\n\n"), f"Ghost double newline detected at end of {csv_path.name}"
    assert not content.endswith(b"\r\n\r\n"), f"Ghost double CRLF detected at end of {csv_path.name}"
    assert content.endswith(b"\n"), f"File {csv_path.name} must terminate with a newline"

