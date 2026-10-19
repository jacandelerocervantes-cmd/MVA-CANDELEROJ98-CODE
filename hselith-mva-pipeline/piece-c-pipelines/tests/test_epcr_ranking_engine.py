"""Suite de pruebas unitarias para el motor de ordenamiento EPCR y generación de submission."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import pytest

from src.modules.sv_calling.compound_het_detector import CandidatePair
from src.modules.sv_calling.epcr_ranking_engine import (
    order_epcr_rows,
    generate_submission_files,
    validate_submission_rows,
    write_submission_csv,
)

pytestmark = pytest.mark.track1


@pytest.fixture
def sample_candidates() -> list[CandidatePair]:
    """Genera 10 candidatos sintéticos con scores decrecientes."""
    scores = [1.0, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.6, 0.5, 0.4]
    candidates = []
    for idx, s in enumerate(scores):
        is_primary = idx < 4
        candidates.append(
            CandidatePair(
                gene="BUB1B" if is_primary else "MAD1L1",
                chrom_1="chr15" if is_primary else "chr7",
                pos_1=40200100 + idx,
                ref_1="T",
                alt_1="G",
                base_score=s,
                finding_type="primary" if is_primary else "secondary",
                notes=f"Test note {idx+1}",
            )
        )
    return candidates


def test_order_epcr_rows_strictly_monotonic(sample_candidates: list[CandidatePair]) -> None:
    """Checks that the EPCR ranking value (not a calibrated probability) is strictly decreasing in (0, 1]."""
    rows = order_epcr_rows(sample_candidates)
    assert len(rows) == 10

    for i in range(len(rows) - 1):
        assert float(rows[i]["epcr"]) > float(rows[i + 1]["epcr"])

    # Probar que falla si hay duplicados o no es monótono
    bad_candidates = sample_candidates[:2] + [sample_candidates[1]]  # Score repetido
    with pytest.raises(ValueError, match="no estrictamente monótono decreciente"):
        order_epcr_rows(bad_candidates)


def test_generate_submission_files(tmp_path: Path, sample_candidates: list[CandidatePair]) -> None:
    """Valida la generación física del CSV oficial de submission."""
    fusion_csv = generate_submission_files(sample_candidates, output_dir=tmp_path)

    assert fusion_csv.exists()

    with open(fusion_csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    assert len(rows) == 10
    # Primeras 4 filas primarias (BUB1B)
    for i in range(4):
        assert rows[i]["finding_type"] == "primary"
        assert rows[i]["chrom_1"] == "chr15"

    # Filas 5 a 10 secundarias
    for i in range(4, 10):
        assert rows[i]["finding_type"] == "secondary"
        assert rows[i]["chrom_1"] == "chr7"


def test_submission_csv_uses_lf_line_endings(tmp_path: Path, sample_candidates: list[CandidatePair]) -> None:
    out = tmp_path / "x.csv"
    write_submission_csv(order_epcr_rows(sample_candidates), out)
    raw = out.read_bytes()
    assert b"\r" not in raw and raw.count(b"\n") == 11
