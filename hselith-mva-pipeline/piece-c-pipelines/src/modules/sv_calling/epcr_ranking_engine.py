"""Motor de ordenamiento EPCR (valor de ordenamiento, no probabilidad calibrada) y generador de archivos de submission oficiales.

Produce submissions/track1/{hf_username}_gatk-manta-fusion.csv conforme a variant_ranking.schema.json.
Garantiza EPCR estrictamente monótono decreciente en (0, 1] y aislamiento de secundarios.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import jsonschema

from src.modules.sv_calling.compound_het_detector import CandidatePair

REPO_ROOT = Path(__file__).resolve().parents[4]
PROJECT_ROOT = REPO_ROOT.parent
SCHEMAS_DIR = REPO_ROOT / "shared-contracts" / "schemas"
OUTPUT_DIR = PROJECT_ROOT / "submissions" / "track1"
SCHEMA_FILE = SCHEMAS_DIR / "variant_ranking.schema.json"

OFFICIAL_COLUMNS = [
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


def validate_submission_rows(rows: list[dict[str, Any]], schema_path: Path = SCHEMA_FILE) -> None:
    """Valida cada fila contra variant_ranking.schema.json con Draft202012Validator."""
    with open(schema_path, "r", encoding="utf-8") as f:
        schema = json.load(f)

    validator = jsonschema.Draft202012Validator(schema)
    for idx, row in enumerate(rows):
        clean_row = {k: v for k, v in row.items() if v is not None and v != ""}
        errors = list(validator.iter_errors(clean_row))
        if errors:
            raise ValueError(f"Fila {idx+1} violó el esquema contractual: {[e.message for e in errors]}")


def order_epcr_rows(candidates: list[CandidatePair]) -> list[dict[str, Any]]:
    """Ordena los candidatos y exige epcr estrictamente decreciente en (0, 1]; genera las filas de la submission."""
    # Verificar y forzar ordenamiento descendente
    sorted_candidates = sorted(candidates, key=lambda c: c.base_score, reverse=True)

    rows: list[dict[str, Any]] = []
    prev_epcr = float("inf")

    for idx, cand in enumerate(sorted_candidates[:10]):
        row = cand.to_submission_dict()
        epcr = float(row["epcr"])

        if not (0.0 < epcr <= 1.0):
            raise ValueError(f"EPCR fuera de rango (0, 1]: {epcr}")

        if epcr >= prev_epcr:
            raise ValueError(f"EPCR no estrictamente monótono decreciente en fila {idx+1}: {epcr} >= {prev_epcr}")

        prev_epcr = epcr
        rows.append(row)

    return rows


def write_submission_csv(rows: list[dict[str, Any]], output_path: Path) -> None:
    """Escribe las filas validadas en el CSV de submission en formato canónico."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=OFFICIAL_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            # Reemplazar None por string vacío para columnas CSV opcionales
            formatted_row = {col: (row.get(col) if row.get(col) is not None else "") for col in OFFICIAL_COLUMNS}
            writer.writerow(formatted_row)


def generate_submission_files(
    candidates: list[CandidatePair],
    output_dir: Path = OUTPUT_DIR,
    hf_username: str = "candeleroj98",
) -> Path:
    """Genera candeleroj98_gatk-manta-fusion.csv como archivo único oficial."""
    rows = order_epcr_rows(candidates)
    validate_submission_rows(rows)

    candelero_gatk_csv = output_dir / f"{hf_username}_gatk-manta-fusion.csv"
    write_submission_csv(rows, candelero_gatk_csv)

    return candelero_gatk_csv
