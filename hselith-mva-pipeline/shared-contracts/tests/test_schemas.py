"""Suite de pruebas para validar esquemas JSON contractuales contra Draft 2020-12."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import jsonschema
from jsonschema import Draft202012Validator
import pytest

SCHEMAS_DIR = Path(__file__).resolve().parents[1] / "schemas"


def load_all_schemas() -> list[tuple[str, dict[str, Any]]]:
    """Carga todos los esquemas JSON de shared-contracts/schemas/."""
    schema_files = list(SCHEMAS_DIR.glob("*.json"))
    schemas = []
    for sf in schema_files:
        with open(sf, "r", encoding="utf-8") as f:
            schemas.append((sf.name, json.load(f)))
    return schemas


@pytest.mark.parametrize("schema_name,schema_data", load_all_schemas())
def test_schema_valid_draft_2020_12(schema_name: str, schema_data: dict[str, Any]) -> None:
    """Verifica que cada esquema contractual sea sintácticamente válido según Draft 2020-12."""
    assert "$schema" in schema_data, f"{schema_name} debe declarar $schema"
    # Validar la estructura del esquema contra la meta-especificación Draft 2020-12
    Draft202012Validator.check_schema(schema_data)


def test_variant_record_schema_instance() -> None:
    """Valida un registro válido e inválido contra variant_record.schema.json."""
    with open(SCHEMAS_DIR / "variant_record.schema.json", "r", encoding="utf-8") as f:
        schema = json.load(f)

    valid_record = {
        "variant_id": "MantaDEL:101",
        "caller": "Manta v1.6.0",
        "chromosome": "chr11",
        "start_position": 95854100,
        "end_position": 95856100,
        "sv_type": "DEL",
        "sv_length": 2000,
        "filter_status": "PASS",
        "genotype": "0/1",
        "qual_score": 999.0,
        "ref_allele": "N",
        "alt_allele": "<DEL>",
        "genes_affected": ["CEP57"],
        "is_mva_panel_gene": True,
        "intersects_cep57_lmn": True,
    }
    jsonschema.validate(instance=valid_record, schema=schema)

    # Invalidez: falta campo requerido
    invalid_record = valid_record.copy()
    del invalid_record["genotype"]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(instance=invalid_record, schema=schema)


def test_hpo_profile_schema_instance() -> None:
    """Valida un registro válido e inválido contra hpo_profile.schema.json."""
    with open(SCHEMAS_DIR / "hpo_profile.schema.json", "r", encoding="utf-8") as f:
        schema = json.load(f)

    valid_profile = {
        "profile_id": "550e8400-e29b-41d4-a716-446655440000",
        "proband_id": "MVA_PROBAND_01",
        "hpo_terms": [
            {"hpo_id": "HP:0002859", "label": "Rhabdomyosarcoma", "is_negated": False},
            {"hpo_id": "HP:0000121", "label": "Nephrocalcinosis", "is_negated": False},
        ],
        "clinical_notes": "Paciente neonato con aneuploidías variegadas y falla medular.",
        "created_at": "2026-09-11T12:00:00Z",
    }
    jsonschema.validate(instance=valid_profile, schema=schema)


def test_dti_result_schema_instance() -> None:
    """Valida un resultado válido contra dti_result.schema.json."""
    with open(SCHEMAS_DIR / "dti_result.schema.json", "r", encoding="utf-8") as f:
        schema = json.load(f)

    valid_dti = {
        "dti_prediction_id": "a0eebc99-9c0b-4ef8-bb6d-6bb9bd380a11",
        "compound_name": "CFI-402257",
        "smiles": "CC1=CN=C(N=C1)N2CCN(CC2)C(=O)C3=C(N4C(=N3)C=CC=C4)NC5=CC(=C(C=C5)F)F",
        "target_name": "TTK",
        "target_uniprot_id": "P33981",
        "target_sequence": "MESEDLSGRELTIDSIMNKVRDIKNKFKNEDLTDELSLNKISADTTDNSGTVNQIMMMANN",
        "model_name": "MPNN_CNN_DAVIS",
        "predicted_pkd": 6.3153,
        "estimated_kd_nm": 483.81,
        "affinity_class": "moderate_affinity",
        "created_at": "2026-09-11T12:00:00Z",
    }
    jsonschema.validate(instance=valid_dti, schema=schema)


def test_stp_result_schema_instance() -> None:
    """Valida un resultado válido contra stp_result.schema.json."""
    with open(SCHEMAS_DIR / "stp_result.schema.json", "r", encoding="utf-8") as f:
        schema = json.load(f)

    valid_stp = {
        "stp_result_id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
        "module_name": "kinetochore_sac",
        "target_node": "TRIP13",
        "optimal_policy": {"state_0": "inhibit_0.5", "state_1": "nominal"},
        "value_function": {"state_0": 0.42, "state_1": 0.05},
        "convergence_iterations": 24,
        "bellman_residual": 0.00045,
        "perturbation_regime": "biallelic_truncation",
        "created_at": "2026-09-11T12:00:00Z",
    }
    jsonschema.validate(instance=valid_stp, schema=schema)

