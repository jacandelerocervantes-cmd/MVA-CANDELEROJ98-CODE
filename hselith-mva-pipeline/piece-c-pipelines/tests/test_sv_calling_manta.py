"""Suite de pruebas automatizadas para SV Calling con Illumina Manta v1.6.0."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import jsonschema
import pytest

from src.modules.sv_calling.call_sv import execute_sv_calling_pipeline
from src.modules.sv_calling.manta_config import MantaConfig
from src.modules.sv_calling.manta_runner import MantaRunner
from src.modules.sv_calling.sv_vcf_parser import (
    filter_mva_panel_variants,
    normalize_chromosome,
    normalize_genotype,
    parse_vcf_record,
    parse_vcf_records,
)

# Cargar el esquema contractual de variant_record
REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = REPO_ROOT / "shared-contracts" / "schemas" / "variant_record.schema.json"


@pytest.fixture
def variant_schema() -> dict[str, Any]:
    with open(SCHEMA_PATH, "r", encoding="utf-8") as f:
        return json.load(f)


def test_manta_config_command_generation(tmp_path: Path) -> None:
    """Verifica que MantaConfig genere comandos CLI con los flags correctos."""
    bam = tmp_path / "patient.bam"
    fasta = tmp_path / "hg38.fa"
    out = tmp_path / "manta_run"

    cfg = MantaConfig(
        bam_path=bam,
        reference_fasta=fasta,
        output_dir=out,
        threads=16,
        memory_gb=64,
        generate_evidence_bam=False,
    )

    config_cmd = cfg.build_config_command()
    assert any(f"--bam={str(bam)}" in arg for arg in config_cmd)
    assert any(f"--referenceFasta={str(fasta)}" in arg for arg in config_cmd)
    assert "--generateEvidenceBam" not in config_cmd
    cfg_with_evidence = MantaConfig(
        bam_path=bam, reference_fasta=fasta, output_dir=out, generate_evidence_bam=True
    )
    assert "--generateEvidenceBam" in cfg_with_evidence.build_config_command()

    run_cmd = cfg.build_run_command()
    assert "-m" in run_cmd and "local" in run_cmd
    assert "-j" in run_cmd and "16" in run_cmd
    assert "-g" in run_cmd and "64" in run_cmd


def test_manta_config_validation_errors(tmp_path: Path) -> None:
    """Verifica errores de validación en hilos y memoria inválidos."""
    with pytest.raises(ValueError, match="threads debe ser >= 1"):
        MantaConfig(Path("a.bam"), Path("r.fa"), threads=0).validate(check_files_exist=False)

    with pytest.raises(ValueError, match="memory_gb debe ser >= 4"):
        MantaConfig(Path("a.bam"), Path("r.fa"), memory_gb=2).validate(check_files_exist=False)


def test_vcf_parsing_and_mva_panel_mapping(variant_schema: dict[str, Any]) -> None:
    """Prueba el parseo de líneas VCF sintéticas y la clasificación hacia el panel MVA."""
    # 1. Deleción dentro del cuerpo génico de CEP57 en GRCh38 (chr11:95789965-95837070). La ventana LMN heredada
    #    (chr11:95854000-95858000) cae fuera del locus real, por eso intersects_cep57_lmn no puede activarse.
    line_cep57_lmn = (
        "chr11\t95800000\tMantaDEL:101\tN\t<DEL>\t999\tPASS\t"
        "SVTYPE=DEL;SVLEN=-2000;END=95802000\tGT:FT:GQ\t0/1:PASS:99"
    )
    rec1 = parse_vcf_record(line_cep57_lmn)
    assert rec1 is not None
    assert rec1["chromosome"] == "chr11"
    assert rec1["sv_type"] == "DEL"
    assert rec1["sv_length"] == 2000
    assert rec1["is_mva_panel_gene"] is True
    assert "CEP57" in rec1["genes_affected"]
    assert rec1["intersects_cep57_lmn"] is False
    assert rec1["genotype"] == "0/1"
    jsonschema.validate(instance=rec1, schema=variant_schema)

    # 2. Duplicación en TRIP13 (chr5:892849-919357)
    line_trip13 = (
        "5\t890000\tMantaDUP:202\tA\t<DUP:TANDEM>\t450\tPASS\t"
        "SVTYPE=DUP;SVLEN=5000;END=895000\tGT:FT\t1/1:PASS"
    )
    rec2 = parse_vcf_record(line_trip13)
    assert rec2 is not None
    assert rec2["chromosome"] == "chr5"
    assert rec2["sv_type"] == "DUP"
    assert rec2["is_mva_panel_gene"] is True
    assert "TRIP13" in rec2["genes_affected"]
    assert rec2["intersects_cep57_lmn"] is False
    assert rec2["genotype"] == "1/1"
    jsonschema.validate(instance=rec2, schema=variant_schema)

    # 3. Inversión fuera de genes panel (chr8)
    line_non_panel = (
        "chr8\t10000000\tMantaINV:303\tN\t<INV>\t120\tMinMQ\t"
        "SVTYPE=INV;END=10015000\tGT\t0/1"
    )
    rec3 = parse_vcf_record(line_non_panel)
    assert rec3 is not None
    assert rec3["is_mva_panel_gene"] is False
    assert len(rec3["genes_affected"]) == 0
    assert rec3["filter_status"] == "MinMQ"
    jsonschema.validate(instance=rec3, schema=variant_schema)


def test_parse_vcf_file_gz_and_filter(tmp_path: Path, variant_schema: dict[str, Any]) -> None:
    """Verifica lectura de archivo comprimido .vcf.gz y filtro de panel."""
    vcf_gz = tmp_path / "test_diploidSV.vcf.gz"
    vcf_content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr15\t40200000\tMantaDEL:1\tG\t<DEL>\t800\tPASS\tSVTYPE=DEL;SVLEN=-500;END=40200500\tGT\t0/1\n"
        "chr1\t1000000\tMantaDEL:2\tA\t<DEL>\t30\tLowQual\tSVTYPE=DEL;SVLEN=-100;END=1000100\tGT\t./.\n"
    )
    with gzip.open(vcf_gz, "wt", encoding="utf-8") as f:
        f.write(vcf_content)

    records = parse_vcf_records(vcf_gz)
    assert len(records) == 2

    # Ambos deben validar contra el esquema contractual
    for r in records:
        jsonschema.validate(instance=r, schema=variant_schema)

    panel_only = filter_mva_panel_variants(records)
    assert len(panel_only) == 1
    assert "BUB1B" in panel_only[0]["genes_affected"]


@patch("subprocess.run")
def test_manta_runner_execution_flow(mock_subproc: MagicMock, tmp_path: Path) -> None:
    """Simula la ejecución completa de MantaRunner con mocks."""
    out_dir = tmp_path / "manta_out"
    workflow_py = out_dir / "runWorkflow.py"
    variants_dir = out_dir / "results" / "variants"
    variants_dir.mkdir(parents=True, exist_ok=True)
    diploid_vcf = variants_dir / "diploidSV.vcf.gz"
    diploid_vcf.touch()

    # Callback para simular creación de runWorkflow.py tras configManta
    def fake_subprocess_run(cmd: list[str], **kwargs: Any) -> MagicMock:
        if "configManta.py" in cmd[1]:
            workflow_py.touch()
        res = MagicMock()
        res.returncode = 0
        res.stdout = "OK"
        res.stderr = ""
        return res

    mock_subproc.side_effect = fake_subprocess_run

    cfg = MantaConfig(
        bam_path=tmp_path / "dummy.bam",
        reference_fasta=tmp_path / "dummy.fa",
        output_dir=out_dir,
    )

    runner = MantaRunner(cfg)
    exec_result = runner.run(validate_inputs=False)

    assert exec_result.success is True
    assert exec_result.diploid_sv_vcf == diploid_vcf
    assert mock_subproc.call_count == 2


@patch("src.modules.sv_calling.call_sv.MantaRunner")
def test_execute_sv_calling_pipeline(mock_runner_cls: MagicMock, tmp_path: Path) -> None:
    """Valida el orquestador principal execute_sv_calling_pipeline."""
    out_dir = tmp_path / "pipe_out"
    variants_dir = out_dir / "results" / "variants"
    variants_dir.mkdir(parents=True, exist_ok=True)
    diploid_vcf = variants_dir / "diploidSV.vcf.gz"

    # Generar VCF sintético con 1 variante en CEP63 y 1 variante control
    vcf_content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr3\t134500000\tMantaDEL:1\tT\t<DEL>\t600\tPASS\tSVTYPE=DEL;SVLEN=-1000;END=134501000\tGT\t0/1\n"
        "chr8\t20000000\tMantaDUP:2\tC\t<DUP:TANDEM>\t500\tPASS\tSVTYPE=DUP;SVLEN=500;END=20000500\tGT\t0/1\n"
    )
    with gzip.open(diploid_vcf, "wt", encoding="utf-8") as f:
        f.write(vcf_content)

    mock_runner_inst = MagicMock()
    mock_runner_inst.run.return_value = MagicMock(
        success=True,
        diploid_sv_vcf=diploid_vcf,
        candidate_sv_vcf=None,
        duration_seconds=1.23,
        error_message=None,
    )
    mock_runner_cls.return_value = mock_runner_inst

    res = execute_sv_calling_pipeline(
        bam_path=tmp_path / "mock.bam",
        reference_fasta=tmp_path / "mock.fa",
        output_dir=out_dir,
        validate_inputs=False,
    )

    assert res["total_sv_count"] == 2
    assert res["pass_sv_count"] == 2
    assert res["panel_sv_count"] == 1
    assert Path(str(res["panel_json_path"])).is_file()
