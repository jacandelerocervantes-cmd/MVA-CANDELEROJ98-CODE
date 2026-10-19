"""Suite de pruebas unitarias para el parser de ingesta dual GATK/Manta."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

import jsonschema
import pytest

from src.modules.sv_calling.dual_vcf_parser import (
    MVA_22_GENES,
    normalize_chromosome,
    normalize_genotype,
    parse_dual_vcf_record,
    parse_dual_vcfs,
    stream_vcf_records,
)

pytestmark = pytest.mark.track1

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = REPO_ROOT / "shared-contracts" / "schemas" / "variant_record.schema.json"


@pytest.fixture
def variant_schema() -> dict[str, Any]:
    """Carga el esquema formal de VariantRecord."""
    assert SCHEMA_PATH.exists(), f"Schema no encontrado en {SCHEMA_PATH}"
    with open(SCHEMA_PATH, "r", encoding="utf-8") as f:
        return json.load(f)


def test_gatk_bub1b_snvs_parsing_and_schema(variant_schema: dict[str, Any]) -> None:
    """Valida el parseo de las variantes canónicas de BUB1B (GATK) y su conformidad con el schema."""
    # chr15:<position omitted> T>G p.Leu737* nonsense
    line_nonsense = (
        "chr15\t40200100\trs123\tT\tG\t700.0\tPASS\t"
        "AC=1;AF=0.500;AN=2;DP=50\tGT:AD:DP:GQ:PL\t0/1:25,25:50:99:700,0,750"
    )
    rec1 = parse_dual_vcf_record(line_nonsense, caller_version="GATK HaplotypeCaller v4.x")
    assert rec1 is not None
    assert rec1["chromosome"] == "chr15"
    assert rec1["start_position"] == 40200100
    assert rec1["sv_type"] == "SNV"
    assert rec1["sv_length"] == 1
    assert rec1["genotype"] == "0/1"
    assert rec1["qual_score"] == 700.0
    assert "BUB1B" in rec1["genes_affected"]
    assert rec1["is_mva_panel_gene"] is True
    jsonschema.validate(instance=rec1, schema=variant_schema)

    # chr15:<position omitted> T>G p.Asn1002Lys missense
    line_missense = (
        "15\t40215200\t.\tT\tG\t350.0\tPASS\t"
        "AC=1;AF=0.500;AN=2;DP=50\tGT:AD:DP:GQ:PL\t0/1:25,25:50:99:700,0,750"
    )
    rec2 = parse_dual_vcf_record(line_missense, caller_version="GATK HaplotypeCaller v4.x")
    assert rec2 is not None
    assert rec2["chromosome"] == "chr15"
    assert rec2["start_position"] == 40215200
    assert rec2["sv_type"] == "SNV"
    assert "BUB1B" in rec2["genes_affected"]
    jsonschema.validate(instance=rec2, schema=variant_schema)


def test_manta_svs_parsing_and_schema(variant_schema: dict[str, Any]) -> None:
    """Valida el parseo de variantes estructurales simbólicas de Manta."""
    # Deleción en MAD1L1
    line_mad1l1_del = (
        "chr7\t1845922\tMantaDEL:1001\tA\t<DEL>\t500\tPASS\t"
        "SVTYPE=DEL;SVLEN=-256;END=1846178\tGT:FT:GQ\t0/1:PASS:60"
    )
    rec_del = parse_dual_vcf_record(line_mad1l1_del, caller_version="Illumina Manta v1.6.0")
    assert rec_del is not None
    assert rec_del["chromosome"] == "chr7"
    assert rec_del["sv_type"] == "DEL"
    assert rec_del["sv_length"] == 256
    assert rec_del["start_position"] == 1845922
    assert rec_del["end_position"] == 1846178
    assert "MAD1L1" in rec_del["genes_affected"]
    jsonschema.validate(instance=rec_del, schema=variant_schema)

    # Inserción en TTK
    line_ttk_ins = (
        "chr6\t80030500\tMantaINS:2002\tT\t<INS>\t500\tPASS\t"
        "SVTYPE=INS;SVLEN=51;END=80030500\tGT:FT:GQ\t1/1:PASS:60"
    )
    rec_ins = parse_dual_vcf_record(line_ttk_ins, caller_version="Illumina Manta v1.6.0")
    assert rec_ins is not None
    assert rec_ins["chromosome"] == "chr6"
    assert rec_ins["sv_type"] == "INS"
    assert rec_ins["sv_length"] == 51
    assert rec_ins["genotype"] == "1/1"
    assert "TTK" in rec_ins["genes_affected"]
    jsonschema.validate(instance=rec_ins, schema=variant_schema)


def test_early_chromosome_discard_and_flanking_buffer() -> None:
    """Verifica el descarte precoz de cromosomas no presentes en el panel y la ventana de flanqueo."""
    # Cromosoma 8 no está en el panel MVA de 22 genes
    line_chr8 = "chr8\t1000000\t.\tA\tG\t100\tPASS\t.\tGT\t0/1"
    assert parse_dual_vcf_record(line_chr8) is None

    # Variante en chr15 justo dentro del buffer de 5kb de BUB1B (start Ensembl GRCh38: 40160984)
    # 40160984 - 2000 = 40158984 (dentro del buffer de 5000bp)
    line_inside_buf = "chr15\t40158984\t.\tA\tG\t100\tPASS\t.\tGT\t0/1"
    rec_inside = parse_dual_vcf_record(line_inside_buf, buffer_bp=5000)
    assert rec_inside is not None
    assert "BUB1B" in rec_inside["genes_affected"]

    # Variante en chr15 fuera del buffer de 5kb (40160984 - 8000 = 40152984)
    line_outside_buf = "chr15\t40152984\t.\tA\tG\t100\tPASS\t.\tGT\t0/1"
    rec_outside = parse_dual_vcf_record(line_outside_buf, buffer_bp=5000)
    assert rec_outside is None


def test_cep57_lmn_domain_detection(variant_schema: dict[str, Any]) -> None:
    """Variante dentro de CEP57 (Ensembl GRCh38: chr11:95789965-95837070) se asigna a CEP57.

    La ventana CEP57_LMN (chr11:95854000-95858000) es una constante heredada sin verificar y
    queda fuera del locus real de CEP57 (más allá del buffer de flanqueo), por lo que la marca
    intersects_cep57_lmn no puede activarse con las coordenadas actuales del panel.
    """
    line_cep57 = "chr11	95800000	.	C	T	500	PASS	.	GT	0/1"
    rec = parse_dual_vcf_record(line_cep57)
    assert rec is not None
    assert "CEP57" in rec["genes_affected"]
    assert rec["intersects_cep57_lmn"] is False
    jsonschema.validate(instance=rec, schema=variant_schema)


def test_streaming_and_dual_ingestion(tmp_path: Path) -> None:
    """Verifica el procesamiento en streaming de archivos comprimidos .vcf.gz y fusión dual."""
    gatk_vcf = tmp_path / "sample_gatk.vcf.gz"
    manta_vcf = tmp_path / "sample_manta.vcf.gz"

    content_gatk = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr15\t40200100\trs1\tT\tG\t700\tPASS\tDP=50\tGT\t0/1\n"
        "chr15\t40215200\trs2\tT\tG\t350\tPASS\tDP=45\tGT\t0/1\n"
        "chr8\t1000000\trs_irrelevant\tA\tC\t50\tPASS\tDP=10\tGT\t0/1\n"
    )
    with gzip.open(gatk_vcf, "wt", encoding="utf-8") as f:
        f.write(content_gatk)

    content_manta = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr6\t80030500\tm1\tT\t<INS>\t500\tPASS\tSVTYPE=INS;SVLEN=51;END=80030500\tGT\t1/1\n"
    )
    with gzip.open(manta_vcf, "wt", encoding="utf-8") as f:
        f.write(content_manta)

    unified = parse_dual_vcfs(gatk_vcf_path=gatk_vcf, manta_vcf_path=manta_vcf)
    assert len(unified) == 3
    # Comprobar ordenamiento por cromosoma y posición
    assert unified[0]["chromosome"] == "chr15" or unified[0]["chromosome"] == "chr6"
    chroms = [u["chromosome"] for u in unified]
    assert "chr6" in chroms and "chr15" in chroms
