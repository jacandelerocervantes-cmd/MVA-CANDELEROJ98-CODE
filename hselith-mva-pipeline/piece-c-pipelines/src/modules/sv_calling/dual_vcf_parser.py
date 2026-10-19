"""Parser de ingesta dual en streaming para variantes puntuales (GATK) y estructurales (Manta).

Cumple estrictamente con shared-contracts/schemas/variant_record.schema.json y límites de AGENTS.md.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class GeneRegion:
    """Región genómica canónica de un gen en GRCh38."""

    symbol: str
    chromosome: str
    start: int
    end: int
    is_mva_tier1: bool = False


# Panel curado de 22 genes MVA (Nivel 1 clínico + Nivel 2 STRING). Coordenadas = cuerpo del gen en GRCh38.p14, Ensembl release 116 (REST /lookup/symbol, consultado 2026-10-05; verificado con NCBI Gene). CCDC84 = CENATAC (nombre actual).
MVA_22_GENES: list[GeneRegion] = [
    # Nivel 1: Evidencia clínica directa de MVA
    GeneRegion("BUB1B", "chr15", 40160984, 40221137, is_mva_tier1=True),
    GeneRegion("CEP57", "chr11", 95789965, 95837070, is_mva_tier1=True),
    GeneRegion("TRIP13", "chr5", 892849, 919357, is_mva_tier1=True),
    GeneRegion("CCDC84", "chr11", 118998051, 119015811, is_mva_tier1=True),
    GeneRegion("SLF2", "chr10", 100912963, 100965134, is_mva_tier1=True),
    GeneRegion("SMC5", "chr9", 70258270, 70354874, is_mva_tier1=True),
    GeneRegion("MAD1L1", "chr7", 1815787, 2233243, is_mva_tier1=True),
    GeneRegion("MAD2L1BP", "chr6", 43629494, 43640960, is_mva_tier1=True),
    GeneRegion("CEP192", "chr18", 12991283, 13125053, is_mva_tier1=True),
    GeneRegion("CENPE", "chr4", 103105349, 103198463, is_mva_tier1=True),
    GeneRegion("CDC20", "chr1", 43358907, 43363212, is_mva_tier1=True),
    GeneRegion("BUB1", "chr2", 110635468, 110678098, is_mva_tier1=True),
    # Nivel 2: Vecinos directos de la red STRING (S >= 0.40)
    GeneRegion("MAD2L2", "chr1", 11657230, 11692605, is_mva_tier1=False),
    GeneRegion("SHLD3", "chr5", 65624733, 65631064, is_mva_tier1=False),
    GeneRegion("CEP63", "chr3", 134485699, 134587789, is_mva_tier1=False),
    GeneRegion("ZWINT", "chr10", 56352710, 56361392, is_mva_tier1=False),
    GeneRegion("SKA1", "chr18", 50374996, 50394176, is_mva_tier1=False),
    GeneRegion("CDCA5", "chr11", 65066299, 65084164, is_mva_tier1=False),
    GeneRegion("CDCA3", "chr12", 6844723, 6852353, is_mva_tier1=False),
    GeneRegion("FAM83D", "chr20", 38926312, 38953106, is_mva_tier1=False),
    GeneRegion("ANAPC15", "chr11", 72106370, 72113688, is_mva_tier1=False),
    GeneRegion("TTK", "chr6", 80003887, 80043580, is_mva_tier1=False),
]

PANEL_CHROMOSOMES: frozenset[str] = frozenset(g.chromosome for g in MVA_22_GENES)
CEP57_LMN_START = 95854000
CEP57_LMN_END = 95858000


def normalize_chromosome(chrom: str) -> str:
    """Asegura el prefijo canónico 'chr' en el cromosoma."""
    clean = chrom.strip()
    return clean if clean.startswith("chr") else f"chr{clean}"


def normalize_genotype(raw_gt: str | None) -> str:
    """Normaliza el genotipo a uno de los patrones aceptados por el esquema contractual."""
    if not raw_gt or raw_gt in (".", "./.", ".|."):
        return "./."
    gt = raw_gt.strip()
    if gt in ("0/1", "1/0", "0|1", "1|0"):
        return "0/1"
    if gt in ("1/1", "1|1"):
        return "1/1"
    if gt in ("0/0", "0|0"):
        return "0/0"
    return "./."


def parse_vcf_info(info_str: str) -> dict[str, str]:
    """Parsea el campo INFO de un registro VCF en un diccionario llave-valor."""
    info_dict: dict[str, str] = {}
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            info_dict[k] = v
        else:
            info_dict[item] = "true"
    return info_dict


def parse_dual_vcf_record(
    line: str,
    caller_version: str = "GATK HaplotypeCaller",
    buffer_bp: int = 5000,
) -> dict[str, Any] | None:
    """Parsea una línea VCF y evalúa intersección contra el panel de 22 genes con buffer.

    Descarta precozmente variantes fuera de cromosomas del panel para optimización de RAM.
    """
    if line.startswith("#") or not line.strip():
        return None

    parts = line.strip().split("\t")
    if len(parts) < 8:
        return None

    chrom = normalize_chromosome(parts[0])
    # Early discard por cromosoma: filtra el ~90% de cromosomas irrelevantes en streaming
    if chrom not in PANEL_CHROMOSOMES:
        return None

    try:
        pos = int(parts[1])
    except ValueError:
        return None

    var_id = parts[2]
    ref = parts[3]
    alt = parts[4]
    qual = float(parts[5]) if parts[5] != "." else None
    filter_status = parts[6]
    info_dict = parse_vcf_info(parts[7])

    # Inferencia de sv_type y coordenadas
    if "SVTYPE" in info_dict:
        sv_type = info_dict["SVTYPE"].upper()
    elif alt.startswith("<") and alt.endswith(">"):
        sv_type = alt[1:-1].split(":")[0].upper()
    elif len(ref) == 1 and len(alt) == 1:
        sv_type = "SNV"
    else:
        sv_type = "INDEL"

    end_pos = int(info_dict["END"]) if "END" in info_dict else pos + max(len(ref), 1) - 1

    if "SVLEN" in info_dict:
        try:
            sv_length = abs(int(info_dict["SVLEN"].split(",")[0]))
        except ValueError:
            sv_length = abs(end_pos - pos + 1)
    elif sv_type == "SNV":
        sv_length = 1
    elif sv_type == "INDEL":
        sv_length = abs(len(alt) - len(ref))
    else:
        sv_length = abs(end_pos - pos + 1)

    # Genotipo normalizado
    genotype = "./."
    if len(parts) >= 10:
        fmt_keys = parts[8].split(":")
        fmt_vals = parts[9].split(":")
        if "GT" in fmt_keys:
            gt_idx = fmt_keys.index("GT")
            if gt_idx < len(fmt_vals):
                genotype = normalize_genotype(fmt_vals[gt_idx])

    # Evaluación de intersección génica con buffer de flanqueo
    genes_affected: list[str] = []
    intersects_cep57_lmn = False

    for gene in MVA_22_GENES:
        if gene.chromosome == chrom:
            start_buf = max(1, gene.start - buffer_bp)
            end_buf = gene.end + buffer_bp
            if max(pos, start_buf) <= min(end_pos, end_buf):
                genes_affected.append(gene.symbol)
                if gene.symbol == "CEP57":
                    if max(pos, CEP57_LMN_START) <= min(end_pos, CEP57_LMN_END):
                        intersects_cep57_lmn = True

    if not genes_affected:
        return None

    canonical_id = var_id if var_id != "." else f"{sv_type}:{chrom}:{pos}_{ref}>{alt}"

    return {
        "variant_id": canonical_id,
        "caller": caller_version,
        "chromosome": chrom,
        "start_position": pos,
        "end_position": end_pos,
        "sv_type": sv_type,
        "sv_length": sv_length,
        "filter_status": filter_status,
        "genotype": genotype,
        "qual_score": qual,
        "ref_allele": ref if ref != "." else None,
        "alt_allele": alt if alt != "." else None,
        "genes_affected": genes_affected,
        "is_mva_panel_gene": True,
        "intersects_cep57_lmn": intersects_cep57_lmn,
    }


def stream_vcf_records(
    vcf_path: Path,
    caller_version: str = "GATK HaplotypeCaller",
    buffer_bp: int = 5000,
) -> Iterator[dict[str, Any]]:
    """Iterador en streaming sobre un archivo VCF o VCF.GZ que emite registros del panel."""
    is_gz = vcf_path.suffix == ".gz" or vcf_path.name.endswith(".vcf.gz")
    open_fn = gzip.open if is_gz else open
    mode = "rt" if is_gz else "r"

    with open_fn(vcf_path, mode, encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = parse_dual_vcf_record(line, caller_version=caller_version, buffer_bp=buffer_bp)
            if rec is not None:
                yield rec


def parse_dual_vcfs(
    gatk_vcf_path: Path | None = None,
    manta_vcf_path: Path | None = None,
    buffer_bp: int = 5000,
) -> list[dict[str, Any]]:
    """Ingesta y unifica variantes de GATK y Manta sobre el panel MVA en una sola lista."""
    unified: list[dict[str, Any]] = []

    if gatk_vcf_path and gatk_vcf_path.exists():
        for rec in stream_vcf_records(
            gatk_vcf_path,
            caller_version="Sentieon Haplotyper 202308.02 + GATK VariantFiltration",
            buffer_bp=buffer_bp,
        ):
            unified.append(rec)

    if manta_vcf_path and manta_vcf_path.exists():
        for rec in stream_vcf_records(
            manta_vcf_path,
            caller_version="Illumina Manta v1.6.0",
            buffer_bp=buffer_bp,
        ):
            unified.append(rec)

    # Ordenar por cromosoma y posición para procesamiento determinista
    unified.sort(key=lambda r: (r["chromosome"], r["start_position"]))
    return unified

