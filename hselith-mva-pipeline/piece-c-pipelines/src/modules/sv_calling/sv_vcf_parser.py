"""Parser de VCF emitido por Illumina Manta y clasificador del panel de 22 genes MVA."""

from __future__ import annotations

import gzip
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

# Región específica del dominio LMN (residuos 257-264) en CEP57 (GRCh38 chr11)
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
    if gt in ("0/1", "1/0"):
        return "0/1"
    if gt in ("1/1", "1|1"):
        return "1/1"
    if gt in ("0/0", "0|0"):
        return "0/0"
    if gt in ("0|1", "1|0"):
        return "0|1"
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


def parse_vcf_record(line: str, caller_version: str = "Manta v1.6.0") -> dict[str, Any] | None:
    """Parsea y valida una línea de variante VCF conforme a variant_record.schema.json."""
    if line.startswith("#") or not line.strip():
        return None

    parts = line.strip().split("\t")
    if len(parts) < 8:
        return None

    chrom = normalize_chromosome(parts[0])
    try:
        pos = int(parts[1])
    except ValueError:
        return None

    var_id = parts[2]
    ref = parts[3]
    alt = parts[4]
    try:
        qual = float(parts[5]) if parts[5] != "." else None
    except ValueError:
        qual = None
    filter_status = parts[6]
    info = parse_vcf_info(parts[7])

    # Inferencia de sv_type
    sv_type = info.get("SVTYPE")
    if not sv_type:
        if "<DEL" in alt:
            sv_type = "DEL"
        elif "<DUP" in alt:
            sv_type = "DUP"
        elif "<INV" in alt:
            sv_type = "INV"
        elif "<INS" in alt:
            sv_type = "INS"
        elif "[" in alt or "]" in alt:
            sv_type = "BND"
        else:
            sv_type = "DEL" if len(ref) > len(alt) else ("INS" if len(alt) > len(ref) else "SNV")

    # Inferencia de sv_length
    sv_len: int | None = None
    if "SVLEN" in info:
        try:
            sv_len = abs(int(info["SVLEN"]))
        except ValueError:
            sv_len = None

    # Inferencia de end_position
    end_pos = pos
    if "END" in info:
        try:
            end_pos = int(info["END"])
        except ValueError:
            end_pos = pos + (sv_len if sv_len else 1)
    elif sv_len:
        end_pos = pos + sv_len
    else:
        end_pos = pos + max(len(ref), 1)

    # Extracción de Genotipo (GT)
    genotype = "./."
    if len(parts) >= 10:
        fmt_keys = parts[8].split(":")
        sample_vals = parts[9].split(":")
        if "GT" in fmt_keys:
            gt_idx = fmt_keys.index("GT")
            if gt_idx < len(sample_vals):
                genotype = normalize_genotype(sample_vals[gt_idx])

    # Cruce con genes diana
    genes_affected: list[str] = []
    intersects_cep57_lmn = False

    for gene in MVA_22_GENES:
        if gene.chromosome == chrom:
            if max(pos, gene.start) <= min(end_pos, gene.end):
                genes_affected.append(gene.symbol)
                if gene.symbol == "CEP57":
                    if max(pos, CEP57_LMN_START) <= min(end_pos, CEP57_LMN_END):
                        intersects_cep57_lmn = True

    is_mva_panel = len(genes_affected) > 0

    return {
        "variant_id": var_id if var_id != "." else f"Manta{sv_type}:{chrom}:{pos}-{end_pos}",
        "caller": caller_version,
        "chromosome": chrom,
        "start_position": pos,
        "end_position": end_pos,
        "sv_type": sv_type,
        "sv_length": sv_len,
        "filter_status": filter_status,
        "genotype": genotype,
        "qual_score": qual,
        "ref_allele": ref if ref != "." else None,
        "alt_allele": alt if alt != "." else None,
        "genes_affected": genes_affected,
        "is_mva_panel_gene": is_mva_panel,
        "intersects_cep57_lmn": intersects_cep57_lmn,
    }


def parse_vcf_records(vcf_path: Path, caller_version: str = "Manta v1.6.0") -> list[dict[str, Any]]:
    """Lee y parsea un archivo VCF (.vcf o .vcf.gz) retornando registros normalizados."""
    records: list[dict[str, Any]] = []
    is_gz = vcf_path.suffix == ".gz" or vcf_path.name.endswith(".vcf.gz")

    open_fn = gzip.open if is_gz else open
    mode = "rt" if is_gz else "r"

    with open_fn(vcf_path, mode, encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = parse_vcf_record(line, caller_version=caller_version)
            if rec is not None:
                records.append(rec)

    return records


def filter_mva_panel_variants(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Filtra y extrae únicamente aquellas variantes que intersectan los 22 genes del panel MVA."""
    return [r for r in records if r.get("is_mva_panel_gene", False)]

