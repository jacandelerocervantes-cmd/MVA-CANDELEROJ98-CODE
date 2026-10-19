"""Convierte el mejor candidato genome-wide en las filas candidatas de la submission.

El gen y las variantes del par primario provienen del ranking calculado sobre el VCF; aquí no hay
ningún gen ni posición fija.
"""

from __future__ import annotations

import math
from typing import Any

from src.modules.sv_calling.compound_het_detector import (
    GENE_DISEASE_NOTES,
    NO_NOTE,
    CandidatePair,
    _quality_gate,
)
from src.modules.sv_calling.genome_wide_search import ABSENT_AF_FLOOR

_PHASE_NOTE = "phase undetermined (unphased short reads)"
# Convención del autor: un gen sin anotación HPO tiene fenotipo desconocido (no evidencia en contra) y puntúa
# este mínimo, para que su candidato pueda listarse y se ordene después de los de fenotipo conocido.
PHENOTYPE_UNKNOWN_FLOOR = 0.01


def pair_rarity_score(af_1: float | None, af_2: float | None, n_candidate_genes: int) -> tuple[float, float]:
    """(log10 de la rareza conjunta, valor de ordenamiento en (0, 1]) para un par de alelos raros.

    Rareza conjunta R = 1 / (4 f1 f2), el inverso de la probabilidad de portar por azar ambos alelos; no es un
    factor de Bayes contra una hipótesis causal: ignora la severidad de las variantes, la fase y el número de
    genes examinados (solo entra como prior 1 / nº de genes con variantes raras). El valor de ordenamiento
    R·p / (R·p + 1 - p) satura en 1.0 y no es una probabilidad calibrada de causalidad.
    """
    f1 = max(af_1, ABSENT_AF_FLOOR) if af_1 is not None else ABSENT_AF_FLOOR
    f2 = max(af_2, ABSENT_AF_FLOOR) if af_2 is not None else ABSENT_AF_FLOOR
    joint_rarity = 1.0 / (4.0 * f1 * f2)
    prior = 1.0 / max(n_candidate_genes, 2)
    ranking_value = (joint_rarity * prior) / (joint_rarity * prior + (1.0 - prior))
    return math.log10(joint_rarity), ranking_value


def _label(variant: dict[str, Any]) -> str:
    return variant.get("hgvsp") or variant["consequence"]


def _homozygous_primary(top: dict[str, Any], n_candidate_genes: int) -> list[CandidatePair]:
    variant = top["variants"][0]
    af = variant.get("af")
    f = max(af, ABSENT_AF_FLOOR) if af is not None else ABSENT_AF_FLOOR
    prior = 1.0 / max(n_candidate_genes, 2)
    joint_rarity = 1.0 / (f if top["type"] == "hemizygous" else f * f)
    ranking_value = (joint_rarity * prior) / (joint_rarity * prior + (1.0 - prior))
    evidence = GENE_DISEASE_NOTES.get(top["gene"], NO_NOTE)
    return [CandidatePair(
        gene=top["gene"], chrom_1=variant["chrom"], pos_1=variant["pos"], ref_1=variant["ref"], alt_1=variant["alt"],
        architecture="homozygous", finding_type="primary", base_score=round(max(0.01, min(1.0, ranking_value)), 2),
        notes=(f"{top['gene']} {top['type']} {variant['consequence']} {_label(variant)} {variant['chrom']}:{variant['pos']} "
               f"{variant['ref']}>{variant['alt']}, top-ranked recessive candidate genome-wide (rank {top['rank']}). "
               f"Gene-disease note: {evidence}."),
    )]


def build_primary_candidates(top: dict[str, Any], n_candidate_genes: int) -> list[CandidatePair]:
    """Fila primaria, orientación recíproca y componentes individuales del mejor candidato recesivo."""
    if top["type"] in ("homozygous", "hemizygous"):
        return _homozygous_primary(top, n_candidate_genes)
    a, b = sorted(top["variants"], key=lambda v: (v["impact"] != "HIGH", v["pos"]))
    gene = top["gene"]
    log_bf, ranking_value = pair_rarity_score(a.get("af"), b.get("af"), n_candidate_genes)
    if ranking_value >= 0.99:
        primary_score = round(min(1.0, ranking_value), 2)
        reciprocal_score = round(primary_score * 0.95, 2)
    else:
        primary_score = round(max(0.01, ranking_value), 2)
        reciprocal_score = round(primary_score * 0.50, 2)
    evidence = GENE_DISEASE_NOTES.get(gene, NO_NOTE)
    pair_text = f"{_label(a)} ({a['consequence']}) + {_label(b)} ({b['consequence']})"

    phenotype = top.get("phenotype_score", 0.0)

    def single(v: dict[str, Any]) -> CandidatePair:
        record = {"qual_score": v["qual"], "gq": v.get("gq"), "ref_allele": v["ref"], "alt_allele": v["alt"]}
        return CandidatePair(
            gene=gene, chrom_1=v["chrom"], pos_1=v["pos"], ref_1=v["ref"], alt_1=v["alt"],
            architecture="single_hit", finding_type="secondary",
            base_score=phenotype * _quality_gate(record),
            notes=(f"{gene} {v['consequence']} {_label(v)} {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']} "
                   f"(QUAL {v['qual']}, GQ {v.get('gq')}, GT {v['gt']}); single-variant component of the primary pair, listed separately."),
        )

    primary = CandidatePair(
        gene=gene, chrom_1=a["chrom"], pos_1=a["pos"], ref_1=a["ref"], alt_1=a["alt"],
        chrom_2=b["chrom"], pos_2=b["pos"], ref_2=b["ref"], alt_2=b["alt"],
        architecture="compound_heterozygous", finding_type="primary", base_score=primary_score,
        notes=(f"{gene} compound-heterozygous pair, top-ranked recessive candidate genome-wide "
               f"(rank {top['rank']}; log10 joint rarity {log_bf:.2f}; phenotype score {top.get('phenotype_score', 0.0):.3f}): {a['chrom']}:{a['pos']} {a['ref']}>{a['alt']} "
               f"{_label(a)} + {b['chrom']}:{b['pos']} {b['ref']}>{b['alt']} {_label(b)} "
               f"({_PHASE_NOTE}). Gene-disease note: {evidence}."),
    )
    reciprocal = CandidatePair(
        gene=gene, chrom_1=b["chrom"], pos_1=b["pos"], ref_1=b["ref"], alt_1=b["alt"],
        chrom_2=a["chrom"], pos_2=a["pos"], ref_2=a["ref"], alt_2=a["alt"],
        architecture="compound_heterozygous", finding_type="secondary", base_score=reciprocal_score,
        notes=f"{gene} compound-heterozygous pair listed in reciprocal variant order ({pair_text}): alternate representation of the primary pair.",
    )
    singles = [single(v) for v in (a, b)] if phenotype > 0 else []
    return [primary, reciprocal, *singles]


def build_secondary_candidate(candidate: dict[str, Any]) -> CandidatePair | None:
    """Fila secundaria para un candidato genome-wide que no es el primario.

    Puntaje = similitud fenotípica HPO del gen (mínimo PHENOTYPE_UNKNOWN_FLOOR si el gen no tiene anotación)
    x puerta de calidad del peor llamado, de modo que epcr es siempre > 0.
    """
    variants = candidate["variants"]
    gene = candidate["gene"]
    phenotype = candidate.get("phenotype_score", 0.0)
    gate = min(_quality_gate({"qual_score": v["qual"], "gq": v.get("gq"), "ref_allele": v["ref"], "alt_allele": v["alt"]}) for v in variants)
    first, second = (variants + [None])[:2]
    labels = " + ".join(f"{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']} {_label(v)} ({v['consequence']}, QUAL {v['qual']}, GQ {v.get('gq')}, GT {v['gt']})" for v in variants)
    return CandidatePair(
        gene=gene, chrom_1=first["chrom"], pos_1=first["pos"], ref_1=first["ref"], alt_1=first["alt"],
        chrom_2=second["chrom"] if second else None, pos_2=second["pos"] if second else None,
        ref_2=second["ref"] if second else None, alt_2=second["alt"] if second else None,
        architecture="compound_heterozygous" if second else "homozygous", finding_type="secondary",
        base_score=max(phenotype, PHENOTYPE_UNKNOWN_FLOOR) * gate,
        notes=(f"{gene} {candidate['type'].replace('_', '-')}, genome-wide rank {candidate['rank']} "
               f"(variant-evidence rank {candidate['variant_rank']}; HPO phenotype score {phenotype:.3f}): {labels}."),
    )
