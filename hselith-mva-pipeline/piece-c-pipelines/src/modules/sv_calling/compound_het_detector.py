"""Hallazgos estructurales secundarios y representación de las filas candidatas de la submission.

Los SV de Manta se clasifican dentro del panel curado de 22 genes MVA (no hay frecuencia poblacional
de SV con la que establecer rareza en todo el genoma). El puntaje de cada fila secundaria es
calculado: similitud fenotípica HPO del gen (phenotype_scoring.py) x puerta de calidad de la llamada.
No hay genes, posiciones ni pesos escritos a mano.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class CandidatePair:
    """Representa una predicción candidata, ya sea par compound-het o variante individual."""

    gene: str
    chrom_1: str
    pos_1: int
    ref_1: str
    alt_1: str
    chrom_2: str | None = None
    pos_2: int | None = None
    ref_2: str | None = None
    alt_2: str | None = None
    architecture: str = "compound_heterozygous"  # compound_heterozygous, homozygous, single_hit
    finding_type: str = "primary"  # primary, secondary
    base_score: float = 1.0
    notes: str = ""

    def to_submission_dict(self, proband_id: str = "PROBAND01") -> dict[str, Any]:
        """Convierte a diccionario compatible con variant_ranking.schema.json."""
        d: dict[str, Any] = {
            "proband_id": proband_id,
            "chrom_1": self.chrom_1,
            "pos_1": self.pos_1,
            "ref_1": self.ref_1,
            "alt_1": self.alt_1,
            "chrom_2": self.chrom_2,
            "pos_2": self.pos_2,
            "ref_2": self.ref_2,
            "alt_2": self.alt_2,
            "epcr": round(self.base_score, 4),
            "finding_type": self.finding_type,
            "notes": self.notes,
        }
        return d


STRUCTURAL_SV_TYPES: frozenset[str] = frozenset({"DEL", "INS", "DUP", "INV"})
MIN_STRUCTURAL_SV_LENGTH_BP: int = 50
MAX_SECONDARY_STRUCTURAL_FINDINGS: int = 6
SECONDARY_SCORE_CEILING: float = 0.80
SECONDARY_SCORE_STEP_EPSILON: float = 1e-4

# Notas descriptivas gen-enfermedad (OMIM verificado en MedGen/NCBI; artículos verificados en PubMed).
# Son solo texto para las notas de la submission: NO intervienen en ningún puntaje.
GENE_DISEASE_NOTES: dict[str, str] = {
    "BUB1B": "MVA1, OMIM #257300 (Hanks 2004; Matsuura 2006; Garcia-Castillo 2008)",
    "CEP57": "MVA2, OMIM #614114 (Snape 2011)",
    "TRIP13": "MVA3, OMIM #617598 (Yost 2017)",
    "CENATAC": "MVA4, OMIM #620153 (de Wolf 2021)",
    "CCDC84": "alias of CENATAC: MVA4, OMIM #620153 (de Wolf 2021)",
    "SLF2": "Atelis syndrome 1 / MVA5, OMIM #620184 (Grange 2022)",
    "SMC5": "Atelis syndrome 2 / MVA6, OMIM #620185 (Grange 2022)",
    "BUB1": "MCPH30, OMIM #620183 (Carvalhal 2022)",
    "CENPE": "MCPH13, OMIM #616051 (Mirzaa 2014); microcephaly with CNS malformation",
    "MAD1L1": "MVA7, OMIM #620189 (Villarroya-Beltri 2022); limited case evidence",
    "MAD2L1BP": "biallelic MAD2L1BP mutation with mosaic aneuploidy (Abdel-Salam 2023, JCI Insight)",
    "CEP192": "MVA with tetraploidy (Guo 2024, HGG Advances)",
    "MAD2L2": "Fanconi anemia complementation group V, OMIM #617243; Fanconi phenotype, not MVA",
}
NO_NOTE = "no entry in the curated gene-disease notes"

_HOMOPOLYMER_RE = re.compile(r"(.)\1{9,}")
_TANDEM_UNIT_RE = re.compile(r"(.{1,4}?)\1{2,}")
_REPEAT_COVERAGE_FRACTION: float = 0.5


def _is_repetitive_context(ref: str, alt: str) -> bool:
    """Heurística de contexto de repetición (microsatélite/homopolímero).

    Illumina de lectura corta + Manta sobreestiman la confianza (QUAL) y generan llamadas espurias
    en este tipo de región (stutter, mapeo ambiguo). Se penaliza multiplicativamente, no se descarta.
    Usa la cobertura acumulada de todas las repeticiones cortas encontradas, porque los STR reales en
    datos Illumina vienen fragmentados por ruido de ensamblaje en islas cortas, no en un solo tramo.
    """
    seq = alt if len(alt) >= len(ref) else ref
    if _HOMOPOLYMER_RE.search(seq):
        return True
    covered = sum(len(m.group(0)) for m in _TANDEM_UNIT_RE.finditer(seq))
    return covered >= _REPEAT_COVERAGE_FRACTION * len(seq)


def _quality_gate(record: dict[str, Any]) -> float:
    """Puerta de calidad (0,1]: confianza del caller, nunca importancia biológica.

    SNV/indel de GATK: GQ/99 (el GQ del genotipo, techo 99; el QUAL de GATK satura y no discrimina).
    SV de Manta: QUAL normalizado a su techo empírico (999.0). En ambos casos x0.35 si la llamada cae
    en contexto repetitivo. Es un multiplicador, no un componente aditivo del puntaje fenotípico.
    """
    if record.get("gq") is not None:
        qual_norm = min(1.0, float(record["gq"]) / 99.0)
    else:
        qual_norm = min(1.0, float(record.get("qual_score") or 0.0) / 999.0)
    if _is_repetitive_context(str(record.get("ref_allele") or ""), str(record.get("alt_allele") or "")):
        qual_norm *= 0.35
    return qual_norm


def secondary_score(gene: str, phenotype_scores: dict[str, float], *records: dict[str, Any]) -> float:
    """Puntaje de una fila secundaria = similitud fenotípica HPO del gen x puerta de calidad.

    Para un par, la credibilidad la fija el eslabón más débil (mínimo de las puertas de calidad):
    la hipótesis bialélica requiere que ambos hallazgos sean reales.
    """
    return phenotype_scores.get(gene, 0.0) * min(_quality_gate(r) for r in records)


def _enforce_strictly_descending(
    candidates: list[CandidatePair], ceiling: float = SECONDARY_SCORE_CEILING
) -> list[CandidatePair]:
    """Re-emite la lista con base_score estrictamente descendente, preservando el orden real.

    Los puntajes crudos pueden empatar (Manta reporta QUAL 999.0 en muchos eventos PASS); el motor de
    ranking exige monotonía estricta en el CSV final. Solo se resuelven empates con un épsilon mínimo.
    """
    ordered = sorted(candidates, key=lambda c: c.base_score, reverse=True)
    resolved: list[CandidatePair] = []
    prev = ceiling + SECONDARY_SCORE_STEP_EPSILON
    for cand in ordered:
        score = min(cand.base_score, prev - SECONDARY_SCORE_STEP_EPSILON)
        prev = score
        resolved.append(
            CandidatePair(
                gene=cand.gene, chrom_1=cand.chrom_1, pos_1=cand.pos_1, ref_1=cand.ref_1, alt_1=cand.alt_1,
                chrom_2=cand.chrom_2, pos_2=cand.pos_2, ref_2=cand.ref_2, alt_2=cand.alt_2,
                architecture=cand.architecture, finding_type=cand.finding_type,
                base_score=round(score, 4), notes=cand.notes,
            )
        )
    return resolved


def _repeat_warning(*records: dict[str, Any]) -> str:
    if any(_is_repetitive_context(r["ref_allele"], r["alt_allele"]) for r in records):
        return (" WARNING: falls in repetitive/microsatellite context, elevated artifact risk, "
                "no orthogonal validation (long-read/PCR) performed.")
    return ""


NON_CODING_SV_FACTOR = 0.5  # convención declarada: un SV sin solape con la CDS MANE vale la mitad (no es una probabilidad estimada)


def _coding_factor(coding_overlap: Callable[[str, int, int], bool] | None, *records: dict[str, Any]) -> float:
    if coding_overlap is None:
        return 1.0
    hit = any(coding_overlap(r["chromosome"], r["start_position"], r.get("end_position") or r["start_position"]) for r in records)
    return 1.0 if hit else NON_CODING_SV_FACTOR


def _coding_note(coding_overlap: Callable[[str, int, int], bool] | None, *records: dict[str, Any]) -> str:
    if coding_overlap is None:
        return ""
    overlaps = any(coding_overlap(r["chromosome"], r["start_position"], r.get("end_position") or r["start_position"]) for r in records)
    return f" Overlaps MANE Select coding sequence: {'yes' if overlaps else 'no'}."


def detect_secondary_structural_architectures(
    records: list[dict[str, Any]],
    phenotype_scores: dict[str, float],
    coding_overlap: Callable[[str, int, int], bool] | None = None,
) -> list[CandidatePair]:
    """Clasifica alteraciones estructurales secundarias desde registros Manta reales.

    Filtra por tipo/longitud de SV (con secuencia expandida, no símbolos <DEL>/<INS>), agrupa por gen
    del panel, empareja heterocigotos en compound-het cuando hay >=2 por gen y usa ref/alt/QUAL
    literales del registro. Un panel sin hallazgos produce una lista vacía, no datos inventados. Las
    filas con puntaje fenotípico 0 (gen sin anotación HPO) se omiten porque epcr debe ser > 0.
    """
    sv_records = [
        r for r in records
        if r.get("sv_type") in STRUCTURAL_SV_TYPES
        and int(r.get("sv_length") or 0) >= MIN_STRUCTURAL_SV_LENGTH_BP
        and r.get("ref_allele") and r.get("alt_allele")
        and r.get("filter_status", "PASS") == "PASS"
        and not str(r.get("alt_allele", "")).startswith("<")
    ]
    by_gene: dict[str, list[dict[str, Any]]] = {}
    for record in sv_records:
        for gene in record.get("genes_affected", []):
            by_gene.setdefault(gene, []).append(record)

    secondaries: list[CandidatePair] = []
    for gene, gene_records in by_gene.items():
        gene_records = sorted(gene_records, key=lambda r: r["start_position"])
        heterozygous = [r for r in gene_records if r.get("genotype") == "0/1"]
        note = GENE_DISEASE_NOTES.get(gene, NO_NOTE)
        pheno = phenotype_scores.get(gene, 0.0)

        if len(heterozygous) >= 2:
            r1, r2 = heterozygous[0], heterozygous[1]
            score = secondary_score(gene, phenotype_scores, r1, r2) * _coding_factor(coding_overlap, r1, r2)
            if score > 0:
                secondaries.append(CandidatePair(
                    gene=gene, chrom_1=r1["chromosome"], pos_1=r1["start_position"], ref_1=r1["ref_allele"], alt_1=r1["alt_allele"],
                    chrom_2=r2["chromosome"], pos_2=r2["start_position"], ref_2=r2["ref_allele"], alt_2=r2["alt_allele"],
                    architecture="compound_heterozygous", finding_type="secondary", base_score=score,
                    notes=(f"{gene} two heterozygous SV calls in the same gene, phase undetermined (compound-heterozygous only if in trans) ({r1['sv_type']} {r1['sv_length']}bp + "
                           f"{r2['sv_type']} {r2['sv_length']}bp; size per Manta SVLEN), real Manta calls. "
                           f"HPO phenotype score {pheno:.3f}. Gene-disease note: {note}."
                           + _coding_note(coding_overlap, r1, r2) + _repeat_warning(r1, r2)),
                ))
            gene_records = [r for r in gene_records if r is not r1 and r is not r2]

        for record in gene_records:
            genotype = record.get("genotype", "./.")
            architecture = "homozygous" if genotype == "1/1" else "single_hit"
            score = secondary_score(gene, phenotype_scores, record) * _coding_factor(coding_overlap, record)
            if score <= 0:
                continue
            secondaries.append(CandidatePair(
                gene=gene, chrom_1=record["chromosome"], pos_1=record["start_position"], ref_1=record["ref_allele"], alt_1=record["alt_allele"],
                architecture=architecture, finding_type="secondary", base_score=score,
                notes=(f"{gene} {architecture} {record['sv_type']} {record['sv_length']}bp "
                       f"(QUAL {record.get('qual_score')}, GT {genotype}; size per Manta SVLEN), real Manta call. "
                       f"HPO phenotype score {pheno:.3f}. Gene-disease note: {note}."
                       + _coding_note(coding_overlap, record) + _repeat_warning(record)),
            ))

    secondaries.sort(key=lambda c: c.base_score, reverse=True)
    return _enforce_strictly_descending(secondaries[:MAX_SECONDARY_STRUCTURAL_FINDINGS])
