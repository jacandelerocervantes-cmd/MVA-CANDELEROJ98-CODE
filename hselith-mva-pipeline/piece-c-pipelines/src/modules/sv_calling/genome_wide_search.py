"""Búsqueda genome-wide de candidatos recesivos a partir del VCF de SNV/indel.

Flujo (todo calculado, ninguna lista de genes ni de variantes de entrada):
1. Modelos de transcrito MANE Select (GTF) -> regiones codificantes y sitios de splicing.
2. Escaneo del VCF: variantes PASS dentro de esas regiones, con genotipo y marca dbSNP.
3. Consecuencia calculada por variante (consequence_calculator).
4. Filtros de calidad: rareza establecida (gnomAD v4), descarte de pares a menos de 150 pb y de
   familias génicas hipervariables. Se añadieron tras una primera pasada sin ellos (ver el reporte, 4b).
5. Candidatos por gen: compuesto heterocigoto u homocigoto; ranking por (nº de variantes HIGH,
   rareza conjunta).
Las funciones de red (secuencias CDS y frecuencias gnomAD) están en genome_wide_io.py.
"""

from __future__ import annotations

import bisect
import collections
import gzip
import math
import re
from pathlib import Path
from typing import Any

from src.modules.sv_calling.consequence_calculator import classify_allele

RARE_AF_THRESHOLD: float = 1e-3
ABSENT_AF_FLOOR: float = 1e-6
CIS_PROXIMITY_BP: int = 150
HYPERVARIABLE_GENES = re.compile(r"^(HLA-|KIR\d|OR\d|MUC\d|KRTAP|PRB\d|TAS2R)")
# Regiones pseudoautosómicas de chrX en GRCh38 (inclusive); fuera de ellas un varón es hemicigoto.
PAR_X_GRCH38 = ((10_001, 2_781_479), (155_701_383, 156_030_895))
HEMIZYGOUS_MIN_X_CALLS = 20
HEMIZYGOUS_MAX_HET_FRACTION = 0.25  # varón: 49 de 295 (17 %) en este VCF; una mujer diploide ronda 35-50 %


def is_x_nonpar(chrom: str, pos: int) -> bool:
    return chrom == "chrX" and not any(lo <= pos <= hi for lo, hi in PAR_X_GRCH38)


def infer_hemizygous_x(hits: list[dict[str, Any]]) -> bool:
    """True si el patrón de genotipos en chrX fuera de las PAR indica un probando con un solo cromosoma X.

    Un varón casi no tiene llamadas heterocigotas allí (las que hay son errores de mapeo); se exige un mínimo
    de llamadas y una fracción heterocigota baja. Se deduce de los datos, no se declara.
    """
    calls = [h["gt"].replace("|", "/").split("/") for h in hits if is_x_nonpar(h["chrom"], h["pos"])]
    calls = [c for c in calls if len(c) == 2 and all(x.isdigit() for x in c)]
    if len(calls) < HEMIZYGOUS_MIN_X_CALLS:
        return False
    het = sum(1 for c in calls if c[0] != c[1])
    return het / len(calls) <= HEMIZYGOUS_MAX_HET_FRACTION
_ATTR = re.compile(r'(\w+) "([^"]*)"')


def load_mane_transcripts(gtf_path: Path) -> dict[str, dict[str, Any]]:
    """Transcritos MANE Select codificantes: hebra, exones, CDS y codón de parada."""
    transcripts: dict[str, dict[str, Any]] = {}
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt", encoding="utf-8") as handle:  # type: ignore[operator]
        for line in handle:
            if line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] not in ("exon", "CDS", "stop_codon") or 'tag "MANE_Select"' not in fields[8]:
                continue
            attrs = dict(_ATTR.findall(fields[8]))
            entry = transcripts.setdefault(
                attrs["transcript_id"],
                {"chrom": fields[0], "strand": fields[6], "gene": attrs.get("gene_name"), "exons": [], "cds": [], "stop": []},
            )
            start, end = int(fields[3]), int(fields[4])
            if fields[2] == "exon":
                entry["exons"].append((start, end))
            elif fields[2] == "CDS":
                entry["cds"].append((start, end, int(fields[7])))
            else:
                entry["stop"].append((start, end))
    return {tid: t for tid, t in transcripts.items() if t["cds"]}


def _regions_index(transcripts: dict[str, dict[str, Any]]) -> tuple[dict[str, list[tuple[int, int, str]]], dict[str, list[int]], dict[str, int]]:
    regions: dict[str, list[tuple[int, int, str]]] = collections.defaultdict(list)
    for tid, t in transcripts.items():
        for start, end, _ in t["cds"]:
            regions[t["chrom"]].append((start, end, tid))
        for start, end in t["stop"]:
            regions[t["chrom"]].append((start, end, tid))
        exons = sorted(t["exons"])
        for i in range(len(exons) - 1):
            regions[t["chrom"]].append((exons[i][1] + 1, exons[i][1] + 2, tid))
            regions[t["chrom"]].append((exons[i + 1][0] - 2, exons[i + 1][0] - 1, tid))
    for chrom in regions:
        regions[chrom].sort()
    starts = {c: [r[0] for r in v] for c, v in regions.items()}
    max_len = {c: max(r[1] - r[0] for r in v) for c, v in regions.items()}
    return regions, starts, max_len


def coding_overlap_checker(transcripts: dict[str, dict[str, Any]]):
    """Devuelve f(chrom, start, end) -> bool: ¿el intervalo toca CDS (o codón de parada) de algún MANE Select?"""
    by_chrom: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    for t in transcripts.values():
        for start, end, _ in t["cds"]:
            by_chrom[t["chrom"]].append((start, end))
        for start, end in t["stop"]:
            by_chrom[t["chrom"]].append((start, end))
    merged: dict[str, list[tuple[int, int]]] = {}
    for chrom, intervals in by_chrom.items():
        out: list[tuple[int, int]] = []
        for start, end in sorted(intervals):
            if out and start <= out[-1][1] + 1:
                out[-1] = (out[-1][0], max(out[-1][1], end))
            else:
                out.append((start, end))
        merged[chrom] = out
    starts = {c: [i[0] for i in v] for c, v in merged.items()}

    def overlaps(chrom: str, start: int, end: int) -> bool:
        intervals = merged.get(chrom)
        if not intervals:
            return False
        j = bisect.bisect_right(starts[chrom], end) - 1
        return j >= 0 and intervals[j][1] >= start

    return overlaps


def _format_gq(format_field: str, sample_field: str) -> int | None:
    """GQ (calidad del genotipo, techo 99 en GATK) del primer campo de muestra, o None si no existe."""
    keys = format_field.split(":")
    values = sample_field.strip().split(":")
    if "GQ" in keys and keys.index("GQ") < len(values) and values[keys.index("GQ")].isdigit():
        return int(values[keys.index("GQ")])
    return None


def scan_vcf_coding(vcf_path: Path, transcripts: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    """Variantes PASS que caen en CDS o a +-2 pb de una unión exón-exón de algún transcrito MANE."""
    regions, starts, max_len = _regions_index(transcripts)
    hits: list[dict[str, Any]] = []
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8") as handle:  # type: ignore[operator]
        for line in handle:
            if line[0] == "#":
                continue
            parts = line.split("\t", 10)
            if parts[6] != "PASS":
                continue
            chrom = parts[0] if parts[0].startswith("chr") else "chr" + parts[0]
            if chrom not in regions:
                continue
            pos = int(parts[1])
            end = pos + len(parts[3]) - 1
            found: set[str] = set()
            j = bisect.bisect_right(starts[chrom], end) - 1
            while j >= 0 and starts[chrom][j] >= pos - max_len[chrom] - 1:
                start_r, end_r, tid = regions[chrom][j]
                if start_r <= end and end_r >= pos:
                    found.add(tid)
                j -= 1
            if found:
                genotype = parts[9].split(":")[0].strip() if len(parts) > 9 else "./."
                gq = _format_gq(parts[8], parts[9]) if len(parts) > 9 else None
                hits.append({
                    "chrom": chrom, "pos": pos, "ref": parts[3], "alt": parts[4], "qual": float(parts[5]),
                    "gt": genotype, "gq": gq, "in_dbsnp": "DB" in parts[7].split(";"), "tids": sorted(found),
                })
    return hits


def compute_effects(hits: list[dict[str, Any]], transcripts: dict[str, dict[str, Any]], cds_sequences: dict[str, str]) -> list[dict[str, Any]]:
    """Una fila por (variante, alelo alternativo presente en el genotipo, transcrito)."""
    effects: list[dict[str, Any]] = []
    for hit in hits:
        alts = hit["alt"].split(",")
        indexes = [int(x) for x in hit["gt"].replace("|", "/").split("/") if x.isdigit() and int(x) > 0]
        for alt_index in sorted(set(indexes)):
            allele = alts[alt_index - 1]
            if allele == "*" or allele.startswith("<"):
                continue
            for tid in hit["tids"]:
                if tid not in cds_sequences:
                    continue
                effect = classify_allele(hit["pos"], hit["ref"], allele, transcripts[tid], cds_sequences[tid])
                effects.append({
                    "chrom": hit["chrom"], "pos": hit["pos"], "ref": hit["ref"], "alt": allele, "qual": hit["qual"],
                    "gt": hit["gt"], "gq": hit.get("gq"), "copies": indexes.count(alt_index), "in_dbsnp": hit["in_dbsnp"], "tid": tid,
                    "gene": transcripts[tid]["gene"], **effect,
                })
    return effects


def variant_key(effect: dict[str, Any]) -> str:
    """Clave chrom:pos:ref:alt tal como figura en el VCF."""
    return f"{effect['chrom']}:{effect['pos']}:{effect['ref']}:{effect['alt']}"


def is_rare_established(effect: dict[str, Any], frequencies: dict[str, dict[str, Any]]) -> bool:
    """Rareza establecida con gnomAD v4 consultado por alelo: ausente de gnomAD, o AF máximo < umbral.

    Una variante cuya frecuencia no se consultó no se considera rara: la rareza no se supone.
    """
    freq = frequencies.get(variant_key(effect))
    if not freq or "present" not in freq:
        return False
    if not freq["present"]:
        return True
    return bool(freq["af"] < RARE_AF_THRESHOLD)


def _allele_prob(effect: dict[str, Any]) -> float:
    af = effect.get("af")
    return max(af, ABSENT_AF_FLOOR) if af is not None else ABSENT_AF_FLOOR


def rank_recessive_candidates(
    effects: list[dict[str, Any]], frequencies: dict[str, dict[str, Any]], hemizygous_x: bool = False
) -> list[dict[str, Any]]:
    """Candidatos recesivos por gen, ordenados por (nº HIGH, log10 de la rareza conjunta).

    Con `hemizygous_x`, en chrX fuera de las PAR solo cuenta el genotipo homocigoto (hemicigoto): la
    probabilidad del genotipo es f, no f², y las llamadas heterocigotas se descartan por inverosímiles.
    """
    by_gene: dict[str, dict[tuple[int, str, str], dict[str, Any]]] = collections.defaultdict(dict)
    for effect in effects:
        if effect["impact"] not in ("HIGH", "MODERATE") or not effect["gene"]:
            continue
        if HYPERVARIABLE_GENES.match(effect["gene"]) or not is_rare_established(effect, frequencies):
            continue
        if hemizygous_x and is_x_nonpar(effect["chrom"], effect["pos"]) and effect["copies"] < 2:
            continue
        key = (effect["pos"], effect["ref"], effect["alt"])
        freq = frequencies.get(variant_key(effect), {})
        entry = by_gene[effect["gene"]].setdefault(key, {**effect, "af": freq.get("af") if freq.get("present") else None})
        entry["copies"] = max(entry["copies"], effect["copies"])
    ranked: list[dict[str, Any]] = []
    for gene, variants in by_gene.items():
        pool = list(variants.values())
        best: tuple[tuple[int, float], dict[str, Any], dict[str, Any]] | None = None
        for i in range(len(pool)):
            for j in range(i + 1, len(pool)):
                a, b = pool[i], pool[j]
                same_site_trans = a["pos"] == b["pos"] and a["alt"] != b["alt"]  # genotipo 1/2: dos alelos del mismo sitio, en trans por definición
                if abs(a["pos"] - b["pos"]) < CIS_PROXIMITY_BP and not same_site_trans:
                    continue
                n_high = (a["impact"] == "HIGH") + (b["impact"] == "HIGH")
                log_bf = math.log10(1.0 / (4 * _allele_prob(a) * _allele_prob(b)))
                if best is None or (n_high, log_bf) > best[0]:
                    best = ((n_high, log_bf), a, b)
        if best:
            ranked.append({"gene": gene, "type": "compound_het", "n_high": best[0][0], "log10_bf": best[0][1], "variants": [best[1], best[2]]})
        for a in pool:
            if a["copies"] >= 2:
                hemizygous = hemizygous_x and is_x_nonpar(a["chrom"], a["pos"])
                ranked.append({"gene": gene, "type": "hemizygous" if hemizygous else "homozygous", "n_high": int(a["impact"] == "HIGH"),
                               "log10_bf": math.log10(1.0 / (_allele_prob(a) if hemizygous else _allele_prob(a) ** 2)), "variants": [a]})
    ranked.sort(key=lambda c: (c["n_high"], c["log10_bf"]), reverse=True)
    for rank, candidate in enumerate(ranked, 1):
        candidate["rank"] = rank
    return ranked


def qualifying_gene_count(effects: list[dict[str, Any]], frequencies: dict[str, dict[str, Any]]) -> int:
    """Genes con al menos una variante HIGH/MODERATE rara (tras los mismos filtros del ranking)."""
    return len({
        e["gene"] for e in effects
        if e["impact"] in ("HIGH", "MODERATE") and e["gene"]
        and not HYPERVARIABLE_GENES.match(e["gene"]) and is_rare_established(e, frequencies)
    })
