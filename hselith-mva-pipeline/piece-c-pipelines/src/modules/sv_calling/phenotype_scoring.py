"""Concordancia fenotipo-gen a partir de HPO (ontología y anotaciones públicas).

Segunda etapa estándar de la priorización diagnóstica: tras seleccionar candidatos por evidencia de
variantes, se ordenan por la similitud entre los términos HPO del probando y los términos HPO
anotados a cada gen (cobertura ponderada por contenido de información, un solo sentido: del probando
hacia el gen). No hay ninguna lista de genes ni pesos escritos a mano.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


def load_hpo_ontology(obo_path: Path) -> tuple[dict[str, set[str]], str]:
    """Ancestros (incluido el propio término) por id de HPO, resolviendo alt_id; devuelve también la versión."""
    parents: dict[str, list[str]] = {}
    alt_to_main: dict[str, str] = {}
    version = ""
    current: str | None = None
    obsolete = False
    for line in obo_path.read_text(encoding="utf-8").splitlines():
        if line.startswith("data-version:"):
            version = line.split(":", 1)[1].strip()
        elif line == "[Term]":
            current, obsolete = None, False
        elif line.startswith("id: HP:"):
            current = line[4:].strip()
            parents.setdefault(current, [])
        elif current and line.startswith("alt_id: HP:"):
            alt_to_main[line[8:].strip()] = current
        elif current and line.startswith("is_obsolete: true"):
            obsolete = True
        elif current and line.startswith("is_a: HP:") and not obsolete:
            parents[current].append(line[6:].split("!")[0].strip())
    ancestors: dict[str, set[str]] = {}

    def closure(term: str) -> set[str]:
        if term not in ancestors:
            ancestors[term] = {term}  # evita recursión en ciclos inesperados
            for parent in parents.get(term, []):
                ancestors[term] |= closure(parent)
        return ancestors[term]

    for term in list(parents):
        closure(term)
    for alt, main in alt_to_main.items():
        ancestors[alt] = ancestors[main]
    return ancestors, version


def load_gene_annotations(genes_to_phenotype_path: Path) -> dict[str, set[str]]:
    """Términos HPO anotados directamente a cada símbolo de gen."""
    annotations: dict[str, set[str]] = defaultdict(set)
    with open(genes_to_phenotype_path, encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            annotations[fields[1]].add(fields[2])
    return dict(annotations)


def information_content(ancestors: dict[str, set[str]], annotations: dict[str, set[str]]) -> dict[str, float]:
    """IC(t) = -log2(fracción de genes anotados a t o a un descendiente de t)."""
    counts: dict[str, int] = defaultdict(int)
    for terms in annotations.values():
        covered: set[str] = set()
        for term in terms:
            covered |= ancestors.get(term, {term})
        for term in covered:
            counts[term] += 1
    total = len(annotations)
    return {term: -math.log2(n / total) for term, n in counts.items()}


def phenotype_scores(
    proband_terms: list[str], ancestors: dict[str, set[str]], annotations: dict[str, set[str]], ic: dict[str, float],
    genes: list[str] | None = None,
) -> dict[str, float]:
    """Puntuación [0, 1] por gen: cobertura ponderada por IC de los términos del probando."""
    weights = {t: ic.get(t, 0.0) for t in proband_terms}
    denominator = sum(weights.values())
    if denominator == 0:
        raise ValueError("Ningún término HPO del probando tiene contenido de información (¿ids obsoletos?).")
    scores: dict[str, float] = {}
    for gene in (genes if genes is not None else annotations):
        terms = annotations.get(gene, set())
        gene_ancestors: set[str] = set()
        for term in terms:
            gene_ancestors |= ancestors.get(term, {term})
        total = 0.0
        for term, weight in weights.items():
            if weight == 0.0:
                continue
            best = max((ic.get(a, 0.0) for a in ancestors.get(term, {term}) if a in gene_ancestors), default=0.0)
            total += weight * min(best / weight, 1.0)
        scores[gene] = total / denominator
    return scores


def load_proband_terms(path: Path) -> list[str]:
    """Ids HPO del probando desde reference_data/proband_hpo_terms.json."""
    return [t["id"] for t in json.loads(path.read_text(encoding="utf-8"))["terms"]]


def rank_by_phenotype(ranked: list[dict[str, Any]], scores: dict[str, float]) -> list[dict[str, Any]]:
    """Orden final de los candidatos: similitud fenotípica, y a igualdad la regla de variantes."""
    for candidate in ranked:
        candidate["variant_rank"] = candidate["rank"]
        candidate["phenotype_score"] = scores.get(candidate["gene"], 0.0)
    final = sorted(ranked, key=lambda c: (c["phenotype_score"], c["n_high"], c["log10_bf"]), reverse=True)
    for position, candidate in enumerate(final, 1):
        candidate["rank"] = position
    return final


def phenotype_rank_among_all_genes(gene: str, scores: dict[str, float]) -> int | None:
    """Posición de un gen entre todos los genes anotados, solo por fenotipo (informativo)."""
    if gene not in scores:
        return None
    return 1 + sum(1 for s in scores.values() if s > scores[gene])


def leave_one_term_out(
    proband_terms: list[str], genes: list[str], ancestors: dict[str, set[str]],
    annotations: dict[str, set[str]], ic: dict[str, float],
) -> dict[str, dict[str, Any]]:
    """Puntuación de los genes dados quitando un término del probando cada vez (análisis de robustez)."""
    report: dict[str, dict[str, Any]] = {}
    for dropped in proband_terms:
        scores = phenotype_scores([t for t in proband_terms if t != dropped], ancestors, annotations, ic)
        ordered = sorted(genes, key=lambda g: scores.get(g, 0.0), reverse=True)
        report[dropped] = {"first": ordered[0], "scores": {g: round(scores.get(g, 0.0), 4) for g in genes}}
    return report


def omission_robustness(
    proband_terms: list[str], leader: str, challenger: str, ancestors: dict[str, set[str]],
    annotations: dict[str, set[str]], ic: dict[str, float], max_omitted: int = 3,
) -> dict[str, Any]:
    """Prueba todos los subconjuntos de hasta `max_omitted` términos omitidos y cuenta cuántas veces el
    segundo candidato empata o supera al primero; informa qué términos están en todos esos subconjuntos."""
    from itertools import combinations

    tested = 0
    flipped: list[tuple[str, ...]] = []
    for size in range(1, max_omitted + 1):
        for dropped in combinations(proband_terms, size):
            tested += 1
            scores = phenotype_scores([t for t in proband_terms if t not in dropped], ancestors, annotations, ic, genes=[leader, challenger])
            if scores.get(challenger, 0.0) >= scores.get(leader, 0.0):
                flipped.append(dropped)
    common = set(proband_terms) if flipped else set()
    for dropped in flipped:
        common &= set(dropped)
    containing = {t: sum(1 for d in combinations_all(proband_terms, max_omitted) if t in d) for t in common}
    return {
        "leader": leader, "challenger": challenger, "max_omitted": max_omitted, "subsets_tested": tested,
        "subsets_where_challenger_ties_or_leads": len(flipped), "terms_present_in_every_such_subset": sorted(common),
        "subsets_omitting_each_such_term": containing,
    }


def combinations_all(terms: list[str], max_omitted: int) -> list[tuple[str, ...]]:
    """Todos los subconjuntos de 1 a `max_omitted` términos."""
    from itertools import combinations

    return [c for size in range(1, max_omitted + 1) for c in combinations(terms, size)]
