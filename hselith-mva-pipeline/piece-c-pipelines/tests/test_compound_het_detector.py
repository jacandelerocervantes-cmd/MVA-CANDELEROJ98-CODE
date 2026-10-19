"""Pruebas de los hallazgos estructurales secundarios y de las filas candidatas (sin genes fijos)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import jsonschema
import pytest

from src.modules.sv_calling.compound_het_detector import (
    CandidatePair,
    _is_repetitive_context,
    detect_secondary_structural_architectures,
)
from src.modules.sv_calling.genome_wide_candidates import build_primary_candidates, build_secondary_candidate
from src.modules.sv_calling.genome_wide_search import coding_overlap_checker

pytestmark = pytest.mark.track1

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = REPO_ROOT / "shared-contracts" / "schemas" / "variant_ranking.schema.json"
# Contexto no repetitivo y repetitivo (secuencias reales de Manta de este caso)
PLAIN_ALT = "CCCCAGGCCACAGCACAGAGGAACTCGGGAGGAAGGGCCAGGCAGGGGGCCTTGTGAGACCA"
TA_ALT = "TGTGTATATATATATATATATATATATATACATATATATATATATATATATT"


def _sv(gene: str, start: int, genotype: str, qual: float = 999.0, alt: str = PLAIN_ALT, ref: str = "C") -> dict[str, Any]:
    return {"chromosome": "chr1", "start_position": start, "end_position": start + 1, "sv_type": "INS", "sv_length": len(alt) - len(ref),
            "qual_score": qual, "genotype": genotype, "ref_allele": ref, "alt_allele": alt, "genes_affected": [gene]}


@pytest.fixture
def ranking_schema() -> dict[str, Any]:
    return json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))


def test_secondary_score_is_phenotype_times_quality_gate() -> None:
    rows = detect_secondary_structural_architectures([_sv("GENEA", 1000, "1/1", qual=500.0)], {"GENEA": 0.4})
    assert len(rows) == 1 and rows[0].gene == "GENEA" and rows[0].architecture == "homozygous"
    assert rows[0].base_score == pytest.approx(0.4 * 500.0 / 999.0, abs=1e-4)


def test_secondary_architectures_derive_from_input_and_skip_genes_without_phenotype() -> None:
    records = [_sv("GENEA", 1000, "1/1"), _sv("GENEB", 2000, "1/1")]
    rows = detect_secondary_structural_architectures(records, {"GENEA": 0.3})  # GENEB sin puntaje fenotípico
    assert [r.gene for r in rows] == ["GENEA"]
    assert detect_secondary_structural_architectures([], {"GENEA": 0.3}) == []


def test_compound_het_pair_and_strictly_descending_unique_scores() -> None:
    records = [_sv("GENEA", 1000, "0/1", qual=700.0), _sv("GENEA", 9000, "0/1", qual=900.0), _sv("GENEA", 20000, "1/1"), _sv("GENEC", 5000, "1/1")]
    rows = detect_secondary_structural_architectures(records, {"GENEA": 0.3, "GENEC": 0.3})
    pair = next(r for r in rows if r.chrom_2 is not None)
    assert pair.architecture == "compound_heterozygous" and (pair.pos_1, pair.pos_2) == (1000, 9000)
    assert pair.base_score == pytest.approx(0.3 * 700.0 / 999.0, abs=1e-4)  # eslabón más débil
    scores = [r.base_score for r in rows]
    assert scores == sorted(scores, reverse=True) and len(set(scores)) == len(scores)
    assert all(not r.alt_1.startswith("<") for r in rows)


def test_symbolic_and_short_structural_variants_are_excluded() -> None:
    symbolic = _sv("GENEA", 1000, "1/1", alt="<DEL>")
    short = {**_sv("GENEA", 2000, "1/1"), "sv_length": 10}
    assert detect_secondary_structural_architectures([symbolic, short], {"GENEA": 0.5}) == []


def test_repetitive_context_lowers_the_quality_gate() -> None:
    plain = detect_secondary_structural_architectures([_sv("GENEA", 1000, "1/1")], {"GENEA": 0.5})[0]
    repeat = detect_secondary_structural_architectures([_sv("GENEA", 1000, "1/1", alt=TA_ALT, ref="TA")], {"GENEA": 0.5})[0]
    assert repeat.base_score == pytest.approx(plain.base_score * 0.35, abs=1e-3) and "WARNING" in repeat.notes


def test_coding_overlap_is_computed_from_transcript_models_and_noted() -> None:
    transcripts = {"T1": {"chrom": "chr1", "strand": "+", "gene": "GENEA", "exons": [(900, 1100)], "cds": [(950, 1050, 0)], "stop": [(1051, 1053)]}}
    overlaps = coding_overlap_checker(transcripts)
    assert overlaps("chr1", 1000, 1000) and overlaps("chr1", 1052, 1060) and not overlaps("chr1", 100, 200) and not overlaps("chr2", 1000, 1000)
    inside = detect_secondary_structural_architectures([_sv("GENEA", 1000, "1/1")], {"GENEA": 0.5}, overlaps)[0]
    outside = detect_secondary_structural_architectures([_sv("GENEA", 5000, "1/1")], {"GENEA": 0.5}, overlaps)[0]
    assert "coding sequence: yes" in inside.notes and "coding sequence: no" in outside.notes


def test_is_repetitive_context_against_real_manta_sequences() -> None:
    assert _is_repetitive_context("TA", "TGTATATATATATACATATATGTGTATATATATGTGTGTATATATATATACATATATGTATATATATGTGTGTATATATATATACATATATGTATATATATGTGTG")
    assert _is_repetitive_context("TA", TA_ALT)
    assert not _is_repetitive_context("TCTTCCCAGCACCTTGGAAAAGTTCATATCCTGTGTACATTTTGCATTTCTGATTAAAAAGTATAAAGTGAAACAAAAACATTTTCTCCCACTTATCTTGTTCACTACAGT", "TTGATATTA")
    assert not _is_repetitive_context("C", PLAIN_ALT)


def _top(gene: str, phenotype: float, kind: str = "compound_het") -> dict[str, Any]:
    a = {"chrom": "chr9", "pos": 5000, "ref": "T", "alt": "G", "impact": "HIGH", "consequence": "stop_gained", "hgvsp": "p.Leu1*", "qual": 700.0, "gt": "0/1", "af": 1e-5}
    b = {"chrom": "chr9", "pos": 9000, "ref": "T", "alt": "G", "impact": "MODERATE", "consequence": "missense_variant", "hgvsp": "p.Asn2Lys", "qual": 400.0, "gt": "0/1", "af": None}
    return {"gene": gene, "type": kind, "variants": [a, b] if kind == "compound_het" else [a], "rank": 1, "variant_rank": 3, "phenotype_score": phenotype, "n_high": 1, "log10_bf": 9.5}


def test_primary_rows_are_built_for_whatever_gene_ranks_first(ranking_schema: dict[str, Any]) -> None:
    rows = build_primary_candidates(_top("ZZZ9", 0.5), n_candidate_genes=188)
    assert [r.finding_type for r in rows] == ["primary", "secondary", "secondary", "secondary"]
    assert all(r.gene == "ZZZ9" for r in rows) and rows[0].base_score == 1.0
    assert "phase undetermined" in rows[0].notes and "trans" not in rows[0].notes
    assert rows[2].base_score == pytest.approx(0.5 * 700.0 / 999.0, abs=1e-3)
    for row in rows:
        jsonschema.validate(instance=row.to_submission_dict(), schema=ranking_schema)


def test_without_phenotype_support_only_pair_rows_are_emitted() -> None:
    rows = build_primary_candidates(_top("ZZZ9", 0.0), n_candidate_genes=188)
    assert [r.finding_type for r in rows] == ["primary", "secondary"]


def test_homozygous_primary_and_secondary_candidates() -> None:
    rows = build_primary_candidates(_top("HOMG", 0.4, "homozygous"), n_candidate_genes=188)
    assert len(rows) == 1 and rows[0].architecture == "homozygous" and rows[0].finding_type == "primary"
    secondary = build_secondary_candidate({**_top("HOMG", 0.4, "homozygous"), "rank": 2, "type": "homozygous"})
    assert isinstance(secondary, CandidatePair) and secondary.base_score == pytest.approx(0.4 * 700.0 / 999.0, abs=1e-3)
    unknown = build_secondary_candidate({**_top("HOMG", 0.0, "homozygous")})
    assert unknown is not None and unknown.base_score == pytest.approx(0.01 * 700.0 / 999.0, abs=1e-3)
