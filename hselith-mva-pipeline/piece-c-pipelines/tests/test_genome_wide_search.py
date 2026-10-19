"""Pruebas de la búsqueda genome-wide con transcritos sintéticos (sin red ni datos de paciente)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from src.modules.sv_calling.consequence_calculator import classify_allele
from src.modules.sv_calling.genome_wide_candidates import pair_rarity_score
from src.modules.sv_calling.genome_wide_search import (
    infer_hemizygous_x,
    load_mane_transcripts,
    rank_recessive_candidates,
    scan_vcf_coding,
    variant_key,
)

# Gen en hebra +: CDS "ATG AAA TTT GGC" + parada "TAA" en chr1:100-114; exón UTR extra en 200-230.
PLUS = {
    "chrom": "chr1", "strand": "+", "gene": "PLUSG",
    "exons": [(90, 114), (200, 230)], "cds": [(100, 111, 0)], "stop": [(112, 114)],
}
PLUS_CDS = "ATGAAATTTGGCTAA"
# Gen en hebra -: misma proteína, leída desde chr1:314 hacia 300 (parada en 300-302).
MINUS = {
    "chrom": "chr1", "strand": "-", "gene": "MINUSG",
    "exons": [(290, 314)], "cds": [(303, 314, 0)], "stop": [(300, 302)],
}


def test_snv_consequences_plus_strand() -> None:
    assert classify_allele(103, "A", "T", PLUS, PLUS_CDS)["consequence"] == "stop_gained"
    missense = classify_allele(103, "A", "G", PLUS, PLUS_CDS)
    assert (missense["consequence"], missense["impact"], missense["hgvsp"]) == ("missense_variant", "MODERATE", "p.Lys2Glu")
    assert classify_allele(105, "A", "G", PLUS, PLUS_CDS)["consequence"] == "synonymous_variant"
    assert classify_allele(100, "A", "G", PLUS, PLUS_CDS)["consequence"] == "start_lost"
    assert classify_allele(100, "C", "G", PLUS, PLUS_CDS)["consequence"] == "reference_mismatch"


def test_snv_consequences_minus_strand() -> None:
    assert classify_allele(314, "T", "C", MINUS, PLUS_CDS)["consequence"] == "start_lost"
    assert classify_allele(311, "T", "A", MINUS, PLUS_CDS)["consequence"] == "stop_gained"


def test_indel_and_splice_consequences() -> None:
    assert classify_allele(104, "A", "AT", PLUS, PLUS_CDS)["consequence"] == "frameshift_variant"
    assert classify_allele(103, "AAA", "A", PLUS, PLUS_CDS)["consequence"] == "frameshift_variant"
    assert classify_allele(104, "A", "ATTT", PLUS, PLUS_CDS)["consequence"] == "inframe_insertion"
    assert classify_allele(115, "G", "A", PLUS, PLUS_CDS)["consequence"] == "splice_donor_variant"
    assert classify_allele(199, "G", "A", PLUS, PLUS_CDS)["consequence"] == "splice_acceptor_variant"
    assert classify_allele(150, "G", "A", PLUS, PLUS_CDS)["impact"] == "MODIFIER"


def test_indel_boundary_cases() -> None:
    # deleción del codón de inicio (ancla en la UTR 5')
    assert classify_allele(99, "NATG".replace("N", "C"), "C", PLUS, PLUS_CDS)["consequence"] == "start_lost"
    # inserción entre la UTR y la primera base de la CDS: no altera la proteína
    assert classify_allele(99, "C", "CGCC", PLUS, PLUS_CDS)["impact"] == "MODIFIER"
    # deleción dentro de la UTR 5' (el ancla no cuenta como base afectada)
    assert classify_allele(92, "AT", "A", PLUS, PLUS_CDS)["impact"] == "MODIFIER"
    # deleción en marco del codón de parada conservado -> stop retenido; deleción que lo elimina -> stop_lost
    assert classify_allele(111, "CTAA", "C", PLUS, PLUS_CDS)["consequence"] == "stop_lost"
    # inserción que crea una parada prematura en marco
    assert classify_allele(105, "A", "ATAA", PLUS, PLUS_CDS)["consequence"] == "stop_gained"
    # la referencia se comprueba también en indels
    assert classify_allele(103, "GGG", "G", PLUS, PLUS_CDS)["consequence"] == "reference_mismatch"


def test_mnv_is_translated_not_assumed_missense() -> None:
    assert classify_allele(103, "AAA", "AAG", PLUS, PLUS_CDS)["consequence"] == "synonymous_variant"
    assert classify_allele(103, "AA", "TA", PLUS, PLUS_CDS)["consequence"] == "stop_gained"
    assert classify_allele(103, "AA", "GA", PLUS, PLUS_CDS)["consequence"] == "missense_variant"
    assert classify_allele(103, "AC", "GG", PLUS, PLUS_CDS)["consequence"] == "reference_mismatch"


ABSENT = {"present": False, "af": None, "homozygotes": 0}


def _effect(gene: str, pos: int, impact: str, copies: int = 1, in_dbsnp: bool = False, alt: str = "G") -> dict[str, Any]:
    return {"gene": gene, "chrom": "chr1", "pos": pos, "ref": "A", "alt": alt, "impact": impact, "consequence": "x",
            "copies": copies, "in_dbsnp": in_dbsnp, "gt": "0/1", "qual": 500.0, "hgvsp": None}


def _absent(*effects: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {variant_key(e): ABSENT for e in effects}


def test_ranking_prefers_pair_with_high_impact_variant() -> None:
    effects = [_effect("A", 1000, "HIGH"), _effect("A", 5000, "MODERATE"),
               _effect("B", 1000, "MODERATE"), _effect("B", 5000, "MODERATE"),
               _effect("C", 1000, "MODERATE", copies=2)]
    ranked = rank_recessive_candidates(effects, _absent(*effects))
    assert ranked[0]["gene"] == "A" and ranked[0]["type"] == "compound_het" and ranked[0]["rank"] == 1


def test_pairs_within_cis_distance_are_not_compound_het() -> None:
    effects = [_effect("A", 1000, "HIGH"), _effect("A", 1100, "MODERATE")]
    assert rank_recessive_candidates(effects, _absent(*effects)) == []


def test_two_alleles_at_same_site_are_in_trans() -> None:
    effects = [_effect("A", 1000, "HIGH", alt="G"), _effect("A", 1000, "HIGH", alt="C")]
    ranked = rank_recessive_candidates(effects, _absent(*effects))
    assert [c["gene"] for c in ranked] == ["A"] and ranked[0]["type"] == "compound_het"


def test_rarity_must_be_established() -> None:
    unqueried = _effect("A", 1000, "HIGH", copies=2)
    common = _effect("B", 2000, "HIGH", copies=2)
    rare_known = _effect("C", 3000, "HIGH", copies=2)
    absent = _effect("D", 4000, "HIGH", copies=2)
    freqs = {variant_key(common): {"present": True, "af": 0.2, "homozygotes": 9},
             variant_key(rare_known): {"present": True, "af": 1e-5, "homozygotes": 0},
             variant_key(absent): ABSENT}
    genes = {c["gene"] for c in rank_recessive_candidates([unqueried, common, rare_known, absent], freqs)}
    assert genes == {"C", "D"}


def test_hypervariable_gene_families_are_excluded() -> None:
    effect = _effect("HLA-DQA1", 1000, "HIGH", copies=2)
    assert rank_recessive_candidates([effect], _absent(effect)) == []


def test_pair_rarity_score_saturates_for_rare_pairs_only() -> None:
    log_bf_rare, post_rare = pair_rarity_score(7.9e-5, None, 212)
    _, post_common = pair_rarity_score(0.2, 0.3, 212)
    assert log_bf_rare > 9 and post_rare > 0.999
    assert post_common < 0.05


def test_scan_vcf_keeps_only_pass_variants_in_coding_regions(tmp_path: Path) -> None:
    attrs = 'gene_id "G1"; transcript_id "T1"; gene_name "PLUSG"; tag "MANE_Select";'
    gtf_lines = [
        ("exon", 90, 114, "."), ("exon", 200, 230, "."), ("CDS", 100, 111, "0"), ("stop_codon", 112, 114, "0"),
    ]
    gtf = tmp_path / "mane.gtf"
    gtf.write_text("\n".join(f"chr1\tsrc\t{f}\t{s}\t{e}\t.\t+\t{fr}\t{attrs}" for f, s, e, fr in gtf_lines) + "\n", encoding="utf-8")
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "1\t103\t.\tA\tT\t500\tPASS\tDP=30\tGT\t0/1\n"      # codificante, PASS
        "1\t105\t.\tA\tG\t500\tLowQual\tDP=30\tGT\t0/1\n"   # codificante, no PASS
        "1\t150\t.\tA\tG\t500\tPASS\tDP=30\tGT\t0/1\n"      # intrón profundo
        "1\t115\trs1\tG\tA\t500\tPASS\tDB;DP=30\tGT\t1/1\n",  # sitio de splicing, en dbSNP
        encoding="utf-8",
    )
    hits = scan_vcf_coding(vcf, load_mane_transcripts(gtf))
    assert [(h["pos"], h["gt"], h["in_dbsnp"]) for h in hits] == [(103, "0/1", False), (115, "1/1", True)]


def _x_effect(pos: int, copies: int, alt: str = "G") -> dict[str, Any]:
    return {**_effect("XG", pos, "HIGH", copies=copies, alt=alt), "chrom": "chrX"}


def test_hemizygous_x_uses_single_allele_probability_and_drops_hets() -> None:
    hom = _x_effect(50_000_000, 2)
    het = _x_effect(60_000_000, 1)
    freqs = {variant_key(hom): {"present": True, "af": 1e-4, "homozygotes": 0}, variant_key(het): ABSENT}
    diploid = rank_recessive_candidates([hom], freqs)
    male = rank_recessive_candidates([hom, het], freqs, hemizygous_x=True)
    assert diploid[0]["type"] == "homozygous" and male[0]["type"] == "hemizygous"
    assert diploid[0]["log10_bf"] == pytest.approx(8.0) and male[0]["log10_bf"] == pytest.approx(4.0)
    assert len(male) == 1  # la llamada heterocigota en chrX de un varón no forma par


def test_par_regions_are_not_hemizygous() -> None:
    par = _x_effect(1_000_000, 2)
    freqs = {variant_key(par): ABSENT}
    assert rank_recessive_candidates([par], freqs, hemizygous_x=True)[0]["type"] == "homozygous"


def test_hemizygous_x_inferred_from_genotypes() -> None:
    male = [{"chrom": "chrX", "pos": 50_000_000 + i, "gt": "1/1"} for i in range(30)]
    female = male[:15] + [{"chrom": "chrX", "pos": 50_000_100 + i, "gt": "0/1"} for i in range(15)]
    assert infer_hemizygous_x(male) is True
    assert infer_hemizygous_x(female) is False
    assert infer_hemizygous_x(male[:5]) is False  # sin llamadas suficientes no se asume
