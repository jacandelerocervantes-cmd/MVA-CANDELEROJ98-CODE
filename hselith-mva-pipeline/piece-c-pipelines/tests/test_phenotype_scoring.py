"""Pruebas de la concordancia fenotipo-gen con una ontología y anotaciones sintéticas."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from src.modules.sv_calling.genome_wide_search import rank_recessive_candidates
from src.modules.sv_calling.phenotype_scoring import (
    information_content,
    leave_one_term_out,
    omission_robustness,
    load_gene_annotations,
    load_hpo_ontology,
    load_proband_terms,
    phenotype_rank_among_all_genes,
    phenotype_scores,
    rank_by_phenotype,
)

OBO = """format-version: 1.2
data-version: hp/releases/test

[Term]
id: HP:0000001
name: All

[Term]
id: HP:0000002
name: Branch
is_a: HP:0000001 ! All

[Term]
id: HP:0000003
name: Leaf under branch
alt_id: HP:0000099
is_a: HP:0000002 ! Branch

[Term]
id: HP:0000004
name: Sibling leaf
is_a: HP:0000001 ! All

[Term]
id: HP:0000005
name: Obsolete
is_obsolete: true
is_a: HP:0000004 ! Sibling leaf
"""
GENES = "ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id\n" \
        "1\tGENEA\tHP:0000003\tx\t-\tOMIM:1\n" \
        "2\tGENEB\tHP:0000004\tx\t-\tOMIM:2\n" \
        "3\tGENEC\tHP:0000002\tx\t-\tOMIM:3\n"


@pytest.fixture()
def hpo(tmp_path: Path):
    obo, genes = tmp_path / "hp.obo", tmp_path / "g2p.txt"
    obo.write_text(OBO, encoding="utf-8")
    genes.write_text(GENES, encoding="utf-8")
    ancestors, version = load_hpo_ontology(obo)
    annotations = load_gene_annotations(genes)
    return ancestors, version, annotations, information_content(ancestors, annotations)


def test_ontology_closure_alt_ids_and_version(hpo) -> None:
    ancestors, version, _, _ = hpo
    assert version == "hp/releases/test"
    assert ancestors["HP:0000003"] == {"HP:0000003", "HP:0000002", "HP:0000001"}
    assert ancestors["HP:0000099"] == ancestors["HP:0000003"]  # alt_id resuelto


def test_information_content_follows_gene_frequency(hpo) -> None:
    _, _, _, ic = hpo
    assert ic["HP:0000001"] == pytest.approx(0.0)
    assert ic["HP:0000002"] == pytest.approx(-math.log2(2 / 3))
    assert ic["HP:0000003"] == pytest.approx(-math.log2(1 / 3))


def test_scores_exact_ancestor_and_unrelated_matches(hpo) -> None:
    ancestors, _, annotations, ic = hpo
    scores = phenotype_scores(["HP:0000003"], ancestors, annotations, ic)
    assert scores["GENEA"] == pytest.approx(1.0)                                  # término exacto
    assert scores["GENEC"] == pytest.approx(ic["HP:0000002"] / ic["HP:0000003"])  # solo un ancestro: crédito parcial
    assert scores["GENEB"] == pytest.approx(0.0)                                  # rama sin relación


def test_descendant_annotation_covers_ancestor_term(hpo) -> None:
    ancestors, _, annotations, ic = hpo
    scores = phenotype_scores(["HP:0000002"], ancestors, annotations, ic)
    assert scores["GENEA"] == pytest.approx(1.0)  # anotado a un descendiente del término del probando


def test_rank_by_phenotype_orders_candidates_and_breaks_ties_by_variant_rule() -> None:
    ranked = [
        {"gene": "X", "rank": 1, "n_high": 1, "log10_bf": 12.0, "type": "homozygous"},
        {"gene": "Y", "rank": 2, "n_high": 1, "log10_bf": 9.5, "type": "compound_het"},
        {"gene": "Z", "rank": 3, "n_high": 0, "log10_bf": 8.0, "type": "compound_het"},
    ]
    final = rank_by_phenotype(ranked, {"X": 0.4, "Y": 0.6, "Z": 0.6})
    assert [c["gene"] for c in final] == ["Y", "Z", "X"]          # Y y Z empatan en fenotipo; Y gana por n_high
    assert [c["variant_rank"] for c in final] == [2, 3, 1]        # se conserva el rango por variantes
    assert [c["rank"] for c in final] == [1, 2, 3]


def test_phenotype_rank_among_all_genes_and_leave_one_out(hpo) -> None:
    ancestors, _, annotations, ic = hpo
    scores = phenotype_scores(["HP:0000003", "HP:0000004"], ancestors, annotations, ic)
    assert phenotype_rank_among_all_genes("GENEA", scores) in (1, 2)
    assert phenotype_rank_among_all_genes("NOPE", scores) is None
    loo = leave_one_term_out(["HP:0000003", "HP:0000004"], ["GENEA", "GENEB"], ancestors, annotations, ic)
    assert loo["HP:0000004"]["first"] == "GENEA" and loo["HP:0000003"]["first"] == "GENEB"


def test_shipped_proband_terms_file_has_eight_hpo_ids() -> None:
    path = Path(__file__).resolve().parents[1] / "src" / "modules" / "sv_calling" / "reference_data" / "proband_hpo_terms.json"
    if not path.exists():
        pytest.skip("proband_hpo_terms.json is patient-derived and not published")
    terms = load_proband_terms(path)
    assert len(terms) == 8 and all(t.startswith("HP:") for t in terms)


def test_variant_ranking_still_available_for_the_phenotype_stage() -> None:
    effects = [{"gene": "A", "chrom": "chr1", "pos": 1000, "ref": "A", "alt": "G", "impact": "HIGH", "consequence": "x",
                "copies": 2, "in_dbsnp": False, "gt": "1/1", "qual": 500.0, "hgvsp": None}]
    freqs = {"chr1:1000:A:G": {"present": False, "af": None, "homozygotes": 0}}
    ranked = rank_recessive_candidates(effects, freqs)
    assert ranked[0]["rank"] == 1 and "variant_rank" not in ranked[0]
    assert rank_by_phenotype(ranked, {"A": 0.5})[0]["variant_rank"] == 1


def test_omission_robustness_identifies_the_term_that_flips_the_order(hpo) -> None:
    ancestors, _, annotations, ic = hpo
    terms = ["HP:0000003", "HP:0000004"]
    result = omission_robustness(terms, "GENEA", "GENEB", ancestors, annotations, ic, max_omitted=1)
    assert result["subsets_tested"] == 2
    assert result["subsets_where_challenger_ties_or_leads"] == 1
    assert result["terms_present_in_every_such_subset"] == ["HP:0000003"]
