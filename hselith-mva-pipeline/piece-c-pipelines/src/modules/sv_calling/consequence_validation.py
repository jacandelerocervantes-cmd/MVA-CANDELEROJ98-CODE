"""Validación del cálculo local de consecuencias contra Ensembl VEP (solo como control externo)."""

from __future__ import annotations

from typing import Any

from src.modules.sv_calling.genome_wide_search import variant_key

_HIGH = {"transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant",
         "stop_lost", "start_lost", "transcript_amplification"}
_MODERATE = {"inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant"}
_LOW = {"splice_region_variant", "splice_donor_5th_base_variant", "splice_donor_region_variant",
        "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant", "start_retained_variant",
        "stop_retained_variant", "synonymous_variant"}
_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}


def vep_impact(terms: list[str]) -> str:
    """Impacto máximo de un conjunto de términos SO de VEP."""
    if any(t in _HIGH for t in terms):
        return "HIGH"
    if any(t in _MODERATE for t in terms):
        return "MODERATE"
    if any(t in _LOW for t in terms):
        return "LOW"
    return "MODIFIER"


def validate_against_vep(effects: list[dict[str, Any]], frequencies: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Concordancia por impacto y por término entre el cálculo local y VEP (una fila por variante-alelo)."""
    seen: set[str] = set()
    compared = exact = impact_ok = ref_mismatch = 0
    called = {"HIGH": 0, "MODERATE": 0}
    missed: list[str] = []      # VEP dice HIGH/MODERATE y el cálculo local no
    extra: list[str] = []       # el cálculo local dice HIGH/MODERATE y VEP no
    for effect in effects:
        key = variant_key(effect)
        if key in seen:
            continue
        seen.add(key)
        if effect["consequence"] == "reference_mismatch":
            ref_mismatch += 1
            continue
        terms = frequencies.get(key, {}).get("vep_terms")
        if not terms:
            continue
        compared += 1
        mine, theirs = effect["impact"], vep_impact(terms)
        exact += effect["consequence"] in terms
        impact_ok += mine == theirs
        if mine in called:
            called[mine] += 1
        if theirs in ("HIGH", "MODERATE") and mine not in ("HIGH", "MODERATE"):
            missed.append(key)
        if mine in ("HIGH", "MODERATE") and theirs not in ("HIGH", "MODERATE"):
            extra.append(key)
    return {
        "variant_alleles_compared": compared, "reference_base_mismatches": ref_mismatch,
        "impact_concordance": impact_ok / compared if compared else None,
        "term_concordance": exact / compared if compared else None,
        "local_high_or_moderate": sum(called.values()),
        "vep_high_or_moderate_missed_locally": len(missed), "local_high_or_moderate_not_in_vep": len(extra),
        "missed_examples": missed[:10], "extra_examples": extra[:10],
    }
