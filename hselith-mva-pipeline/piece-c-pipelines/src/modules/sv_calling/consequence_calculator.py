"""Cálculo de consecuencias funcionales a partir de datos, sin anotación externa.

Entradas: modelo de transcrito (MANE Select: exones, CDS, codón de parada, hebra), secuencia CDS de
referencia y alelos del VCF. Reglas de Sequence Ontology (codón, marco de lectura, codón de inicio y
de parada, sitios de splicing a +-2 pb). Ninguna etiqueta proviene de una tabla por posición.
"""

from __future__ import annotations

from typing import Any

_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE: dict[str, str] = {
    a + b + c: _AMINO[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
AMINO_3: dict[str, str] = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "E": "Glu", "Q": "Gln", "G": "Gly",
    "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
    "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "*",
}


def reverse_complement(sequence: str) -> str:
    """Complemento inverso de una secuencia de ADN."""
    return sequence.translate(_COMPLEMENT)[::-1]


def coding_segments(transcript: dict[str, Any]) -> list[tuple[int, int]]:
    """Segmentos codificantes (CDS + codón de parada) en orden de transcripción."""
    segments = [(s, e) for s, e, _ in transcript["cds"]] + [tuple(x) for x in transcript["stop"]]
    ordered = sorted(set(segments))
    return ordered if transcript["strand"] == "+" else ordered[::-1]


def genomic_to_cds_offset(pos: int, transcript: dict[str, Any]) -> int | None:
    """Posición (base 0) dentro de la CDS para una coordenada genómica; None si no es codificante."""
    offset = 0
    for start, end in coding_segments(transcript):
        if start <= pos <= end:
            return offset + (pos - start if transcript["strand"] == "+" else end - pos)
        offset += end - start + 1
    return None


def splice_site_class(pos: int, transcript: dict[str, Any]) -> str | None:
    """'splice_donor_variant' o 'splice_acceptor_variant' si pos está a 1-2 pb de una unión exón-exón."""
    exons = sorted(transcript["exons"])
    plus = transcript["strand"] == "+"
    for i in range(len(exons) - 1):
        left_end, right_start = exons[i][1], exons[i + 1][0]
        if left_end + 1 <= pos <= left_end + 2:
            return "splice_donor_variant" if plus else "splice_acceptor_variant"
        if right_start - 2 <= pos <= right_start - 1:
            return "splice_acceptor_variant" if plus else "splice_donor_variant"
    return None


def _snv_effect(pos: int, ref: str, alt: str, transcript: dict[str, Any], cds_seq: str | None) -> dict[str, Any]:
    out: dict[str, Any] = {"consequence": "other", "impact": "MODIFIER", "hgvsp": None}
    offset = genomic_to_cds_offset(pos, transcript)
    if offset is None:
        site = splice_site_class(pos, transcript)
        if site:
            out.update(consequence=site, impact="HIGH")
        return out
    if cds_seq is None:
        return out
    plus = transcript["strand"] == "+"
    ref_base = ref if plus else reverse_complement(ref)
    alt_base = alt if plus else reverse_complement(alt)
    if cds_seq[offset] != ref_base:
        out["consequence"] = "reference_mismatch"
        return out
    codon_index = offset // 3
    ref_codon = cds_seq[codon_index * 3: codon_index * 3 + 3]
    alt_codon = ref_codon[: offset % 3] + alt_base + ref_codon[offset % 3 + 1:]
    ref_aa, alt_aa = CODON_TABLE.get(ref_codon, "?"), CODON_TABLE.get(alt_codon, "?")
    label = f"p.{AMINO_3.get(ref_aa, '?')}{codon_index + 1}{AMINO_3.get(alt_aa, '?')}"
    if codon_index == 0 and ref_codon == "ATG" and alt_codon != "ATG":
        out.update(consequence="start_lost", impact="HIGH", hgvsp="p.Met1?")
    elif ref_aa == alt_aa:
        out.update(consequence="synonymous_variant", impact="LOW")
    elif alt_aa == "*":
        out.update(consequence="stop_gained", impact="HIGH", hgvsp=label)
    elif ref_aa == "*":
        out.update(consequence="stop_lost", impact="HIGH", hgvsp=label)
    else:
        out.update(consequence="missense_variant", impact="MODERATE", hgvsp=label)
    return out


def normalize_indel(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Representación mínima: quita el prefijo y sufijo comunes (ref vacío = inserción antes de pos)."""
    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    ref, alt, pos = ref[i:], alt[i:], pos + i
    j = 0
    while j < len(ref) and j < len(alt) and ref[-1 - j] == alt[-1 - j]:
        j += 1
    return pos, ref[: len(ref) - j], alt[: len(alt) - j]


def _first_stop_codon(sequence: str) -> int | None:
    for index in range(0, len(sequence) - 2, 3):
        if CODON_TABLE.get(sequence[index: index + 3]) == "*":
            return index // 3
    return None


def _translate_indel(offset: int, ref_len: int, alt_seq: str, cds_seq: str) -> dict[str, Any]:
    """Aplica el indel sobre la CDS (orientación del transcrito) y traduce el resultado."""
    out: dict[str, Any] = {"consequence": "other", "impact": "MODIFIER", "hgvsp": None}
    mutated = cds_seq[:offset] + alt_seq + cds_seq[offset + ref_len:]
    diff = len(alt_seq) - ref_len
    stop_index = len(cds_seq) // 3 - 1
    touches_start = offset < 3 and cds_seq[:3] == "ATG"
    if touches_start and mutated[:3] != "ATG":
        out.update(consequence="start_lost", impact="HIGH")
    elif offset + ref_len > len(cds_seq) - 3:  # toca el codón de parada
        new_stop = _first_stop_codon(mutated)
        if new_stop is None:
            out.update(consequence="stop_lost", impact="HIGH")
        elif diff % 3 == 0 and new_stop == stop_index + diff // 3:
            out.update(consequence="stop_retained_variant", impact="LOW")
        elif new_stop < stop_index + diff // 3:
            out.update(consequence="stop_gained", impact="HIGH")
        else:
            out.update(consequence="stop_lost", impact="HIGH")
    elif diff % 3 != 0:
        out.update(consequence="frameshift_variant", impact="HIGH")
    else:
        new_stop = _first_stop_codon(mutated)
        if new_stop is not None and new_stop < stop_index:
            out.update(consequence="stop_gained", impact="HIGH")
        else:
            out.update(consequence="inframe_insertion" if diff > 0 else "inframe_deletion", impact="MODERATE")
    return out


def _indel_effect(pos: int, ref: str, alt: str, transcript: dict[str, Any], cds_seq: str | None) -> dict[str, Any]:
    out: dict[str, Any] = {"consequence": "other", "impact": "MODIFIER", "hgvsp": None}
    pos, ref, alt = normalize_indel(pos, ref, alt)
    plus = transcript["strand"] == "+"
    if ref:  # deleción o sustitución compleja: bases afectadas pos..pos+len(ref)-1
        bases = list(range(pos, pos + len(ref)))
    else:  # inserción entre pos-1 y pos
        bases = [pos - 1, pos]
    site = next((c for c in (splice_site_class(b, transcript) for b in bases) if c), None)
    offsets = [genomic_to_cds_offset(b, transcript) for b in bases]
    coding = [o for o in offsets if o is not None]
    contiguous = len(coding) == len(offsets) and all(abs(offsets[i + 1] - offsets[i]) == 1 for i in range(len(offsets) - 1))
    if cds_seq is None or not coding:
        if site:
            out.update(consequence=site, impact="HIGH")
        return out
    if site and ref:
        out.update(consequence=site, impact="HIGH")
        return out
    if not contiguous:  # cruza el borde CDS/UTR o un intrón
        if site:
            out.update(consequence=site, impact="HIGH")
        elif not ref:  # inserción justo en el límite entre UTR y CDS: no altera la proteína
            return out
        elif any(o < 3 for o in coding) and cds_seq[:3] == "ATG":
            out.update(consequence="start_lost", impact="HIGH")
        elif any(o >= len(cds_seq) - 3 for o in coding):
            out.update(consequence="stop_lost", impact="HIGH")
        else:
            out.update(consequence="coding_sequence_variant", impact="MODERATE")
        return out
    if ref:
        start = min(offsets) if not plus else offsets[0]
        ref_t = ref if plus else reverse_complement(ref)
        alt_t = alt if plus else reverse_complement(alt)
        if cds_seq[start: start + len(ref_t)] != ref_t:
            out["consequence"] = "reference_mismatch"
            return out
        return _translate_indel(start, len(ref_t), alt_t, cds_seq)
    insert_at = offsets[1] if plus else offsets[0]  # índice en la CDS donde se inserta
    return _translate_indel(insert_at, 0, alt if plus else reverse_complement(alt), cds_seq)


def _mnv_effect(pos: int, ref: str, alt: str, transcript: dict[str, Any], cds_seq: str | None) -> dict[str, Any]:
    """Sustitución de varias bases: se aplica a la CDS y se comparan los codones afectados (no se asume 'missense')."""
    out: dict[str, Any] = {"consequence": "other", "impact": "MODIFIER", "hgvsp": None}
    plus = transcript["strand"] == "+"
    changes: dict[int, str] = {}
    for i, (r, a) in enumerate(zip(ref, alt)):
        offset = genomic_to_cds_offset(pos + i, transcript)
        if offset is None:
            site = splice_site_class(pos + i, transcript)
            if site and r != a:
                out.update(consequence=site, impact="HIGH")
            continue
        if cds_seq is None:
            return out
        r_c, a_c = (r, a) if plus else (reverse_complement(r), reverse_complement(a))
        if cds_seq[offset] != r_c:
            return {**out, "consequence": "reference_mismatch"}
        if r_c != a_c:
            changes[offset] = a_c
    if not changes or cds_seq is None:
        return out
    results = []
    for ci in sorted({o // 3 for o in changes}):
        ref_codon = cds_seq[ci * 3: ci * 3 + 3]
        alt_codon = "".join(changes.get(ci * 3 + k, ref_codon[k]) for k in range(len(ref_codon)))
        ref_aa, alt_aa = CODON_TABLE.get(ref_codon, "?"), CODON_TABLE.get(alt_codon, "?")
        label = f"p.{AMINO_3.get(ref_aa, '?')}{ci + 1}{AMINO_3.get(alt_aa, '?')}"
        if ci == 0 and ref_codon == "ATG" and alt_codon != "ATG":
            results.append(("start_lost", "HIGH", "p.Met1?"))
        elif alt_aa == "*" and ref_aa != "*":
            results.append(("stop_gained", "HIGH", label))
        elif ref_aa == "*" and alt_aa != "*":
            results.append(("stop_lost", "HIGH", label))
        elif ref_aa != alt_aa:
            results.append(("missense_variant", "MODERATE", label))
    order = {"HIGH": 0, "MODERATE": 1}
    if results:
        c, imp, label = min(results, key=lambda r: order[r[1]])
        out.update(consequence=c, impact=imp, hgvsp=label)
    elif out["impact"] == "MODIFIER":
        out.update(consequence="synonymous_variant", impact="LOW")
    return out


def classify_allele(pos: int, ref: str, alt: str, transcript: dict[str, Any], cds_seq: str | None) -> dict[str, Any]:
    """Consecuencia calculada de un alelo (SNV o indel) respecto del transcrito."""
    if len(ref) == 1 and len(alt) == 1:
        return _snv_effect(pos, ref, alt, transcript, cds_seq)
    if len(ref) == len(alt):  # sustitución de varias bases
        return _mnv_effect(pos, ref, alt, transcript, cds_seq)
    return _indel_effect(pos, ref, alt, transcript, cds_seq)
