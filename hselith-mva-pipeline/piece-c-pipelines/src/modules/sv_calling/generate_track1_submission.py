"""Generator of the official Hackathon Track 1 predictions file.

Produces the ranked candidate variants CSV conforming strictly to
shared-contracts/schemas/variant_ranking.schema.json.

The primary candidate is the top-ranked recessive candidate of a genome-wide search over all
MANE Select coding regions of the input VCF (consequences are computed, not looked up); structural
variants from Manta are ranked as secondary findings within the curated 22-gene MVA panel.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
import sys

# Ensure piece-c-pipelines root is in sys.path automatically
_PIECE_C_ROOT = Path(__file__).resolve().parents[3]
if str(_PIECE_C_ROOT) not in sys.path:
    sys.path.insert(0, str(_PIECE_C_ROOT))

from src.modules.sv_calling.compound_het_detector import _enforce_strictly_descending, detect_secondary_structural_architectures
from src.modules.sv_calling.consequence_validation import validate_against_vep
from src.modules.sv_calling.dual_vcf_parser import parse_dual_vcfs
from src.modules.sv_calling.epcr_ranking_engine import generate_submission_files
from src.modules.sv_calling.genome_wide_candidates import build_primary_candidates, build_secondary_candidate
from src.modules.sv_calling.genome_wide_io import fetch_cds_sequences, fetch_gnomad_frequencies, fetch_vep_terms, validate_cds_lengths
from src.modules.sv_calling.genome_wide_search import (
    coding_overlap_checker,
    compute_effects,
    load_mane_transcripts,
    infer_hemizygous_x,
    qualifying_gene_count,
    rank_recessive_candidates,
    scan_vcf_coding,
    variant_key,
)
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

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[4]
PROJECT_ROOT = REPO_ROOT.parent
OUTPUT_DIR = PROJECT_ROOT / "submissions" / "track1"
REFERENCE_DIR = PROJECT_ROOT / "data" / "reference"
DEFAULT_MANE_GTF = REFERENCE_DIR / "MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz"
DEFAULT_CACHE_DIR = REFERENCE_DIR / "cache"
RANKING_TOP_N = 25
MAX_SUBMISSION_ROWS = 10
DEFAULT_HPO_OBO = REFERENCE_DIR / "hp.obo"
DEFAULT_HPO_GENES = REFERENCE_DIR / "genes_to_phenotype.txt"
PROBAND_TERMS_PATH = Path(__file__).with_name("reference_data") / "proband_hpo_terms.json"


def write_ranking_csv(ranked: list[dict], path: Path, top_n: int = RANKING_TOP_N) -> None:
    """Ranking genome-wide de candidatos recesivos (archivo de soporte, no es la submission)."""
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["rank", "variant_rank", "gene", "type", "n_high_impact", "log10_joint_rarity", "phenotype_score", "variants"])
        for candidate in ranked[:top_n]:
            variants = " | ".join(
                f"{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']} {v['consequence']} "
                f"AF={'absent' if v.get('af') is None else format(v['af'], '.3g')} GT={v['gt']}"
                for v in candidate["variants"]
            )
            writer.writerow([candidate["rank"], candidate["variant_rank"], candidate["gene"], candidate["type"], candidate["n_high"],
                             f"{candidate['log10_bf']:.2f}", f"{candidate['phenotype_score']:.3f}", variants])


def run_track1_pipeline(
    gatk_vcf: Path,
    manta_vcf: Path,
    mane_gtf: Path = DEFAULT_MANE_GTF,
    cache_dir: Path = DEFAULT_CACHE_DIR,
    hpo_obo: Path = DEFAULT_HPO_OBO,
    hpo_genes: Path = DEFAULT_HPO_GENES,
    output_dir: Path = OUTPUT_DIR,
    hf_username: str = "candeleroj98",
) -> Path:
    """Ejecuta el pipeline unificado de Track 1: búsqueda genome-wide + SV secundarios."""
    logger.info("Iniciando pipeline de Track 1...")

    # 1. Búsqueda genome-wide de candidatos recesivos (consecuencias calculadas)
    transcripts = load_mane_transcripts(mane_gtf)
    hits = scan_vcf_coding(gatk_vcf, transcripts)
    cds = fetch_cds_sequences(sorted({t for h in hits for t in h["tids"]}), cache_dir / "cds_sequences.json")
    bad = validate_cds_lengths(transcripts, cds)
    if bad:
        raise ValueError(f"CDS de referencia inconsistente con el GTF en {len(bad)} transcritos (p. ej. {bad[:3]}).")
    effects = compute_effects(hits, transcripts, cds)
    mismatches = [e for e in effects if e["consequence"] == "reference_mismatch"]
    if mismatches:
        raise ValueError(f"{len(mismatches)} variantes con base de referencia distinta del VCF (ensamblado incorrecto?).")
    keys = sorted({variant_key(e) for e in effects})
    candidate_keys = sorted({variant_key(e) for e in effects if e["impact"] in ("HIGH", "MODERATE")})
    frequencies = fetch_gnomad_frequencies(candidate_keys, cache_dir / "gnomad_v4_frequencies.json")
    vep_terms = fetch_vep_terms(keys, cache_dir / "vep_terms.json")
    validation = validate_against_vep(effects, vep_terms)
    logger.info("Validación del cálculo local vs VEP: %s", {k: v for k, v in validation.items() if not k.endswith("examples")})
    hemizygous_x = infer_hemizygous_x(hits)
    logger.info("chrX fuera de las PAR: genotipos %s", "hemicigotos (probando con un solo X)" if hemizygous_x else "diploides")
    ranked = rank_recessive_candidates(effects, frequencies, hemizygous_x=hemizygous_x)

    # 1b. Segunda etapa: similitud fenotípica (HPO) entre los candidatos por variantes
    ancestors, hpo_version = load_hpo_ontology(hpo_obo)
    annotations = load_gene_annotations(hpo_genes)
    ic = information_content(ancestors, annotations)
    if not PROBAND_TERMS_PATH.exists():
        raise FileNotFoundError(
            f"{PROBAND_TERMS_PATH.name} not found: copy proband_hpo_terms.example.json to {PROBAND_TERMS_PATH.name} "
            "and list the proband's HPO terms (the file is patient-derived and not published)."
        )
    proband_terms = load_proband_terms(PROBAND_TERMS_PATH)
    pheno = phenotype_scores(proband_terms, ancestors, annotations, ic)
    candidate_genes = sorted({c["gene"] for c in ranked})
    robustness = leave_one_term_out(proband_terms, candidate_genes, ancestors, annotations, ic)
    ranked = rank_by_phenotype(ranked, pheno)
    top = ranked[0]
    primary_rows = build_primary_candidates(top, n_candidate_genes=qualifying_gene_count(effects, frequencies))
    logger.info("Candidato primario genome-wide: %s (rango final %d de %d; rango por variantes %d; fenotipo %.3f)",
                top["gene"], top["rank"], len(ranked), top["variant_rank"], top["phenotype_score"])

    # 2. Hallazgos estructurales secundarios (Manta, panel curado de 22 genes)
    sv_records = parse_dual_vcfs(gatk_vcf_path=None, manta_vcf_path=manta_vcf)
    secondaries = detect_secondary_structural_architectures(sv_records, pheno, coding_overlap_checker(transcripts))

    # 2b. Candidatos genome-wide siguientes: compiten por las filas libres con los SV por puntaje, no por orden de llegada
    genome_wide_secondaries = [build_secondary_candidate(c) for c in ranked[1:]]

    # 3. Calibración y generación del CSV oficial
    output_dir.mkdir(parents=True, exist_ok=True)
    # Primario y recíproco fijos; el resto se ordena por puntaje y los empates se resuelven con un épsilon mínimo
    leading = primary_rows[:2] if len(primary_rows) >= 2 else primary_rows
    rest = _enforce_strictly_descending(primary_rows[len(leading):] + secondaries + genome_wide_secondaries)[: MAX_SUBMISSION_ROWS - len(leading)]
    candelero_path = generate_submission_files(leading + rest, output_dir=output_dir, hf_username=hf_username)
    write_ranking_csv(ranked, output_dir / f"{hf_username}_genomewide_ranking.csv")
    summary = {
        "chrx_nonpar_hemizygous": hemizygous_x, "hpo_release": hpo_version, "proband_terms": proband_terms, "genes_with_hpo_annotations": len(annotations),
        "candidate_scores": {c["gene"]: round(c["phenotype_score"], 4) for c in ranked},
        "candidate_phenotype_rank_among_all_genes": {g: phenotype_rank_among_all_genes(g, pheno) for g in candidate_genes},
        "top10_genes_by_phenotype_alone": [[g, round(v, 4)] for g, v in sorted(pheno.items(), key=lambda kv: -kv[1])[:10]],
        "leave_one_term_out": robustness,
        "omission_robustness_up_to_3_terms": (
            omission_robustness(proband_terms, ranked[0]["gene"], ranked[1]["gene"], ancestors, annotations, ic) if len(ranked) > 1 else None
        ),
    }
    for name, payload in (("phenotype_summary", summary), ("consequence_validation", validation)):
        with open(output_dir / f"{hf_username}_{name}.json", "w", encoding="utf-8", newline="\n") as handle:
            handle.write(json.dumps(payload, indent=2) + "\n")
    logger.info("Archivo oficial de submission generado: %s", candelero_path)
    return candelero_path


def main() -> None:
    """Punto de entrada CLI para generar los archivos oficiales de submission."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Generador oficial de predicciones Track 1 para HSelith MVA Hackathon 2026."
    )
    parser.add_argument("--gatk-vcf", type=Path, required=True, help="Ruta al VCF de SNVs/indels")
    parser.add_argument("--manta-vcf", type=Path, required=True, help="Ruta al VCF de Manta SVs")
    parser.add_argument("--mane-gtf", type=Path, default=DEFAULT_MANE_GTF, help="GTF MANE Select (GRCh38)")
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR, help="Caché de secuencias CDS y frecuencias")
    parser.add_argument("--hpo-obo", type=Path, default=DEFAULT_HPO_OBO, help="Ontología HPO (hp.obo)")
    parser.add_argument("--hpo-genes", type=Path, default=DEFAULT_HPO_GENES, help="Anotaciones HPO gen-fenotipo (genes_to_phenotype.txt)")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Directorio destino de submissions")
    parser.add_argument("--username", type=str, default="candeleroj98", help="Usuario de Hugging Face")

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    run_track1_pipeline(
        gatk_vcf=args.gatk_vcf,
        manta_vcf=args.manta_vcf,
        mane_gtf=args.mane_gtf,
        cache_dir=args.cache_dir,
        hpo_obo=args.hpo_obo,
        hpo_genes=args.hpo_genes,
        output_dir=args.output_dir,
        hf_username=args.username,
    )


if __name__ == "__main__":
    main()
