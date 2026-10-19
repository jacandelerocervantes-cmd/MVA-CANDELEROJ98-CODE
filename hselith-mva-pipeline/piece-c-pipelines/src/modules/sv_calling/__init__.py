"""Módulo ejecutor de SV Calling con Illumina Manta v1.6.0, búsqueda genome-wide y clasificación de panel MVA."""

from __future__ import annotations

from src.modules.sv_calling.call_sv import execute_sv_calling_pipeline
from src.modules.sv_calling.manta_config import MantaConfig
from src.modules.sv_calling.manta_runner import MantaExecutionResult, MantaRunner
from src.modules.sv_calling.compound_het_detector import (
    CandidatePair,
    detect_secondary_structural_architectures,
)
from src.modules.sv_calling.consequence_calculator import classify_allele
from src.modules.sv_calling.dual_vcf_parser import (
    parse_dual_vcf_record,
    parse_dual_vcfs,
    stream_vcf_records,
)
from src.modules.sv_calling.epcr_ranking_engine import (
    order_epcr_rows,
    generate_submission_files,
    validate_submission_rows,
)
from src.modules.sv_calling.sv_vcf_parser import (
    MVA_22_GENES,
    GeneRegion,
    filter_mva_panel_variants,
    parse_vcf_record,
    parse_vcf_records,
)

__all__ = [
    "MantaConfig",
    "MantaRunner",
    "MantaExecutionResult",
    "GeneRegion",
    "MVA_22_GENES",
    "parse_vcf_record",
    "parse_vcf_records",
    "filter_mva_panel_variants",
    "execute_sv_calling_pipeline",
    "parse_dual_vcf_record",
    "parse_dual_vcfs",
    "stream_vcf_records",
    "classify_allele",
    "CandidatePair",
    "detect_secondary_structural_architectures",
    "order_epcr_rows",
    "generate_submission_files",
    "validate_submission_rows",
]

