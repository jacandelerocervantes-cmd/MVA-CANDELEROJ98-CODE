"""Entrypoint de contenedor y orquestador para SV_CALLING (Manta v1.6.0)."""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path

from src.modules.sv_calling.manta_config import MantaConfig
from src.modules.sv_calling.manta_runner import MantaRunner
from src.modules.sv_calling.sv_vcf_parser import filter_mva_panel_variants, parse_vcf_records

logger = logging.getLogger("hselith.sv_calling")


def execute_sv_calling_pipeline(
    bam_path: Path,
    reference_fasta: Path,
    output_dir: Path,
    threads: int = 8,
    memory_gb: int = 32,
    generate_evidence_bam: bool = False,
    call_regions_bed: Path | None = None,
    validate_inputs: bool = True,
) -> dict[str, int | str | float | list[dict]]:
    """Ejecuta el ciclo de Manta, extrae variantes diploides y clasifica contra el panel MVA."""
    output_dir.mkdir(parents=True, exist_ok=True)

    config = MantaConfig(
        bam_path=bam_path,
        reference_fasta=reference_fasta,
        output_dir=output_dir,
        threads=threads,
        memory_gb=memory_gb,
        generate_evidence_bam=generate_evidence_bam,
        call_regions_bed=call_regions_bed,
    )

    runner = MantaRunner(config)
    exec_result = runner.run(validate_inputs=validate_inputs)

    if not exec_result.success or exec_result.diploid_sv_vcf is None:
        raise RuntimeError(f"Falla en la ejecución de Manta: {exec_result.error_message}")

    # Parseo y clasificación de variantes estructurales
    records = parse_vcf_records(exec_result.diploid_sv_vcf)
    panel_records = filter_mva_panel_variants(records)
    pass_records = [r for r in records if r.get("filter_status") == "PASS"]

    # Guardar reporte de variantes del panel en JSON
    panel_json_path = output_dir / "mva_panel_variants.json"
    with open(panel_json_path, "w", encoding="utf-8") as f:
        json.dump(panel_records, f, indent=2)

    return {
        "total_sv_count": len(records),
        "pass_sv_count": len(pass_records),
        "panel_sv_count": len(panel_records),
        "duration_seconds": exec_result.duration_seconds,
        "diploid_vcf_path": str(exec_result.diploid_sv_vcf),
        "candidate_vcf_path": str(exec_result.candidate_sv_vcf) if exec_result.candidate_sv_vcf else "",
        "panel_json_path": str(panel_json_path),
        "panel_variants": panel_records,
    }


def run_sv_calling_job(pipeline_job_id: str) -> None:
    """Orquesta el job de base de datos para SV_CALLING."""
    from sqlalchemy import select
    from job_orchestration.job_status_updater import mark_job_completed
    from shared.database import get_session
    from shared.models import PipelineJob, SVCallingResult
    from shared.storage import StorageClient

    storage = StorageClient()

    with get_session() as session:
        job = session.execute(
            select(PipelineJob).where(PipelineJob.pipeline_job_id == pipeline_job_id)
        ).scalar_one()
        input_bam_ref = job.input_storage_ref

    local_dir = Path(os.environ.get("MANTA_WORK_DIR", "/tmp/manta_work"))
    local_dir.mkdir(parents=True, exist_ok=True)

    local_bam = local_dir / "patient.bam"
    local_bai = local_dir / "patient.bam.bai"
    storage.download_to_path(input_bam_ref, str(local_bam))
    storage.download_to_path(f"{input_bam_ref}.bai", str(local_bai))

    ref_fasta_env = os.environ.get("REFERENCE_FASTA_PATH", "/data/references/GRCh38.fa")
    ref_fasta = Path(ref_fasta_env)
    threads = int(os.environ.get("MANTA_THREADS", "8"))
    mem_gb = int(os.environ.get("MANTA_MEMORY_GB", "32"))

    summary = execute_sv_calling_pipeline(
        bam_path=local_bam,
        reference_fasta=ref_fasta,
        output_dir=local_dir / "output",
        threads=threads,
        memory_gb=mem_gb,
    )

    # Subida a Cloud Storage
    base_storage_prefix = input_bam_ref.rsplit("/", 1)[0]
    diploid_dest = f"{base_storage_prefix}/sv_calling/{pipeline_job_id}/diploidSV.vcf.gz"
    vcf_storage_ref = storage.upload_from_path(str(summary["diploid_vcf_path"]), diploid_dest)

    cand_storage_ref = None
    if summary["candidate_vcf_path"]:
        cand_dest = f"{base_storage_prefix}/sv_calling/{pipeline_job_id}/candidateSV.vcf.gz"
        cand_storage_ref = storage.upload_from_path(str(summary["candidate_vcf_path"]), cand_dest)

    tool_versions = json.dumps({"manta": "1.6.0", "duration_seconds": summary["duration_seconds"]})

    with get_session() as session:
        result = SVCallingResult(
            pipeline_job_id=pipeline_job_id,
            caller_used="Manta v1.6.0",
            vcf_storage_ref=vcf_storage_ref,
            candidate_vcf_storage_ref=cand_storage_ref,
            total_sv_count=int(summary["total_sv_count"]),
            pass_sv_count=int(summary["pass_sv_count"]),
            panel_sv_count=int(summary["panel_sv_count"]),
            tool_versions=tool_versions,
        )
        session.add(result)
        session.flush()
        sv_result_id = result.sv_calling_id

    mark_job_completed(
        pipeline_job_id,
        candidate_result_ref=f"sv_calling:{sv_result_id}",
        score=float(summary["panel_sv_count"]),
    )
    logger.info("Job SV_CALLING %s finalizado exitosamente", pipeline_job_id)


def main() -> None:
    """Punto de entrada CLI para ejecución local o en Cloud Batch."""
    parser = argparse.ArgumentParser(description="HSelith SV Calling con Illumina Manta v1.6.0")
    parser.add_argument("--job-id", type=str, default=os.environ.get("PIPELINE_JOB_ID"))
    parser.add_argument("--bam", type=str, help="Ruta local al archivo BAM alineado")
    parser.add_argument("--fasta", type=str, help="Ruta local al FASTA de referencia GRCh38")
    parser.add_argument("--output-dir", type=str, default="manta_output")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--memory-gb", type=int, default=32)
    parser.add_argument("--generate-evidence-bam", action="store_true", default=False)
    args = parser.parse_args()

    if args.job_id:
        from job_orchestration.job_status_updater import mark_job_failed
        try:
            run_sv_calling_job(args.job_id)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Fallo en job SV_CALLING %s", args.job_id)
            mark_job_failed(args.job_id, error_message=str(exc))
            sys.exit(1)
    else:
        if not args.bam or not args.fasta:
            parser.error("Se requieren --bam y --fasta si no se especifica --job-id")
        res = execute_sv_calling_pipeline(
            bam_path=Path(args.bam),
            reference_fasta=Path(args.fasta),
            output_dir=Path(args.output_dir),
            threads=args.threads,
            memory_gb=args.memory_gb,
            generate_evidence_bam=args.generate_evidence_bam,
        )
        print(f"SV Calling completado. Variantes totales: {res['total_sv_count']}, Panel MVA: {res['panel_sv_count']}")


if __name__ == "__main__":
    main()

