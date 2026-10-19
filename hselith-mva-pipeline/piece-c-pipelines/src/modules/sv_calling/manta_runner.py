"""Runner y orquestador de ejecución para Illumina Manta v1.6.0."""

from __future__ import annotations

import logging
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

from src.modules.sv_calling.manta_config import MantaConfig

logger = logging.getLogger("hselith.sv_calling.manta")


@dataclass(frozen=True)
class MantaExecutionResult:
    """Resultado estructurado de la corrida de Manta v1.6.0."""

    success: bool
    output_dir: Path
    diploid_sv_vcf: Path | None
    candidate_sv_vcf: Path | None
    candidate_indels_vcf: Path | None
    duration_seconds: float
    stdout: str
    stderr: str
    error_message: str | None = None


class MantaRunner:
    """Orquestador de subprocessos para Manta v1.6.0."""

    def __init__(self, config: MantaConfig) -> None:
        self.config = config

    def run(self, validate_inputs: bool = True) -> MantaExecutionResult:
        """Ejecuta el ciclo de dos etapas de Manta: configManta.py -> runWorkflow.py."""
        start_time = time.perf_counter()
        logger.info("Iniciando ejecución de Illumina Manta v1.6.0")

        if validate_inputs:
            try:
                self.config.validate(check_files_exist=True)
            except Exception as e:
                logger.error("Falla de pre-requisitos para Manta: %s", e)
                return MantaExecutionResult(
                    success=False,
                    output_dir=self.config.output_dir,
                    diploid_sv_vcf=None,
                    candidate_sv_vcf=None,
                    candidate_indels_vcf=None,
                    duration_seconds=time.perf_counter() - start_time,
                    stdout="",
                    stderr="",
                    error_message=str(e),
                )

        # 1. Fase de configuración
        config_cmd = self.config.build_config_command()
        logger.info("Ejecutando configManta: %s", " ".join(config_cmd))
        try:
            config_proc = subprocess.run(
                config_cmd,
                capture_output=True,
                text=True,
                check=False,
            )
        except Exception as e:
            return MantaExecutionResult(
                success=False,
                output_dir=self.config.output_dir,
                diploid_sv_vcf=None,
                candidate_sv_vcf=None,
                candidate_indels_vcf=None,
                duration_seconds=time.perf_counter() - start_time,
                stdout="",
                stderr="",
                error_message=f"Error invocando configManta.py: {e}",
            )

        if config_proc.returncode != 0:
            err = f"configManta.py falló con código {config_proc.returncode}: {config_proc.stderr}"
            logger.error(err)
            return MantaExecutionResult(
                success=False,
                output_dir=self.config.output_dir,
                diploid_sv_vcf=None,
                candidate_sv_vcf=None,
                candidate_indels_vcf=None,
                duration_seconds=time.perf_counter() - start_time,
                stdout=config_proc.stdout,
                stderr=config_proc.stderr,
                error_message=err,
            )

        workflow_script = self.config.resolve_workflow_script()
        if not workflow_script.is_file():
            err = f"Script de workflow no generado por configManta: {workflow_script}"
            logger.error(err)
            return MantaExecutionResult(
                success=False,
                output_dir=self.config.output_dir,
                diploid_sv_vcf=None,
                candidate_sv_vcf=None,
                candidate_indels_vcf=None,
                duration_seconds=time.perf_counter() - start_time,
                stdout=config_proc.stdout,
                stderr=config_proc.stderr,
                error_message=err,
            )

        # 2. Fase de ejecución del workflow
        run_cmd = self.config.build_run_command()
        logger.info("Ejecutando runWorkflow: %s", " ".join(run_cmd))
        try:
            run_proc = subprocess.run(
                run_cmd,
                capture_output=True,
                text=True,
                check=False,
            )
        except Exception as e:
            return MantaExecutionResult(
                success=False,
                output_dir=self.config.output_dir,
                diploid_sv_vcf=None,
                candidate_sv_vcf=None,
                candidate_indels_vcf=None,
                duration_seconds=time.perf_counter() - start_time,
                stdout=config_proc.stdout,
                stderr="",
                error_message=f"Error invocando runWorkflow.py: {e}",
            )

        total_duration = time.perf_counter() - start_time
        combined_stdout = f"--- CONFIG STDOUT ---\n{config_proc.stdout}\n--- RUN STDOUT ---\n{run_proc.stdout}"
        combined_stderr = f"--- CONFIG STDERR ---\n{config_proc.stderr}\n--- RUN STDERR ---\n{run_proc.stderr}"

        if run_proc.returncode != 0:
            err = f"runWorkflow.py falló con código {run_proc.returncode}: {run_proc.stderr}"
            logger.error(err)
            return MantaExecutionResult(
                success=False,
                output_dir=self.config.output_dir,
                diploid_sv_vcf=None,
                candidate_sv_vcf=None,
                candidate_indels_vcf=None,
                duration_seconds=total_duration,
                stdout=combined_stdout,
                stderr=combined_stderr,
                error_message=err,
            )

        # 3. Detección de archivos de salida esperados
        variants_dir = self.config.output_dir / "results" / "variants"
        diploid_vcf = variants_dir / "diploidSV.vcf.gz"
        candidate_vcf = variants_dir / "candidateSV.vcf.gz"
        indels_vcf = variants_dir / "candidateSmallIndels.vcf.gz"

        diploid_found = diploid_vcf if diploid_vcf.is_file() else None
        candidate_found = candidate_vcf if candidate_vcf.is_file() else None
        indels_found = indels_vcf if indels_vcf.is_file() else None

        if not diploid_found:
            err = f"Manta finalizó exitosamente pero no se encontró {diploid_vcf}"
            logger.warning(err)
            return MantaExecutionResult(
                success=False,
                output_dir=self.config.output_dir,
                diploid_sv_vcf=None,
                candidate_sv_vcf=candidate_found,
                candidate_indels_vcf=indels_found,
                duration_seconds=total_duration,
                stdout=combined_stdout,
                stderr=combined_stderr,
                error_message=err,
            )

        logger.info("Ejecución de Manta completada con éxito en %.2f segundos", total_duration)
        return MantaExecutionResult(
            success=True,
            output_dir=self.config.output_dir,
            diploid_sv_vcf=diploid_found,
            candidate_sv_vcf=candidate_found,
            candidate_indels_vcf=indels_found,
            duration_seconds=total_duration,
            stdout=combined_stdout,
            stderr=combined_stderr,
            error_message=None,
        )

