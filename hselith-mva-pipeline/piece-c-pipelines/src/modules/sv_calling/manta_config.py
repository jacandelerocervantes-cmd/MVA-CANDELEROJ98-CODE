"""Configuración y validación de parámetros para Illumina Manta v1.6.0."""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class MantaConfig:
    """Parámetros de ejecución para Illumina Manta v1.6.0 WGS Germline.

    Diseñado para el BAM WGS 50x del paciente (55.7 GiB) y el genoma de referencia GRCh38.
    Deshabilita por defecto generate_evidence_bam para optimizar tiempo (~30%) y espacio de disco.
    """

    bam_path: Path
    reference_fasta: Path
    output_dir: Path = Path("manta_output")
    threads: int = 8
    memory_gb: int = 32
    generate_evidence_bam: bool = False
    call_regions_bed: Path | None = None
    manta_install_dir: Path | None = None

    def validate(self, check_files_exist: bool = True) -> None:
        """Valida consistencia de tipos, rangos y existencia física de archivos."""
        if self.threads < 1:
            raise ValueError(f"threads debe ser >= 1, recibido {self.threads}")
        if self.memory_gb < 4:
            raise ValueError(f"memory_gb debe ser >= 4 para WGS 50x, recibido {self.memory_gb}")

        if check_files_exist:
            if not self.bam_path.is_file():
                raise FileNotFoundError(f"Archivo BAM no encontrado: {self.bam_path}")
            # Verificación del índice .bai o .bam.bai
            bai_path_1 = self.bam_path.with_suffix(self.bam_path.suffix + ".bai")
            bai_path_2 = self.bam_path.with_suffix(".bai")
            if not (bai_path_1.is_file() or bai_path_2.is_file()):
                raise FileNotFoundError(
                    f"Índice BAM (.bai) no encontrado para {self.bam_path}. Manta requiere índice."
                )

            if not self.reference_fasta.is_file():
                raise FileNotFoundError(f"FASTA de referencia no encontrado: {self.reference_fasta}")
            # Verificación de índice .fai
            fai_path = self.reference_fasta.with_suffix(self.reference_fasta.suffix + ".fai")
            if not fai_path.is_file():
                raise FileNotFoundError(
                    f"Índice FASTA (.fai) no encontrado para {self.reference_fasta}. Manta requiere samtools faidx."
                )

            if self.call_regions_bed is not None and not self.call_regions_bed.is_file():
                raise FileNotFoundError(f"Archivo BED de regiones no encontrado: {self.call_regions_bed}")

    def resolve_config_script(self) -> str:
        """Resuelve la ruta ejecutable al script configManta.py."""
        if self.manta_install_dir:
            script = self.manta_install_dir / "bin" / "configManta.py"
            if script.is_file():
                return str(script)
        script_in_path = shutil.which("configManta.py")
        if script_in_path:
            return script_in_path
        return "configManta.py"

    def resolve_workflow_script(self) -> Path:
        """Ruta al script de ejecución generado por configManta.py dentro de output_dir."""
        return self.output_dir / "runWorkflow.py"

    def build_config_command(self) -> list[str]:
        """Construye la lista de argumentos CLI para el paso de configuración configManta.py."""
        cmd: list[str] = [
            "python",
            self.resolve_config_script(),
            f"--bam={str(self.bam_path)}",
            f"--referenceFasta={str(self.reference_fasta)}",
            f"--runDir={str(self.output_dir)}",
        ]

        if self.generate_evidence_bam:
            cmd.append("--generateEvidenceBam")

        if self.call_regions_bed is not None:
            cmd.append(f"--callRegions={str(self.call_regions_bed)}")

        return cmd

    def build_run_command(self) -> list[str]:
        """Construye la lista de argumentos CLI para la ejecución del workflow runWorkflow.py."""
        workflow_script = self.resolve_workflow_script()
        return [
            "python",
            str(workflow_script),
            "-m",
            "local",
            "-j",
            str(self.threads),
            "-g",
            str(self.memory_gb),
        ]

