"""Root conftest.py for HSelith-MVA-Genomics repository.

Acts as an architectural guardrail preventing naked `pytest` executions
in unconfigured interpreters from crashing with unhelpful tracebacks.
Uses strictly ASCII text to avoid encoding corruption in Windows consoles.
"""

from __future__ import annotations

import sys
import pytest


def pytest_configure(config: pytest.Config) -> None:
    """Validate runtime environment dependencies before running tests."""
    missing_deps: list[str] = []

    try:
        import jsonschema  # noqa: F401
    except ImportError:
        missing_deps.append("jsonschema")

    if missing_deps:
        msg = (
            "\n"
            + "=" * 80
            + "\n"
            + "[ERROR DE ENTORNO EN HSELITH-MVA-GENOMICS]\n"
            + "pytest se esta ejecutando en un interprete Python sin las dependencias requeridas:\n"
            + f"  Interprete actual: {sys.executable}\n"
            + f"  Paquetes faltantes: {', '.join(missing_deps)}\n\n"
            + "Para ejecutar la verificacion con el entorno virtual preconfigurado del repositorio:\n"
            + "  - En Windows PowerShell:   .\\run_track1_verification.ps1\n"
            + "  - En Linux / macOS / Bash: ./run_track1_verification.sh\n\n"
            + "O bien ejecuta pytest directamente usando el interprete del entorno virtual:\n"
            + "  & .\\hselith-mva-pipeline\\.venv_x64\\Scripts\\python.exe -m pytest\n\n"
            + "Si deseas utilizar tu interprete Python actual, instala las dependencias con:\n"
            + "  pip install -r requirements-track1.txt\n"
            + "=" * 80
        )
        pytest.exit(msg, returncode=1)
