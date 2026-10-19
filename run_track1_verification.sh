#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"
PYTHON="${PYTHON:-python}"
# La configuración de pytest.ini (raíz) ya limita las pruebas a las de Track 1 y los esquemas compartidos.
"$PYTHON" -m pytest -v
