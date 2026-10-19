"""El script de descarga debe usar versiones fijas y verificar SHA-256 sin necesitar red si el archivo ya es correcto."""

from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "download_reference_data.py"
spec = importlib.util.spec_from_file_location("download_reference_data", SCRIPT)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_urls_are_pinned_and_hashes_are_sha256() -> None:
    for name, (url, expected) in module.REFERENCE_FILES.items():
        assert "latest" not in url and "/current/" not in url, f"{name}: la URL debe apuntar a una versión fija"
        assert len(expected) == 64 and int(expected, 16) >= 0


def test_existing_file_with_correct_hash_is_not_downloaded(tmp_path: Path) -> None:
    content = b"contenido de prueba"
    (tmp_path / "ref.txt").write_bytes(content)
    status = module.ensure_file("ref.txt", "http://invalid.invalid/ref.txt", hashlib.sha256(content).hexdigest(), tmp_path)
    assert status == "ok"  # no hubo intento de red
    assert module.sha256_file(tmp_path / "ref.txt") == hashlib.sha256(content).hexdigest()
