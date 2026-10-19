"""Descarga los datos públicos de referencia de Track 1 con versión fija y verifica su SHA-256.

Uso (desde la raíz del repositorio):
    python hselith-mva-pipeline/piece-c-pipelines/scripts/download_reference_data.py [--dest data/reference]

Un archivo ya presente con el hash correcto no se vuelve a descargar. Un hash distinto aborta el
proceso: significa que la fuente cambió y los resultados dejarían de ser reproducibles.
No descarga datos de pacientes (los VCF del reto se obtienen del dataset oficial).
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.request
from pathlib import Path

REFERENCE_FILES: dict[str, tuple[str, str]] = {
    "MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz": (
        "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz",
        "71d4b3c89c9d948683bf0db5d81b00d9ae9f3b943177b74cbd87e68f08e34d66",
    ),
    "hp.obo": (
        "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/hp.obo",
        "93dace952fcb3ec4728818857f6ba76bc2d5312f4d83266519b8694b8e798f22",
    ),
    "genes_to_phenotype.txt": (
        "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/genes_to_phenotype.txt",
        "507a17bff9c49e6329fbd88b1f91734fa7e9fdea5c06304044c0892eb6ab248c",
    ),
}


def sha256_file(path: Path) -> str:
    """SHA-256 de un archivo, leído en bloques."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def ensure_file(name: str, url: str, expected: str, dest: Path) -> str:
    """Devuelve 'ok' (ya estaba) o 'descargado'; lanza ValueError si el hash no coincide."""
    target = dest / name
    if target.exists() and sha256_file(target) == expected:
        return "ok"
    dest.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    with urllib.request.urlopen(url, timeout=300) as response, open(partial, "wb") as out:  # noqa: S310 - URL fija
        while chunk := response.read(1 << 20):
            out.write(chunk)
    found = sha256_file(partial)
    if found != expected:
        partial.unlink()
        raise ValueError(f"{name}: SHA-256 {found} distinto del esperado {expected}; la fuente cambió.")
    partial.replace(target)
    return "descargado"


def main() -> int:
    parser = argparse.ArgumentParser(description="Descarga y verifica los datos públicos de referencia de Track 1.")
    parser.add_argument("--dest", type=Path, default=Path("data/reference"), help="Directorio destino")
    args = parser.parse_args()
    for name, (url, expected) in REFERENCE_FILES.items():
        print(f"{name}: {ensure_file(name, url, expected, args.dest)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
