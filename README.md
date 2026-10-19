# HSelith — Track 1 (Variant Prediction), MVA Hackathon 2026

Deterministic pipeline that searches a whole-genome sequencing VCF for recessive candidate variants in a proband with
Mosaic Variegated Aneuploidy (MVA), computes every functional consequence locally, and writes the 10-row submission file.
This repository contains **Track 1 only**: the code, the tests, the submission file, the written report and the methods description. It does **not**
contain patient data (reads, VCFs, alignments), the patient's clinical terms, the supporting patient-derived outputs (genome-wide ranking and support JSONs) or Track 2 material.

| | |
|---|---|
| Hugging Face user | `candeleroj98` |
| Proband | `PROBAND01` (sample identifier omitted), GRCh38 |
| Primary prediction | `BUB1B` compound-heterozygous pair p.Leu737\* + p.Asn1002Lys (positions in the submission CSV); phase undetermined |
| Submission file | [`submissions/track1/candeleroj98_gatk-manta-fusion.csv`](submissions/track1/candeleroj98_gatk-manta-fusion.csv) |
| Report / methods | [`report (md)`](submissions/track1/candeleroj98_track1_report.md), [`methods (md)`](submissions/track1/candeleroj98_methods_description.md), [`methods form (xlsx)`](submissions/track1/candeleroj98_methods_description_form.xlsx) |
| Known limits | section 4.1 of the report |

## What the pipeline does

1. Scans every PASS SNV/indel of the GATK VCF that falls in a coding region or within 2 bp of a splice junction of the 19,293 MANE Select transcripts (MANE v1.5, GRCh38).
2. Computes the Sequence Ontology consequence of each allele from the transcript model and the reference coding sequence (Ensembl REST). Ensembl VEP is used only as an external check.
3. Queries gnomAD v4.1 allele by allele (public GraphQL API) to keep rare variants (absent or AF < 1e-3), pairs them into recessive candidates (compound heterozygous or homozygous/hemizygous) and ranks them.
4. Scores the candidates against the proband's 8 HPO terms (HPO release 2026-09-01) as a consistency check.
5. Adds structural-variant rows (Manta) inside a curated 22-gene MVA panel and writes the CSV (primary row + secondary rows, at most 10).

Every number in the CSV is computed from the inputs; no gene, position or annotation table decides the primary row (the code holds only the structural-variant panel, descriptive disease notes that enter no score, and the list of excluded hypervariable gene families). The `epcr` column is a **ranking value, not a calibrated probability**.

## Repository layout

```
.
├── README.md, LICENSE, THIRD_PARTY_LICENSES.md
├── pytest.ini, conftest.py                    test configuration (69 tests pass in a clean clone; 1 needs the local proband HPO terms and is skipped)
├── requirements-track1.txt / .lock            dependencies (lock = exact versions used)
├── run_track1_verification.sh / .ps1          runs the test suite
├── submissions/track1/                        submitted CSV, report and methods (supporting patient-derived outputs are generated locally and git-ignored)
└── hselith-mva-pipeline/
    ├── piece-c-pipelines/
    │   ├── src/modules/sv_calling/            the pipeline (generate_track1_submission.py is the entry point)
    │   ├── scripts/download_reference_data.py pinned download of MANE and HPO with SHA-256 checks
    │   └── tests/
    └── shared-contracts/schemas/              JSON Schema (Draft 2020-12) for the submission rows
```

## Requirements

- **Python 3.10 or newer** (tested with 3.10.11) and **Git**. No compiled dependencies for the main pipeline.
- Internet access on the first run (Ensembl REST, gnomAD API, NCBI and GitHub downloads). Later runs can be offline once the cache exists.
- About 2 GB of free disk for reference data and caches (more if you also reproduce the alignment, see Path C).

## Tested environment and versions

Everything below was executed by the author on **Windows 11 (10.0.26200), Python 3.10.11, pip 23.0.1, setuptools 78.1.0, Git 2.55**, from a clean copy of this repository (all 70 tests passed with the local patient-derived files present; the pipeline regenerated the four output files byte for byte). The code contains no operating-system-specific logic (paths use `pathlib`, files are written with explicit LF), so the Linux and macOS commands follow the same code path, but the author did **not** execute them on those systems.

**Python packages** (exact pins in [`requirements-track1.lock`](requirements-track1.lock); the looser `requirements-track1.txt` accepts `jsonschema>=4.20,<5` and `pytest>=8,<10`):

| Package | Version | Role |
|---|---|---|
| pytest | 9.1.1 | test runner (pulls pluggy 1.6.0, iniconfig 2.3.0, packaging 26.3, pygments 2.21.0; colorama 0.4.6 on Windows only; exceptiongroup 1.3.1 and tomli 2.4.1 on Python < 3.11) |
| jsonschema | 4.26.0 | validates every submission row against JSON Schema Draft 2020-12 (pulls attrs 26.1.0, jsonschema-specifications 2025.9.1, referencing 0.37.0, rpds-py 0.30.0, typing_extensions 4.16.0) |

All of them are pure-Python wheels except `rpds-py`, a Rust extension that publishes prebuilt wheels for CPython 3.10+ on Windows, macOS (x86_64 and arm64) and Linux (manylinux): no compiler is needed. Apart from `jsonschema` (schema validation of the output rows), the pipeline uses only the Python standard library. Check what pip will install with `pip install --dry-run -r requirements-track1.lock` and inspect the result with `pip list` or `pip freeze`. The interpreter version is also recorded in `.python-version`.

**Data and services** (all versions pinned or recorded):

| Component | Version | Used for |
|---|---|---|
| Reference genome | GRCh38 `GCA_000001405.15_GRCh38_no_alt_analysis_set` (with `chrEBV`) | alignment for Manta; the challenge VCF header names the `plus_hs38d1_maskedGRC_exclusions_v2` variant |
| MANE Select | v1.5 (`MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz`, NCBI release_1.5), SHA-256 verified | transcript models |
| Ensembl REST | release 116 (GRCh38.p14) | reference coding sequences; VEP terms for validation only |
| gnomAD | v4.1 (`gnomad_r4` dataset, public GraphQL API, queried 2026-10-05) | allele frequencies (rarity) |
| Human Phenotype Ontology | release 2026-09-01 (`hp.obo`, `genes_to_phenotype.txt`), SHA-256 verified | phenotype-similarity score |
| SNV/indel calls (challenge VCF) | Sentieon Haplotyper 202308.02 + GATK VariantFiltration (from the VCF header) | input variants |
| Structural variants | Illumina Manta 1.6.0 | input structural variants |
| Alignment (Path C only) | BWA-MEM2 2.2.1, samtools, Python 2 for Manta's `runWorkflow.py` | reproducing the Manta VCF |

## Path A — Verify the code (no data needed, about 1 minute)

Works the same on all three systems; only the shell commands differ.

### Linux and macOS (bash or zsh)

```bash
git clone https://github.com/jacandelerocervantes-cmd/MVA-CANDELEROJ98-CODE.git hselith-track1
cd hselith-track1

python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements-track1.lock      # or requirements-track1.txt for the latest compatible versions

python -m pytest -v                          # 69 passed, 1 skipped (the skipped test needs the local proband HPO terms)
```

macOS: if `python3` is missing, install it with Homebrew (`brew install python`) or from python.org.

### Windows (PowerShell)

```powershell
git clone https://github.com/jacandelerocervantes-cmd/MVA-CANDELEROJ98-CODE.git hselith-track1
cd hselith-track1

py -3 -m venv .venv
# If activation is blocked by the execution policy, run this once for the current window only:
#   Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements-track1.lock

python -m pytest -v                          # 69 passed, 1 skipped (the skipped test needs the local proband HPO terms)
```

Git Bash on Windows also works with the Linux commands above, activating with `source .venv/Scripts/activate`.

The helper scripts `./run_track1_verification.sh` (Linux/macOS/Git Bash) and `.\run_track1_verification.ps1` (PowerShell) run the same test suite with the interpreter on your `PATH` (override with the `PYTHON` environment variable).

## Path B — Regenerate the submission from the VCFs

You need two input VCFs and the public reference data. Use the layout below (the pipeline takes explicit paths, so any layout works).

```
data/
├── dataset_mva_hackathon/<challenge_snv_indel_vcf>.vcf.gz     SNV/indel VCF from the challenge dataset
├── manta_results/diploidSV.vcf.gz                            Manta v1.6.0 structural variants (see Path C)
└── reference/                                                created by the download script below
```

`data/` is ignored by git on purpose: the VCFs come from the challenge dataset under its data use agreement and must not be republished.

### 1. Download the public reference data (pinned versions, SHA-256 verified)

Linux / macOS:

```bash
python hselith-mva-pipeline/piece-c-pipelines/scripts/download_reference_data.py --dest data/reference
```

Windows (PowerShell):

```powershell
python hselith-mva-pipeline\piece-c-pipelines\scripts\download_reference_data.py --dest data\reference
```

This fetches MANE v1.5 (`MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz`) and HPO release 2026-09-01 (`hp.obo`, `genes_to_phenotype.txt`). A hash mismatch aborts the download: it means the source changed and results would no longer be reproducible.

### 2. Get the challenge SNV/indel VCF

Access to the dataset is governed by the challenge terms. After accepting them:

```bash
pip install huggingface_hub
hf download SageBio/mva-hackathon-2026-data --repo-type dataset \
    --include "*.vcf.gz" "*.vcf.gz.tbi" --local-dir data/dataset_mva_hackathon
```

(The full dataset also holds about 80 GB of FASTQ files; the `--include` filter downloads only the VCF. Check the file name against the checksum table below. On older `huggingface_hub` versions the command is `huggingface-cli download`.)

**Proband HPO terms.** The patient's clinical terms are not published. Copy `hselith-mva-pipeline/piece-c-pipelines/src/modules/sv_calling/reference_data/proband_hpo_terms.example.json` to `proband_hpo_terms.json` in the same folder (git-ignored) and list the proband's HPO terms from the Hackathon clinical document.

### 3. Run the generator

Linux / macOS:

```bash
python hselith-mva-pipeline/piece-c-pipelines/src/modules/sv_calling/generate_track1_submission.py \
  --gatk-vcf data/dataset_mva_hackathon/<challenge_snv_indel_vcf>.vcf.gz \
  --manta-vcf data/manta_results/diploidSV.vcf.gz \
  --mane-gtf data/reference/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
  --hpo-obo data/reference/hp.obo \
  --hpo-genes data/reference/genes_to_phenotype.txt \
  --output-dir out
```

Windows (PowerShell; the backtick continues the line):

```powershell
python hselith-mva-pipeline\piece-c-pipelines\src\modules\sv_calling\generate_track1_submission.py `
  --gatk-vcf data\dataset_mva_hackathon\<challenge_snv_indel_vcf>.vcf.gz `
  --manta-vcf data\manta_results\diploidSV.vcf.gz `
  --mane-gtf data\reference\MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz `
  --hpo-obo data\reference\hp.obo `
  --hpo-genes data\reference\genes_to_phenotype.txt `
  --output-dir out
```

Writing to `out/` keeps the submitted file in `submissions/track1/` untouched so you can compare. The generator has no default input: without real VCF paths it stops with an error instead of producing a placeholder.

**Timing.** The first run downloads reference coding sequences (Ensembl), gnomAD v4 frequencies and VEP validation terms into `data/reference/cache/`: about **60 minutes** (63 min measured on 2026-10-05; the gnomAD service rate-limits and may answer "Service overloaded", which the code retries and never reads as "variant absent"). With the cache present the run takes **under a minute** and needs no network. The cache holds the proband's variant list, so it is not published; keep it local.

### 4. Check the result

Compare your output with the submitted file (the byte-for-byte match holds for the VCF → CSV step):

```bash
# Linux
sha256sum out/candeleroj98_gatk-manta-fusion.csv submissions/track1/candeleroj98_gatk-manta-fusion.csv
# macOS
shasum -a 256 out/candeleroj98_gatk-manta-fusion.csv submissions/track1/candeleroj98_gatk-manta-fusion.csv
```

```powershell
# Windows
Get-FileHash -Algorithm SHA256 out\candeleroj98_gatk-manta-fusion.csv, submissions\track1\candeleroj98_gatk-manta-fusion.csv
```

Expected: identical hashes, identical to the table in the next section. Line endings are LF on every system (`.gitattributes` enforces it).

## Path C — Reproduce the Manta VCF from the reads (optional, heavy)

The structural-variant VCF is generated from the proband's reads and cannot be redistributed. Re-running alignment and calling is **not guaranteed to reproduce it bit for bit**; the checksums below apply to the VCF → CSV step only.

Requirements: Linux with at least 12 vCPUs, 32 GB RAM and 350 GB of SSD (128 GB RAM to build the BWA index; measured 25 min on 12 vCPUs). **Manta and BWA-MEM2 do not run natively on Windows or macOS**: on Windows use WSL2 (Ubuntu 22.04) or Docker Desktop; on macOS use a Linux VM or Docker (`--platform linux/amd64`, much slower on Apple Silicon). The commands below run inside Linux.

```bash
# tools
sudo apt-get update && sudo apt-get install -y build-essential zlib1g-dev libbz2-dev liblzma-dev samtools bzip2 python3-pip
wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
tar -xjf bwa-mem2-2.2.1_x64-linux.tar.bz2 && sudo cp bwa-mem2-2.2.1_x64-linux/bwa-mem2* /usr/local/bin/
wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2
tar -xjf manta-1.6.0.centos6_x86_64.tar.bz2 && sudo mv manta-1.6.0.centos6_x86_64 /opt/manta

# reads (about 80 GB; accept the challenge terms first)
hf download SageBio/mva-hackathon-2026-data --repo-type dataset --local-dir data/fastq

# reference GRCh38 no-alt analysis set, index
mkdir -p data/reference data/alignment data/sv_calling && cd data/reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
bwa-mem2 index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# alignment (about 14 h on 12 vCPUs)
cd ../alignment
bwa-mem2 mem -t 12 -R '@RG\tID:PROBAND01\tSM:PROBAND01\tPL:ILLUMINA' \
  ../reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  <(zcat ../fastq/*_R1_*.fastq.gz) <(zcat ../fastq/*_R2_*.fastq.gz) | \
  samtools sort -@ 8 -m 2G -o patient_aligned.sorted.bam
samtools index -@ 8 patient_aligned.sorted.bam

# structural variants (about 1.5 h)
cd ../sv_calling
/opt/manta/bin/configManta.py --bam ../alignment/patient_aligned.sorted.bam \
  --referenceFasta ../reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --runDir ./manta_workspace
python2 ./manta_workspace/runWorkflow.py -m local -j 12 -g 30
# result: data/sv_calling/manta_workspace/results/variants/diploidSV.vcf.gz  -> use it as --manta-vcf
```

## Checksums (VCF → CSV step)

Run `sha256sum` / `shasum -a 256` / `Get-FileHash` on each file. Four output files are listed with their hashes; only the submission CSV is versioned, and the other three (patient-derived supporting outputs) are not published but their hashes let you verify a regeneration. The rest are inputs and local caches.

| File | SHA-256 |
|---|---|
| Challenge SNV/indel VCF (GATK) | `bec64c68d3ad10d05b8071d07cdf87d217013d66d3ba02e24b505035c9d42175` |
| `diploidSV.vcf.gz` (Manta v1.6.0) | `bb12959d2ebd99496cb6b353e53801e64008b00ee568a185331e2d1ef70e8db1` |
| `submissions/track1/candeleroj98_gatk-manta-fusion.csv` (output) | `264a0a43f75997ca04603a3c6292cc8aa8ef305e16297ec90fdf5a96fbf51c92` |
| `candeleroj98_genomewide_ranking.csv` (output, not published) | `0da25bcd3c01d0401270229234a4e3e7f1242af28f6dfc2b7c48bfa8bf8b606b` |
| `candeleroj98_phenotype_summary.json` (output, not published) | `82485593d9b5d1464cb5aff8a7e2a9a07604f5e657a2a12de06f0586a1b33511` |
| `candeleroj98_consequence_validation.json` (output, not published) | `4880ac10e1b08d02abc581ec96187c4c26a0f5e9204e4a362951cb06ea1d4821` |
| `MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz` (NCBI MANE v1.5, input) | `71d4b3c89c9d948683bf0db5d81b00d9ae9f3b943177b74cbd87e68f08e34d66` |
| `hp.obo` (HPO release 2026-09-01, input) | `93dace952fcb3ec4728818857f6ba76bc2d5312f4d83266519b8694b8e798f22` |
| `genes_to_phenotype.txt` (HPO, input) | `507a17bff9c49e6329fbd88b1f91734fa7e9fdea5c06304044c0892eb6ab248c` |
| `data/reference/cache/cds_sequences.json` (Ensembl reference CDS cache) | `bf7563dd4b53431adc18ec78ff1a1ca36cc5f7b394c0dd6f57d29702dedd8b5c` |
| `data/reference/cache/gnomad_v4_frequencies.json` (gnomAD v4 frequency cache) | `0b5b1e42c9357df768d492ea4800886c964911fb0db88e79f8e382ef2791b7f0` |
| `data/reference/cache/vep_terms.json` (Ensembl VEP consequence terms, validation only) | `377fd7deadd7a6d99b85cc9e4360bcb0bdbe40b53d3a8356edfaf0c7d93442f7` |
| `src/modules/sv_calling/reference_data/proband_hpo_terms.json` (proband HPO terms; not published) | `a77a6014003c23601a71c78104f6a27118dc97ef579e3f13589389a28eb26ada` |

SV sizes in the CSV notes are Manta `SVLEN` (they exclude the VCF anchor base). Without the caches, the Ensembl and gnomAD versions served at run time may differ from those used here.

## Troubleshooting

| Symptom | Likely cause and fix |
|---|---|
| `ModuleNotFoundError: jsonschema` | The virtual environment is not active or the requirements are not installed (`pip install -r requirements-track1.lock`). |
| `...ps1 cannot be loaded because running scripts is disabled` (Windows) | `Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass`, then activate again. |
| Download script stops with a hash error | The remote file changed; do not continue, results would not be reproducible. |
| gnomAD query fails or says "Service overloaded" | The code retries; if it still fails, rerun later (the cache keeps what was already fetched). A failed query is never treated as "variant absent". |
| The CSV hash differs by line endings | Make sure git did not rewrite them: `git config core.autocrlf false`, then check the files out again (`.gitattributes` already forces LF). |
| `Manta` or `bwa-mem2` not found on Windows/macOS | They are Linux-only; use WSL2, Docker or a Linux VM (Path C). |

## Human Phenotype Ontology

This service/product uses the Human Phenotype Ontology (release 2026-09-01). Find out more at http://www.human-phenotype-ontology.org. The HPO files are used unaltered, under the HPO licence (https://hpo.jax.org/license).

## Licence

As the Hackathon rules require, the submission (report, tables, other submitted files and the source code) is released under CC BY 4.0 (file `LICENSE`; https://creativecommons.org/licenses/by/4.0/). Third-party resources keep their own licences (see `THIRD_PARTY_LICENSES.md`); where a resource is CC BY-SA, derived tables that reproduce its content keep that licence.

## Data privacy and compliance

Raw reads, VCFs, alignments and caches derived from the proband are covered by the challenge's data use agreement. They are ignored by git (`data/`, `*.bam`, `*.fastq.gz`, `*.vcf.gz`) and must be deleted by the deadline set in the agreement. The only patient-derived output versioned here is the submission CSV (a derived prediction table, which the Hackathon rules allow to be shared). The genome-wide ranking, the support JSONs and the proband's HPO terms are generated or supplied locally and are git-ignored (see `.gitignore`).

## Acknowledgments

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.

We also thank the Wilhelm Foundation for their guidance and insight on the project. No information that could re-identify the patient or family is included in this repository.

## Generative-AI disclosure

See section 5.4 of the report and question 4 of [`submissions/track1/candeleroj98_methods_description.md`](submissions/track1/candeleroj98_methods_description.md).
