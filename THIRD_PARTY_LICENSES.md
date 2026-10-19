# Third-party data and software

Third-party data are not redistributed in this repository (`data/` is ignored); they are downloaded or queried at run time. The licences below are **not verified in this session**: check the terms of each source before reusing it.

| Resource | Version used | Terms |
|---|---|---|
| MANE Select (NCBI and EMBL-EBI) | v1.5, GRCh38 | see the NCBI MANE page; cite Morales et al., Nature 2022 |
| Ensembl (REST service, coding sequences, VEP as an external control) | release 116 | EMBL-EBI terms of use: Ensembl states that it imposes no restrictions on access to or use of its data; some third-party material may carry its own constraints |
| gnomAD (public GraphQL API) | v4.1 (gnomad_r4), queried 2026-10-05 | see the terms on the gnomAD downloads page; cite Chen et al., Nature 2024 |
| Human Phenotype Ontology | release 2026-09-01 | see the HPO licence page (https://hpo.jax.org/app/license); acknowledge and cite the HPO consortium |
| ClinVar, OMIM, UniProt | queried for annotation notes | see each source; OMIM numbers are cited only as identifiers |
| Challenge dataset (SageBio/mva-hackathon-2026-data) | gated | CC BY 4.0 licence file shipped with the dataset, plus the challenge data use agreement; not redistributed |

Python dependencies: `jsonschema` and `pytest` (see `requirements-track1.lock`); `markdown` is optional, for the PDF of the report.
