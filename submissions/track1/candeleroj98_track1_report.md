# Rare-variant pairing and phenotype ranking place two heterozygous *BUB1B* variants first in the patient, with phase undetermined

**MVA Hackathon 2026 · Track 1 – Variant Prediction** · Juan Antonio Candelero Cervantes (single author) · Hugging Face model `candeleroj98_gatk-manta-fusion` · Code: https://github.com/jacandelerocervantes-cmd/MVA-CANDELEROJ98-CODE-TRACK1 · Licence: CC BY 4.0

---

## Summary

The patient carries two rare changes in a gene called *BUB1B*: one is expected to cut the gene's protein short, and the other swaps a single building block of that protein. Changes in both copies of this gene are linked to mosaic variegated aneuploidy, a rare condition in which some of the body's cells carry the wrong number of chromosomes. The first change is already listed in a public medical database as harmful, but the effect of the second change is unknown, and it is not yet known whether the two changes sit on different copies of the gene. *BUB1B* is therefore the clearest lead from this analysis and a starting point for further research into the patient's condition.

---

## 1 Methods

### 1.1 Data

- Reference genome: GRCh38.
- Sequencing reads from the patient.
- A file of single-nucleotide variants and small insertions and deletions (indels), produced with Sentieon Haplotyper 202308.02 [1] and GATK VariantFiltration [2].
- Eight Human Phenotype Ontology (HPO) [3] terms describing the patient's phenotype.

### 1.2 Target transcripts and variant annotation

The analysis was restricted to the 19,293 MANE Select transcripts (Matched Annotation from NCBI and EMBL-EBI, v1.5) [4]. Of 4,740,790 variants that passed the caller's filters, 22,995 lay in coding sequence or within 2 bp of an exon–intron boundary. Records with two alternate alleles in the genotype were split into one allele each, giving 23,035 unique alleles. Protein consequences were computed with custom code, using Sequence Ontology terms [5] and reference coding sequences from Ensembl release 116 [6]. The Ensembl Variant Effect Predictor (VEP) [7] was used only as an independent comparison: agreement was 98.73% for impact class and 98.68% for consequence term, with no reference-allele mismatches. Of the 23,035 alleles, 11,364 were classed as high or moderate impact.

### 1.3 Population-frequency filter

Allele frequencies were taken from the Genome Aggregation Database (gnomAD) v4.1 [8] (807,162 individuals, released 19 April 2024; public API, dataset identifier gnomad_r4; queried on 2026-10-05), using the higher of the exome and genome frequencies. A variant was classed as rare if it was absent from gnomAD or had an allele frequency below 1 in 1,000, a fixed threshold; Whiffin et al. [9] describe a stricter, disease-specific alternative based on the maximum credible allele frequency.

### 1.4 Recessive model and variant pairing

Rare high- and moderate-impact variants were tested against a recessive model: two variants in the same gene, a homozygous variant, or a hemizygous variant on the X chromosome outside the pseudoautosomal regions (criteria in Section 3). Seven highly polymorphic gene families (HLA, KIR, OR, MUC, KRTAP, PRB, TAS2R) were excluded. After the frequency filter and the gene-family exclusion, 209 rare variants in 190 genes remained. Within these, pairs of heterozygous variants less than 150 bp apart were discarded because they probably lie on the same copy of the gene (in cis), and heterozygous calls on the non-pseudoautosomal X chromosome were discarded (Section 3). Five genes met the recessive model.

### 1.5 Candidate ranking

The variant-based rank orders candidates by the number of high-impact variants and then by joint rarity. Joint rarity is −log₁₀ of the probability of carrying the observed genotype, computed from population allele frequencies with a frequency floor of 1 in 1,000,000 (formula in Section 3). Phenotype similarity is the coverage of the patient's eight HPO terms by each gene's HPO annotations, weighted by the information content of each term. The final rank orders candidates by phenotype similarity, with the number of high-impact variants and then joint rarity as tie-breakers.

### 1.6 Structural variant filtering

Reads were re-aligned to GRCh38 with BWA-MEM2 2.2.1 [10] (99.57% of reads mapped; mean autosomal depth 51.66×), and structural variants (SVs) were called with Manta 1.6.0 [11]. Calls were retained if they passed Manta's filters, were at least 50 bp long, had fully resolved sequence and overlapped a panel of 22 genes associated with mosaic variegated aneuploidy (MVA), each extended by 5 kb on either side.

### 1.7 Ranking values

Each of the ten entries in the ranked candidate list carries a ranking value in the challenge's confidence column (`epcr`). The ranking value is ordinal, lies between 0 and 1, and is not a calibrated probability. The value of the top-ranked pair comes from joint rarity and a prior of 1 in 190 (the number of genes with a rare high- or moderate-impact variant); this formula saturates and rounds to 1.0, and it does not reflect variant severity or phase. The same pair in reciprocal order receives 0.95. Every other row receives its gene's phenotype similarity multiplied by call quality, with a factor of 0.5 for SVs outside coding exons, and is capped at 0.80. Annex A1 lists all rows and Section 3 all parameters. Rows 2 to 4 are alternative representations of the top-ranked candidate (the same pair in reciprocal order and each variant alone).

## 2 Results

### 2.1 Candidate genes and variants

The analysis retained seven genes: the five recessive candidates from the genome-wide search and two panel genes with SV calls (Table 1).

**Table 1.** Genes retained by the analysis, with variant class, genotype and classification.

| Gene | Variant class (protein level) | Genotype | Classification |
|---|---|---|---|
| *BUB1B* | Stop-gain p.(Leu737\*); missense p.(Asn1002Lys) | Two heterozygous variants, phase undetermined | Stop-gain: pathogenic/likely pathogenic (ClinVar [12]). Missense: no ClinVar record; a different nucleotide change producing the same amino-acid substitution is listed in ClinVar as of uncertain significance |
| *MAD1L1* | Four SVs outside coding exons: 61 bp insertion, 60 bp deletion, 256 bp deletion, 136 bp insertion | 61 bp insertion homozygous; 60 bp deletion, 256 bp deletion and 136 bp insertion heterozygous | Not assessed; no predicted coding change |
| *CENPE* | 110 bp deletion with 8 bp inserted, outside coding exons | Heterozygous | Not assessed; no predicted coding change |
| *IGBP1* | Missense p.(Met151Leu) | Hemizygous (X chromosome) | Not assessed |
| *DNAH10* | Two missense | Two heterozygous variants, phase undetermined | Not assessed |
| *IZUMO4* | Missense; in-frame insertion | Two heterozygous variants, phase undetermined | Not assessed (no HPO annotation) |
| *ZNF418* | Missense | Homozygous | Not assessed (no HPO annotation) |

*IZUMO4* and *ZNF418* have no HPO annotation, so their phenotype similarity is zero (Table 2) and their ranking value receives the floor of 0.01. *IGBP1* moves from fifth to second place because the final rank is set by phenotype similarity (Section 1.5; Table 2).

**Table 2.** The five recessive candidate genes: variant-based rank and final rank after phenotype ranking.

| Final rank | Variant-based rank | Gene | High-impact variants | Joint rarity | Phenotype similarity |
|---|---|---|---|---|---|
| 1 | 1 | *BUB1B* | 1 | 9.50 | 0.557 |
| 2 | 5 | *IGBP1* | 0 | 4.11 | 0.141 |
| 3 | 3 | *DNAH10* | 0 | 7.42 | 0.024 |
| 4 | 2 | *IZUMO4* | 0 | 8.22 | 0.000 |
| 5 | 4 | *ZNF418* | 0 | 6.32 | 0.000 |

### 2.2 Leading candidate gene: *BUB1B*

MVA (OMIM 257300 [13]) was first linked to *BUB1B* by Hanks et al. [14], who found that affected individuals carried a missense allele together with a truncating or non-expressed allele, inherited from different parents. Matsuura et al. [15] found a monoallelic *BUB1B* mutation in each of seven families with premature chromatid separation (PCS) syndrome, with no second mutation detected but a modest decrease in *BUB1B* transcripts; they regarded PCS syndrome as probably identical to MVA. Sieben et al. [16] summarised the typical pattern in MVA as a nonsense allele combined with a kinase-domain missense allele. The patient's combination of a truncating variant and a kinase-domain missense variant is consistent with this pattern.

#### 2.2.1 Variant characteristics

Both variants are heterozygous (genotype quality 99) and lie on the MANE Select transcript of *BUB1B*, which encodes BUBR1.

The stop-gain variant lies in exon 17 of 23; under the rule summarised by Abou Tayoun et al. [17], the transcript is therefore predicted to undergo nonsense-mediated decay (NMD), although this was not tested. It is rare in gnomAD v4.1 (allele frequency 7.87 × 10⁻⁵, no homozygotes), is classified in ClinVar as pathogenic/likely pathogenic, and lies one codon before the frameshift p.(Ser738fs) reported by Hanks et al. [14] in two MVA families.

The missense variant lies in the kinase domain (UniProt O60566 [18]) and is observed once in gnomAD v4.1. Predictor scores are mixed: CADD 24.5 [19], SIFT 0.01 [20], PolyPhen-2 0.997 [21], AlphaMissense 0.923 [22] and REVEL 0.472 [23]. Under the ClinGen calibration of Pejaver et al. [24], only PolyPhen-2 reaches supporting evidence for pathogenicity, and REVEL, the only one of the four single tools recommended by Pejaver et al. (BayesDel, MutPred2, REVEL, VEST4) that was scored here, gives neither PP3 nor BP4; AlphaMissense was not calibrated. Under ACMG/AMP criteria [25], the variant is of uncertain significance.

#### 2.2.2 Allelic configuration and possible functional impact

Phase is undetermined because the variants are 10,911 bp apart, too far for short reads to span; testing both parents would show whether the two variants lie on different copies of the gene. Whether the missense variant reduces BUBR1 function is not known. Suijkerbuijk et al. [26] showed that BUBR1 catalysis is dispensable for chromosome segregation and that its ATP-binding residues are needed for protein stability; in an earlier study, the same group [27] found that one class of MVA missense alleles, all located in or near the kinase domain, lowers BUBR1 abundance by increasing protein turnover. The patient's missense variant could act in this way. If confirmed experimentally, such an effect would support reclassifying the missense variant from uncertain significance.

### 2.3 Excluded candidates

In the first pass (Section 4.1), *PEX5* and *ADCK5* were classed as rare because the VEP frequency annotation used at that stage failed to match indel alleles. Direct gnomAD v4.1 queries returned an allele frequency above 0.5 for each, and both genes were removed. No SV call retained in the 22-gene panel overlapped a coding exon.

### 2.4 Unresolved findings

- Eight large SV calls were not analysed: five passed Manta's filters but carry symbolic alleles without resolved sequence, and three did not pass the filters (Annex A2).
- No copy-number analysis was performed.
- *TTK* (a panel gene inside one of the large SV calls), *IZUMO4* and *ZNF418* have no HPO annotation; absence of annotation is not evidence against a gene.
- Splicing effects more than 2 bp from an exon–intron boundary were not modelled.
- Exomiser [28], LIRICAL [29] and Genomiser [30] were not run.

## 3 Experimental design

The study is a single-case, retrospective analysis of one patient's sequencing data and eight HPO terms. The analysis was performed blind to the reference variants; because the clinical information already points to the gene involved, the design evaluates prioritisation at the variant level. The pipeline runs the ten stages in Table 3 in a fixed order, from the variant file to the ranked candidate list. The frequency, proximity and gene-family filters were adjusted after a first pass (Section 4.1).

**Table 3.** Stages of the analysis pipeline, with the tool or source, the rule and parameters, and the output of each stage.

| Stage | Tool or source | Rule and parameters | Output |
|---|---|---|---|
| 1 Input | Variant file (Sentieon Haplotyper 202308.02, GATK VariantFiltration); reads re-aligned with BWA-MEM2 2.2.1 | Variants that passed the caller's filters | 4,740,790 variants |
| 2 Target | MANE Select v1.5 (19,293 transcripts) | Coding sequence or within 2 bp of an exon–intron boundary | 22,995 variants; 23,035 alleles |
| 3 Annotation | Custom code; reference coding sequences from Ensembl release 116; VEP as comparison only | High or moderate impact retained | 11,364 alleles |
| 4 Frequency filter | gnomAD v4.1 | Absent, or allele frequency below 1 in 1,000 (higher of exome and genome); a failed query halts the run | — |
| 5 Gene-family exclusion | — | HLA, KIR, OR, MUC, KRTAP, PRB, TAS2R excluded (ZNF not excluded) | 209 variants in 190 genes |
| 6 Recessive model | — | Two heterozygous variants at least 150 bp apart; a homozygous variant; or a hemizygous variant on the non-pseudoautosomal X chromosome (at least 20 calls and at most 25% heterozygous). Heterozygous X calls discarded | 5 candidate genes |
| 7 Variant-based rank | — | Number of high-impact variants, then joint rarity: −log₁₀(4 × f₁ × f₂) for two heterozygous variants (either phase), −log₁₀(f²) for a homozygous variant, −log₁₀(f) for a hemizygous variant; frequency floor 1 in 1,000,000 | Table 2 |
| 8 Phenotype ranking | HPO; eight terms | Coverage of the patient's terms by each gene's annotations, weighted by information content; final rank by phenotype similarity, ties broken by the variant-based rank | Table 2 |
| 9 SV filtering | Manta 1.6.0 | Passed Manta's filters; at least 50 bp; fully resolved sequence; within a 22-gene MVA panel, gene ± 5 kb; at most 6 SV calls ranked | Calls in 2 panel genes |
| 10 Ranked candidate list | — | Top pair: joint rarity with a prior of 1 in 190, saturating at 1.0; reciprocal order 0.95 × top entry. Other entries: phenotype similarity × call quality (genotype quality/99; SV quality/999; × 0.35 in repetitive sequence; × 0.5 for SVs outside coding exons); floor 0.01 for genes without HPO annotation; cap 0.80; tie-break 0.0001 | 10 ranked entries (Annex A1) |

## 4 Limitations and future work

### 4.1 Limitations of the analysis

- This is a single-case analysis, so its performance cannot be generalised.
- The frequency, proximity and gene-family filters were adjusted after a first pass in which the leading candidate ranked fifth of 22 candidate genes; the sensitivity analysis (Annex A3) therefore describes the final pipeline, not its performance.
- Without the frequency filter, the candidate pair falls to eighth place (Annex A3).
- The phenotype-only rank of the leading gene depends strongly on one HPO term, rhabdomyosarcoma (from 14th to 228th of 5,277 genes when it is removed), although the gene stays first among the five candidates (Annex A4).
- No phenotype permutation test was run.
- Only recessive models in MANE Select transcripts were searched; de novo, dominant and low-fraction mosaic variants were not examined, and heterozygous calls on the non-pseudoautosomal X chromosome were discarded.
- Indels were not left-normalised before the gnomAD lookup; this does not affect the two single-nucleotide variants of the leading candidate.
- Ancestry-specific allele frequencies were not used.

### 4.2 Analyses to strengthen the result

Two analyses would strengthen the result most directly: genotyping both parents, or long-read sequencing of the patient, to establish whether the two variants lie on different copies of the gene; and a functional assay of the missense variant. These and further open questions, with the test and data each requires, are listed in Annex A5.

## 5 Reproducibility and data governance

### 5.1 Reproducibility

Installation instructions, and instructions for re-running the pipeline with access to the dataset, are given in the repository README.

### 5.2 Data governance and licences

The patient's data belong to the patient and their family; they were used under the terms of the Hackathon and are not redistributed. All data from the challenge will be deleted by 23 November 2026, within 30 days of the submission deadline of 24 October 2026, as the Official Rules require, and the deletion will be confirmed to the organisers by e-mail. This report, the supporting tables and the source code are released under the Creative Commons Attribution 4.0 International licence (CC BY 4.0; https://creativecommons.org/licenses/by/4.0/). Public resources used in the analysis, such as GRCh38, gnomAD, ClinVar and the HPO, keep their own licences, which are listed in the repository file THIRD_PARTY_LICENSES.md.

### 5.3 Use of generative AI

In compliance with the Mandatory Generative AI Disclosure Guidelines, the following commercially available generative-AI tools were used:

- Provider: Anthropic. Tool: Claude Code. Models: Claude Sonnet 5 and Claude Sonnet 5.5. Plan: Claude Pro (earlier phase), Claude Max (later phase). Setting: model training on my data disabled ("Help improve Claude" off). Use: writing and verifying code, and auditing the repository and the reports.
- Provider: Anthropic. Tool: Claude (claude.ai). Model: Claude Opus 5.5. Plan: Claude Max. Setting: model training on my data disabled ("Help improve Claude" off). Use: editing and drafting the final report.
- Provider: Google. Tool: Gemini. Model: Gemini 3.1 Pro. Plan: Google AI Pro. Setting: activity and model training disabled. Use: generating code.

No dataset file (FASTQ, BAM or VCF) was uploaded or attached to any AI assistant, and no AI assistant analysed the data or produced the ranking. The ranking comes from the deterministic Python pipeline.

## Acknowledgements

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible. We also thank the Wilhelm Foundation.

---

## Annex A. Supplementary tables

### A1 Ranked candidate list

**Table A1.** Ranked candidate list, with ranking values and their derivation.

| Row | Gene | Event | Ranking value | Derivation |
|---|---|---|---|---|
| 1 | *BUB1B* | Pair: stop-gain and missense | 1.0 | Joint rarity with prior 1 in 190; saturates and rounds to 1.0 |
| 2 | *BUB1B* | Same pair, reciprocal order | 0.95 | 1.0 × 0.95 |
| 3 | *BUB1B* | Stop-gain alone | 0.5566 | Phenotype similarity × genotype quality 99/99 |
| 4 | *BUB1B* | Missense alone | 0.5565 | As row 3, minus a 0.0001 tie-break |
| 5 | *CENPE* | 110 bp deletion (8 bp inserted) | 0.1377 | 0.275 × SV quality 999/999 × 0.5 |
| 6 | *MAD1L1* | 61 bp insertion, homozygous | 0.1261 | 0.252 × SV quality 999/999 × 0.5 |
| 7 | *IGBP1* | Missense p.(Met151Leu), hemizygous | 0.0986 | 0.1414 × genotype quality 69/99 |
| 8 | *MAD1L1* | 60 bp deletion | 0.0703 | 0.252 × SV quality 557/999 × 0.5 |
| 9 | *MAD1L1* | 256 bp deletion and 136 bp insertion, both heterozygous, phase undetermined | 0.0342 | 0.252 × lower SV quality/999 × 0.5 |
| 10 | *DNAH10* | Pair of two missense variants | 0.0236 | 0.0236 × genotype quality 99/99 |

Phenotype similarity is shown rounded; products may differ from the listed values in the last digit.

### A2 Large structural variant calls not analysed

Genome-wide, 92 SV calls between 50 bp and 1 Mb that passed Manta's filters overlapped coding sequence; this is a descriptive count, not a pipeline filter. The calls in Table A2 were not analysed for the following reasons. The five calls that passed Manta's filters carry symbolic alleles without resolved sequence, and the pipeline accepts only fully resolved SVs; the three remaining calls did not pass Manta's filters. In addition, megabase-scale calls from short reads are often artefactual, and without a read-depth (copy-number) analysis they could not be confirmed. If any of them were real, for example the 21.7 Mb duplication spanning *CDC20*, the interpretation could change; they are therefore listed as unresolved findings (Section 2.4).

**Table A2.** Large SV calls not analysed (five megabase-scale calls that passed Manta's filters and three flagged calls of unresolved size), with the panel gene each contains and the Manta filter status.

| Type | Size | Gene inside the call | Manta filter status |
|---|---|---|---|
| Duplication | 21.7 Mb | *CDC20* | Passed |
| Duplication | 148.6 Mb | *SHLD3* | Passed |
| Duplication | 41.5 Mb | *MAD2L1BP* | Passed |
| Deletion | 48.8 Mb | *TTK* | Passed (quality 395) |
| Deletion | 43.2 Mb | *CEP192* | Passed (quality 96) |
| Duplication with breakpoints near-identical to the 43.2 Mb deletion | Not available | *CEP192* | Flagged (no paired-read support) |
| Two duplications on chromosome 3 | Not available | *CEP63* | Flagged (ploidy) |

### A3 Sensitivity of the *BUB1B* ranking to each filter

**Table A3.** Rank of the *BUB1B* pair when each filter is changed or removed.

| Change | Variant-based rank | Final rank |
|---|---|---|
| Final configuration | 1 of 5 | 1 |
| Frequency threshold 1 in 100 / 1 in 10,000 | 1 of 26 / 1 of 2 | 1 / 1 |
| Proximity limit 500 bp / 50 bp | 1 of 5 / 2 of 6 | 1 / 1 |
| No gene-family exclusion | 1 of 14 | 1 |
| No proximity filter | 3 of 8 | 1 |
| No frequency filter | 22 of 5,407 | 8 |
| No frequency, proximity or gene-family filter | 69 of 6,263 | 8 |

### A4 Dependence on individual phenotype terms

**Table A4.** Phenotype similarity and phenotype-only rank of *BUB1B* when each of the eight HPO terms is removed in turn. Ranks are among 5,277 HPO-annotated genes; —, not available.

| Term set | *BUB1B* similarity | *BUB1B* phenotype-only rank | *IGBP1* similarity | *DNAH10* similarity | *BUB1B* rank among the five candidates |
|---|---|---|---|---|---|
| All eight terms | 0.557 | 14 | 0.141 | 0.024 | 1 |
| Removal of the rhabdomyosarcoma term, which lowers *BUB1B* most | 0.440 | 228 | — | — | 1 |
| Each of the other seven removals | Above 0.440, up to 0.722 | 6–36 | — | — | 1 |
| Removal that raises *BUB1B* most (one of the seven) | 0.722 | — | 0.184 | 0.031 | 1 |

With all eight terms, the three highest-scoring genes by phenotype alone are *PTCH1* (0.696), *CLCNKB* (0.695) and *ENPP1* (0.666).

### A5 Open questions

**Table A5.** Open questions, with the proposed test and the data required.

| Open question | Proposed test | Data required |
|---|---|---|
| Do the two *BUB1B* variants lie on different copies of the gene? | Parental genotyping, or phasing with reads spanning both variants | Parental DNA; long-read sequencing of the patient |
| Does the missense variant affect BUBR1 function? | Functional assay of the missense variant | Laboratory cell or protein data |
| Are the megabase-scale calls real, and do they affect an MVA gene? | Copy-number analysis | Read depth from the existing alignment |
| Were any rare indels wrongly classed as absent from gnomAD? | Left-normalisation before lookup; count of indels among the 209 rare variants | The patient's variant file (aggregate counts only) |
| Does the ranking depend on chance phenotype overlap? | Phenotype permutation test | HPO annotations |
| Do established prioritisation tools rank the same variants first? | Exomiser or LIRICAL on the same input | The patient's variant file and HPO terms |

## References

1. Freed DN, Aldana R, Weber JA, et al.: The Sentieon Genomics Tools – a fast and accurate solution to variant calling from next-generation sequence data. bioRxiv. 2017; 115717. doi: 10.1101/115717
2. Van der Auwera GA, O'Connor BD: Genomics in the Cloud: Using Docker, GATK, and WDL in Terra. O'Reilly Media. 2020. ISBN 978-1491975190.
3. Gargano MA, Matentzoglu N, Coleman B, et al.: The Human Phenotype Ontology in 2024: phenotypes around the world. Nucleic Acids Res. 2024; 52(D1): D1333–D1346. doi: 10.1093/nar/gkad1005
4. Morales J, Pujar S, Loveland JE, et al.: A joint NCBI and EMBL-EBI transcript set for clinical genomics and research. Nature. 2022; 604(7905): 310–315. doi: 10.1038/s41586-022-04558-8
5. Eilbeck K, Lewis SE, Mungall CJ, et al.: The Sequence Ontology: a tool for the unification of genome annotations. Genome Biol. 2005; 6(5): R44. doi: 10.1186/gb-2005-6-5-r44
6. Dyer SC, Austine-Orimoloye O, Azov AG, et al.: Ensembl 2025. Nucleic Acids Res. 2025; 53(D1): D948–D957. doi: 10.1093/nar/gkae1071
7. McLaren W, Gil L, Hunt SE, et al.: The Ensembl Variant Effect Predictor. Genome Biol. 2016; 17(1): 122. doi: 10.1186/s13059-016-0974-4
8. Chen S, Francioli LC, Goodrich JK, et al.: A genomic mutational constraint map using variation in 76,156 human genomes. Nature. 2024; 625(7993): 92–100. doi: 10.1038/s41586-023-06045-0
9. Whiffin N, Minikel E, Walsh R, et al.: Using high-resolution variant frequencies to empower clinical genome interpretation. Genet Med. 2017; 19(10): 1151–1158. doi: 10.1038/gim.2017.26
10. Vasimuddin M, Misra S, Li H, et al.: Efficient architecture-aware acceleration of BWA-MEM for multicore systems. In: 2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS). IEEE, 2019; 314–324. doi: 10.1109/IPDPS.2019.00041
11. Chen X, Schulz-Trieglaff O, Shaw R, et al.: Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. Bioinformatics. 2016; 32(8): 1220–1222. doi: 10.1093/bioinformatics/btv710
12. Landrum MJ, Chitipiralla S, Kaur K, et al.: ClinVar: updates to support classifications of both germline and somatic variants. Nucleic Acids Res. 2025; 53(D1): D1313–D1321. doi: 10.1093/nar/gkae1090
13. Amberger JS, Bocchini CA, Scott AF, et al.: OMIM.org: leveraging knowledge across phenotype–gene relationships. Nucleic Acids Res. 2019; 47(D1): D1038–D1043. doi: 10.1093/nar/gky1151
14. Hanks S, Coleman K, Reid S, et al.: Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in BUB1B. Nat Genet. 2004; 36(11): 1159–1161. doi: 10.1038/ng1449
15. Matsuura S, Matsumoto Y, Morishima K, et al.: Monoallelic BUB1B mutations and defective mitotic-spindle checkpoint in seven families with premature chromatid separation (PCS) syndrome. Am J Med Genet A. 2006; 140A(4): 358–367. doi: 10.1002/ajmg.a.31069
16. Sieben CJ, Jeganathan KB, Nelson GG, et al.: BubR1 allelic effects drive phenotypic heterogeneity in mosaic-variegated aneuploidy progeria syndrome. J Clin Invest. 2020; 130(1): 171–188. doi: 10.1172/JCI126863
17. Abou Tayoun AN, Pesaran T, DiStefano MT, et al.: Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion. Hum Mutat. 2018; 39(11): 1517–1524. doi: 10.1002/humu.23626
18. UniProt Consortium: UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Res. 2025; 53(D1): D609–D617. doi: 10.1093/nar/gkae1010
19. Schubach M, Maass T, Nazaretyan L, et al.: CADD v1.7: using protein language models, regulatory CNNs and other nucleotide-level scores to improve genome-wide variant predictions. Nucleic Acids Res. 2024; 52(D1): D1143–D1154. doi: 10.1093/nar/gkad989
20. Ng PC, Henikoff S: SIFT: predicting amino acid changes that affect protein function. Nucleic Acids Res. 2003; 31(13): 3812–3814. doi: 10.1093/nar/gkg509
21. Adzhubei IA, Schmidt S, Peshkin L, et al.: A method and server for predicting damaging missense mutations. Nat Methods. 2010; 7(4): 248–249. doi: 10.1038/nmeth0410-248
22. Cheng J, Novati G, Pan J, et al.: Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science. 2023; 381(6664): eadg7492. doi: 10.1126/science.adg7492
23. Ioannidis NM, Rothstein JH, Pejaver V, et al.: REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. Am J Hum Genet. 2016; 99(4): 877–885. doi: 10.1016/j.ajhg.2016.08.016
24. Pejaver V, Byrne AB, Feng BJ, et al.: Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. Am J Hum Genet. 2022; 109(12): 2163–2177. doi: 10.1016/j.ajhg.2022.10.013
25. Richards S, Aziz N, Bale S, et al.: Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genet Med. 2015; 17(5): 405–424. doi: 10.1038/gim.2015.30
26. Suijkerbuijk SJE, van Dam TJP, Karagöz GE, et al.: The vertebrate mitotic checkpoint protein BUBR1 is an unusual pseudokinase. Dev Cell. 2012; 22(6): 1321–1329. doi: 10.1016/j.devcel.2012.03.009
27. Suijkerbuijk SJE, van Osch MHJ, Bos FL, et al.: Molecular causes for BUBR1 dysfunction in the human cancer predisposition syndrome mosaic variegated aneuploidy. Cancer Res. 2010; 70(12): 4891–4900. doi: 10.1158/0008-5472.CAN-09-4319
28. Smedley D, Jacobsen JOB, Jäger M, et al.: Next-generation diagnostics and disease-gene discovery with the Exomiser. Nat Protoc. 2015; 10(12): 2004–2015. doi: 10.1038/nprot.2015.124
29. Robinson PN, Ravanmehr V, Jacobsen JOB, et al.: Interpretable clinical genomics with a likelihood ratio paradigm. Am J Hum Genet. 2020; 107(3): 403–417. doi: 10.1016/j.ajhg.2020.06.021
30. Smedley D, Schubach M, Jacobsen JOB, et al.: A whole-genome analysis framework for effective identification of pathogenic regulatory variants in Mendelian disease. Am J Hum Genet. 2016; 99(3): 595–606. doi: 10.1016/j.ajhg.2016.07.005
