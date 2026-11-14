# NGS Analysis and Quality Control Pipeline

An educational, reproducible next-generation sequencing (NGS) project demonstrating paired-end read alignment, sequencing quality assessment, haploid variant calling, variant normalization, coverage analysis, automated quality reporting, and regression testing.

The project implements two independent workflows using publicly documented training datasets: a mouse ChIP-seq alignment workflow and a simulated bacterial variant-analysis workflow.

## Technologies

- **Languages:** Python, Bash
- **Python libraries:** pandas, pysam, Plotly
- **Read alignment:** BWA-MEM
- **Alignment processing:** SAMtools
- **Variant analysis:** BCFtools
- **Sequencing quality control:** FastQC
- **Environment:** Ubuntu / Windows Subsystem for Linux (WSL)

## Project Workflows

### 1. Mouse ChIP-seq Read Alignment and Quality Control

**Dataset:** Galaxy Training Network paired-end ChIP-seq example.

**Data source:** Wang et al. training dataset, available through Zenodo.

- [Galaxy Training Network tutorial](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/formation_of_super-structures_on_xi/tutorial.html)
- [Zenodo dataset (DOI: 10.5281/zenodo.1324070)](https://doi.org/10.5281/zenodo.1324070)

**Input files:**
- `wt_H3K4me3_read1.fastq.gz`
- `wt_H3K4me3_read2.fastq.gz`

**Workflow:**

FASTQ → FastQC → BWA-MEM → Sorted BAM → BAM Index → Alignment QC

**Quality assessment and troubleshooting:**

FastQC was used to assess paired-end sequencing quality. The investigation identified reads containing ambiguous nucleotide calls, including 318 all-N reads in R2 (0.636% of the 50,000-read dataset).

Alignment-level investigation showed that all 318 all-N R2 reads were unmapped. Of their corresponding R1 mates, 314 mapped and four were unmapped.

**Alignment results:**

- Total alignment records: 100,003
- Primary alignment records: 100,000
- Primary mapped reads: 98,874
- Primary mapping rate: 98.87%
- Properly paired reads: 97.27%
- Mapped reads with MAPQ ≥30: 94.32%

The alignment statistics include separate primary and supplementary alignment records; mapping rates are reported using primary reads where specified.

### 2. Simulated Bacterial Variant Analysis

**Dataset:** Simulated *Staphylococcus aureus* training dataset.

**Data source:** Gladman, Seemann, and Bulach (2017), *Bacterial training dataset for Galaxy training network tutorials on Genome assembly*.

- [Zenodo dataset (DOI: 10.5281/zenodo.582600)](https://doi.org/10.5281/zenodo.582600)

The dataset represents an imaginary bacterium with a miniature reference genome. It is an educational training dataset, not a real clinical isolate.

**Input files:**
- `wildtype.fna`
- `mutant_R1.fastq`
- `mutant_R2.fastq`

**Reference genome length:** 197,394 bp

**Workflow:**

Paired-end FASTQ → BWA-MEM → SAMtools Sorting and Indexing → BCFtools mpileup → Haploid Variant Calling → Preliminary Filtering → Variant Normalization → Variant QC

The pipeline uses BCFtools to call haploid variants, apply preliminary quality and depth filters, and normalize variants against the reference genome.

**Alignment results:**

- Total alignment records: 25,022
- Primary alignment records: 24,960
- Primary mapped reads: 24,960
- Properly paired reads: 99.37%
- Singletons: 0

**Variant-calling results:**

- Total candidate variants: 19
- SNPs: 12
- Indels: 7
- Minimum variant depth: 10
- Minimum QUAL: 185.42
- Minimum alternate allele fraction: 90%

All 19 candidate variants were retained after the pipeline's preliminary filtering.

Reference-allele validation using `bcftools norm -c e` completed successfully.

### 3. Coverage Analysis

Coverage analysis was performed across the 197,394-bp simulated bacterial reference genome.

**Coverage results:**

- Mean sequencing depth: 18.95×
- Reference breadth at ≥1×: 91.52%
- Reference breadth at ≥10×: 90.79%
- Reference breadth at ≥20×: 49.84%
- Bases below 1× coverage: 16,744
- Bases below 10× coverage: 18,181
- Contiguous regions below 10× coverage: 67

These results demonstrate why mean sequencing depth should be evaluated alongside coverage breadth and low-coverage regions.

### 4. Automated Quality Control and Reporting

A Python-based reporting workflow generates an interactive HTML quality-control report using Plotly.

The report includes six visualizations and automated checks covering:

- Mapping and proper-pairing rates
- Variant depth distributions
- Variant QUAL distributions
- Alternate allele fractions
- SNP and indel counts
- Reference coverage and breadth
- Low-coverage regions

The bacterial dataset passed all eight implemented QC checks.

The thresholds are illustrative project-level criteria and are not clinically validated.

## Reproducibility and Regression Testing

An independent execution of the bacterial variant pipeline reproduced all 19 normalized chromosome-position-reference-alternate variant records in the same order.

This confirms reproducibility of the variant identities and ordering for the tested execution; it does not establish independent biological validation.

Three automated regression tests passed:

1. Original dataset accepted.
2. Low-depth test dataset flagged.
3. Low-QUAL test dataset flagged.

The original HTML report was restored after testing.

## Quick start (Python report and regression tests)

Run from the repository root in Ubuntu/WSL or another Python 3.12 environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python scripts/test_qc_pipeline.py
python scripts/create_enhanced_qc_report.py
```

The bundled curated example results let the reporting tests run without raw reads or a reference. The original full analysis requires BWA, SAMtools, BCFtools and FastQC installed separately, and source datasets placed at the paths below. The test runner uses temporary directories and does not overwrite the published HTML report.

**Publication permissions:** The original collaborator's code has no supplied license. Read `LICENSE-NOTICE.md` and obtain explicit permission before making the repository public. See `docs/VALIDATION_AND_RELEASE.md` for outstanding checks.

## Running the Project

### Bacterial Variant Pipeline

From the project directory:

```bash
./scripts/run_variant_pipeline.sh
```

The script performs alignment, BAM processing, variant calling, filtering, normalization, variant statistics, coverage analysis, tabular variant export and HTML QC report generation. Run it from any working directory; it resolves its own repository root.

### Generate the Interactive QC Report

```bash
python scripts/create_enhanced_qc_report.py
```

### Run Regression Tests

```bash
python scripts/test_qc_pipeline.py
```

### Mouse Alignment Pipeline

The mouse alignment workflow uses:

```bash
./scripts/run_alignment_pipeline.sh
```

The mouse FASTQ files, mm10 reference genome, and required BWA reference indexes must be available locally before running this workflow. Both shell scripts resolve their paths relative to their own location.

FastQC is run separately on the paired-end FASTQ files.

## Data Provenance

Both training datasets were verified against the MD5 checksums published by their respective Zenodo records.

**Mouse ChIP-seq dataset:**
- Zenodo record: 1324070
- Both paired-end FASTQ files matched their published MD5 checksums.

**Bacterial variant dataset:**
- Zenodo record: 582600
- Both paired-end FASTQ files and the reference FASTA matched their published MD5 checksums.

Original sequencing data and reference files are retained locally and excluded from the public repository. Check the source records for applicable licensing and redistribution terms before sharing original data files.

## Limitations

- Both workflows use educational training datasets.
- QC thresholds are illustrative and have not undergone clinical validation.
- No clinical-grade analytical validation was performed.
- No independent variant truth-set comparison was performed.
- High mapping rates, sequencing depth, and QUAL scores alone do not establish biological validity.
- Nearby deletions showed read-level co-occurrence evidence, but the discrepancy between read-level counts and VCF allele-depth counts remains unresolved.
- Reproducibility of variant records does not establish biological accuracy.

## Intended Use

This project demonstrates practical bioinformatics skills in:

- NGS data processing and quality assessment
- Paired-end read alignment
- BAM processing and alignment statistics
- Haploid variant calling and normalization
- Sequencing coverage analysis
- Python-based QC reporting
- Workflow automation and reproducibility
- Regression testing and troubleshooting

**This project is intended for educational and portfolio demonstration purposes only. It is not intended for clinical diagnosis or patient-care decisions.**
	
