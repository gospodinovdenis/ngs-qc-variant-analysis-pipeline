# NGS Analysis, Variant Calling, and Quality Control Pipeline

## Overview

A reproducible bioinformatics portfolio project demonstrating paired-end sequencing alignment, alignment quality control, haploid variant calling, reference-fragment coverage analysis, variant normalization, automated reporting, and regression testing.

The project contains two independent workflows:

1. Mouse H3K4me3 ChIP-seq read alignment using the mm10 reference genome.
2. Bacterial DNA variant analysis using a 197,394-bp reference fragment.

The workflows demonstrate practical command-line bioinformatics, Python-based analysis, quality assessment, and reproducible reporting.

## Technologies

### Bioinformatics tools

- BWA-MEM 0.7.19-r1273
- SAMtools 1.22.1
- BCFtools 1.22
- FastQC 0.12.1

### Programming and environment

- Python 3.14.4
- pandas 3.1.0
- pysam 0.24.1
- Plotly 7.1.0
- Bash
- Ubuntu on Windows Subsystem for Linux (WSL)

## Workflow 1: Mouse ChIP-seq Alignment and QC

Paired-end H3K4me3 ChIP-seq reads were aligned against the mouse mm10 reference genome.

**Workflow**

FASTQ → BWA-MEM → SAMtools sorting → BAM validation → BAM indexing → Alignment QC

FastQC was used to assess sequencing-read quality before alignment.

### Read-quality investigation

The paired-end training dataset contained 50,000 read pairs.

The second read file contained:
- 432 reads with ambiguous bases (N).
- 318 reads consisting entirely of N bases (0.636% of reads).
- 0.743% ambiguous bases across all read positions.

All 318 all-N second reads were unmapped. Of their corresponding first-read mates, 314 mapped and four were unmapped.

### Alignment results

| Metric | Result |
|---|---:|
| Total alignment records | 100,003 |
| Primary alignment records | 100,000 |
| Primary mapped reads | 98,874 |
| Primary mapping rate | 98.874% |
| Properly paired reads | 97,266 |
| Proper-pair rate | 97.266% |
| Singleton records | 602 |
| Mapped reads with MAPQ ≥ 30 | 94.322% |

The alignment record total includes three supplementary records. Mapping and pairing metrics are reported separately to make their denominators clear.

## Workflow 2: Bacterial Variant Analysis

A paired-end bacterial DNA training dataset was aligned against a 197,394-bp reference fragment identified by the FASTA record name `Wildtype`.

This is a reference fragment, not a complete bacterial genome.

**Workflow**

FASTQ → BWA-MEM → SAMtools sorting and indexing → BCFtools mpileup → Haploid variant calling → Quality filtering → Variant normalization → Variant statistics

### Alignment results

| Metric | Result |
|---|---:|
| Total alignment records | 25,022 |
| Primary alignment records | 24,960 |
| Primary mapped reads | 24,960 |
| Primary mapping rate | 100.00% |
| Properly paired reads | 24,802 |
| Proper-pair rate | 99.37% |
| Singletons | 0 |

The total record count includes 62 supplementary records.

### Variant calling

Variants were called using BCFtools with haploid ploidy.

The workflow applied the following preliminary filtering criteria:

- Variant QUAL ≥ 30
- Read depth (DP) ≥ 10

Variants were normalized against the supplied reference fragment, including multiallelic representation splitting.

| Variant metric | Result |
|---|---:|
| Candidate variants | 19 |
| SNPs | 12 |
| Indels | 7 |
| Minimum variant depth | 10 |
| Minimum QUAL | 185.42 |
| Minimum alternate allele fraction | 90% |

Reference allele validation using BCFtools normalization completed successfully.

These are candidate variants from an educational training dataset. They have not been validated against an independent truth set.

## Reference-Fragment Coverage Analysis

Coverage was calculated across the 197,394-bp reference fragment.

| Coverage metric | Result |
|---|---:|
| Mean depth | 18.95× |
| Breadth ≥ 1× | 91.52% |
| Breadth ≥ 10× | 90.79% |
| Breadth ≥ 20× | 49.84% |
| Bases below 1× | 16,744 |
| Bases below 10× | 18,181 |
| Bases below 20× | 99,022 |
| Contiguous regions below 10× | 67 |

The difference between mean depth and breadth at higher thresholds illustrates why mean sequencing depth alone is insufficient to characterize coverage uniformity.

## Automated QC Reporting

The enhanced HTML report combines alignment, variant, and coverage metrics into an interactive quality-control summary.

Report contents include:

- Mapping and proper-pairing rates
- Variant depth and QUAL distributions
- Alternate allele fractions
- SNP and indel counts
- Mean depth and coverage breadth
- Low-coverage region summaries
- Automated PASS/WARN checks

The report contains six interactive Plotly charts and was visually inspected.

Eight illustrative QC checks passed for the training dataset:

- Mapping rate
- Proper-pairing rate
- Minimum variant depth
- Minimum variant QUAL
- Minimum variant allele fraction
- Mean coverage
- Coverage breadth at ≥1×
- Coverage breadth at ≥10×

These thresholds are educational examples, not clinically validated acceptance criteria.

## Reproducibility and Regression Testing

An independent pipeline execution reproduced the same 19 chromosome-position-reference-alternate variant records in the same order.

This comparison verifies consistency of variant identity and ordering; it does not establish full VCF equivalence or biological accuracy.

Three automated regression tests passed:

1. Original training dataset accepted.
2. Low-depth test dataset flagged.
3. Low-QUAL test dataset flagged.

The regression test workflow restores the original QC report after testing.

## Repository Structure

```text
ngs-qc-variant-analysis-pipeline/
├── README.md
├── requirements.txt
├── .gitignore
├── scripts/
│   ├── analyze_fastq.py
│   ├── analyze_sequencing_coverage.py
│   ├── check_deletion_phasing.py
│   ├── create_enhanced_qc_report.py
│   ├── create_variant_report.py
│   ├── generate_alignment_qc.py
│   ├── investigate_all_n.py
│   ├── mapping_qc.py
│   ├── run_alignment_pipeline.sh
│   ├── run_variant_pipeline.sh
│   ├── test_qc_pipeline.py
│   └── variant_qc.py
└── results/
    └── variant_pipeline/
        ├── coverage_summary.tsv
        ├── low_coverage_regions.tsv
        └── enhanced_qc_report.html
```

Large sequencing files, BAM files, reference indexes, intermediate VCF files, and generated temporary files are excluded from version control.

## Installation and Requirements

Python dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The `requirements.txt` file contains the direct Python dependencies.

BWA, SAMtools, BCFtools, and FastQC must be installed separately and made available on the system PATH.

### External tool versions used

```text
BWA       0.7.19-r1273
SAMtools  1.22.1
BCFtools  1.22
FastQC    0.12.1
Python    3.14.4
```

## Input Data

The original sequencing data and large mm10 reference genome are not included in the repository.

Before running either workflow, obtain the corresponding training data and place the files at the paths expected by the scripts.

Expected input paths:

```text
data/
├── raw/
│   ├── wt_H3K4me3_read1.fastq.gz
│   └── wt_H3K4me3_read2.fastq.gz
└── variant_demo/
    ├── mutant_R1.fastq
    └── mutant_R2.fastq

reference/
├── mm10.fa
└── variant_demo/
    └── wildtype.fna
```

The reference files and sequencing data must match the datasets used to generate the reported results.

**Dataset source citations and redistribution terms should be added here before public release.**

## Running the Workflows

### Mouse ChIP-seq alignment

After placing the mm10 reference, its BWA indexes, and the paired-end reads in the expected locations:

```bash
./scripts/run_alignment_pipeline.sh
```

### Bacterial variant pipeline

After placing the reference fragment and paired-end reads in the expected locations:

```bash
./scripts/run_variant_pipeline.sh
```

### Generate the enhanced QC report

```bash
python scripts/create_enhanced_qc_report.py
```

### Run regression tests

```bash
python scripts/test_qc_pipeline.py
```

## Limitations

- Educational training datasets, not clinical sequencing data.
- The bacterial reference represents a 197,394-bp fragment, not a complete genome.
- QC thresholds are illustrative and have not undergone clinical validation.
- No independent variant truth-set comparison was performed.
- Nearby deletion calls show read-level co-occurrence evidence, but the difference between read-level counts and VCF allele-depth counts remains unresolved.
- High mapping rates and variant QUAL scores alone do not establish biological validity.
- A zero duplicate-flag count does not demonstrate the absence of duplicate reads.

## Intended Use

This project demonstrates practical bioinformatics skills in sequencing data processing, alignment, variant analysis, coverage assessment, QC automation, troubleshooting, and reproducibility.

It is intended for educational and portfolio demonstration purposes only.

**Not intended for clinical diagnosis or clinical decision-making.**
