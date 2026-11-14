# NGS Pipeline — Cleanup and Validation Report

**Prepared:** 2026-11-14  
**Basis:** Uploaded README, requirements, scripts.zip, and results.zip.

## Changes made

1. Replaced fixed `$HOME/penn_ngs_pipeline` assumptions in both Bash entry points with paths derived from each script location.
2. Extended the bacterial pipeline to export the seven-column `variants.tsv`, compute sequencing coverage, and regenerate the Plotly QC report after variant calling.
3. Refactored the existing three-case QC regression runner to use temporary output directories; tests no longer overwrite the curated HTML report.
4. Added two standard-library unit tests verifying the published variant count, SNP/indel split, reference length, coverage arithmetic and breadth.
5. Added GitHub Actions Python checks, an explicit `.gitignore`, source-data placeholders, and release/validation documentation.
6. Retained selected source-derived outputs (including the interactive HTML report), while excluding binary alignments, intermediate VCFs, original reads, and references. Removed backup Python files and `__pycache__`.
7. Preserved the original README's numerical claims and educational limitations, adding quick-start and publication-permission guidance.

## Checks executed in this build environment

- `python -m compileall -q scripts tests` — passed.
- `bash -n scripts/*.sh` — passed.
- `python -m unittest discover -s tests -v` — 2 tests passed.
- `python scripts/test_qc_pipeline.py` — 3 tests passed (baseline, low depth, low QUAL).

## Not yet validated

- Neither alignment workflow nor BCFtools calling was executed end-to-end against original FASTQ/FASTA inputs in this build environment.
- GitHub-hosted Actions execution has not been observed; CI configuration must be checked after the first push.
- The nearby-deletion read-count versus FORMAT/AD discrepancy remains unresolved.
- No independent variant truth set or clinical validation has been performed.
- Confirm tool versions, full provenance/checksums, and author approval before public release.
- `requirements.txt` retains the originally supplied pins; test dependency compatibility on the intended clean environment.

## Release decision

**Technically prepared for GitHub staging, not yet approved for public release.** Obtain the original author's permission and select a license before publishing. Do not interpret training QC PASS as clinical validity.
