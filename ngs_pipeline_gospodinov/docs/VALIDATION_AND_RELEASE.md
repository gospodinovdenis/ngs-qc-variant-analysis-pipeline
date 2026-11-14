# Validation and release checklist

## Verified in the packaged artifact
- Python syntax compilation for packaged scripts.
- Unit checks on curated variant and coverage output examples.
- QC report regression cases executed in a temporary directory (requires Python dependencies).

## Still required before scientific claims or public release
1. Obtain collaborator approval to publish their code, figures and derived outputs, and select an appropriate license.
2. Download source datasets and validate against the Zenodo MD5 values; the raw reads and reference are not included.
3. Execute both workflows end-to-end on a fresh Ubuntu/WSL environment with BWA, SAMtools, BCFtools and FastQC installed; capture version numbers and logs.
4. Reconcile read-level evidence versus FORMAT/AD for nearby deletions; check read filtering, CIGAR alignment, strand and variant normalization.
5. Confirm reported QC thresholds are illustrative, not clinically validated; do not claim independent biological validation.
6. Review outputs for sensitive metadata and confirm dataset redistribution conditions.
7. Enable GitHub Actions and confirm CI passes after first push.

The CI workflow tests Python report logic and published output consistency. It does **not** rerun alignment or variant calling because input reads/references and system bioinformatics binaries are not included.
