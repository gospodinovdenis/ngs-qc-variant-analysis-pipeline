#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

REFERENCE="reference/mm10.fa"
R1="data/raw/wt_H3K4me3_read1.fastq.gz"
R2="data/raw/wt_H3K4me3_read2.fastq.gz"

BAM="data/processed/wt_H3K4me3.sorted.bam"

mkdir -p data/processed results/qc logs

LOG="logs/alignment_pipeline.log"

exec > >(tee -a "$LOG") 2>&1

echo "======================================"
echo "NGS ALIGNMENT PIPELINE"
echo "Started: $(date)"
echo "======================================"

for tool in bwa samtools tee; do
    command -v "$tool" >/dev/null || { echo "ERROR: Missing tool: $tool" >&2; exit 1; }
done

# Verify inputs
for file in "$REFERENCE" "$R1" "$R2"; do
    if [[ ! -s "$file" ]]; then
        echo "ERROR: Missing or empty input: $file"
        exit 1
    fi
done

# Verify existing BWA reference index
for extension in amb ann bwt pac sa; do
    if [[ ! -s "${REFERENCE}.${extension}" ]]; then
        echo "ERROR: Missing BWA index: ${REFERENCE}.${extension}"
        exit 1
    fi
done

echo "[1/4] Aligning paired-end reads..."

bwa mem -t 2 "$REFERENCE" "$R1" "$R2" \
    | samtools sort -@ 2 -o "$BAM" -

echo "[2/4] Checking BAM integrity..."

samtools quickcheck -v "$BAM"

echo "[3/4] Indexing BAM..."

samtools index "$BAM"

echo "[4/4] Generating alignment QC..."

samtools flagstat "$BAM" \
    > results/qc/alignment_flagstat.txt

samtools stats "$BAM" \
    > results/qc/alignment_stats.txt

echo "======================================"
echo "Pipeline completed: $(date)"
echo "======================================"
