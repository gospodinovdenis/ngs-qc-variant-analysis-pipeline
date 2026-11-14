#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."

REF="reference/variant_demo/wildtype.fna"
R1="data/variant_demo/mutant_R1.fastq"
R2="data/variant_demo/mutant_R2.fastq"

OUT="results/variant_pipeline"
BAM="$OUT/mutant.sorted.bam"

mkdir -p "$OUT" logs

echo "Starting DNA variant pipeline: $(date)"

# Check required programs
for tool in bwa samtools bcftools; do
    command -v "$tool" >/dev/null || {
        echo "ERROR: Missing tool: $tool"
        exit 1
    }
done

# Check input files
for file in "$REF" "$R1" "$R2"; do
    [[ -s "$file" ]] || {
        echo "ERROR: Missing or empty input: $file"
        exit 1
    }
done

echo "[1/7] Checking reference indexes"

[[ -s "${REF}.fai" ]] || samtools faidx "$REF"

if [[ ! -s "${REF}.bwt" ]]; then
    bwa index "$REF"
fi

echo "[2/7] Aligning DNA reads"

bwa mem -t 2 "$REF" "$R1" "$R2" \
    2> "$OUT/bwa.log" \
    | samtools sort -@ 2 -o "$BAM" -

echo "[3/7] Validating and indexing BAM"

samtools quickcheck -v "$BAM"
samtools index "$BAM"

echo "[4/7] Generating alignment QC"

samtools flagstat "$BAM" > "$OUT/flagstat.txt"
samtools stats "$BAM" > "$OUT/alignment_stats.txt"

echo "[5/7] Calling variants"

bcftools mpileup \
    -f "$REF" \
    -q 20 -Q 20 \
    -a FORMAT/DP,FORMAT/AD \
    -Ou "$BAM" \
    | bcftools call \
        --ploidy 1 \
        -mv \
        -Oz \
        -o "$OUT/raw.vcf.gz"

bcftools index -t "$OUT/raw.vcf.gz"

echo "[6/7] Filtering and normalizing variants"

bcftools view \
    -i 'QUAL>=30 && FMT/DP>=10' \
    -Ou "$OUT/raw.vcf.gz" \
    | bcftools norm \
        -f "$REF" \
        -m -any \
        -Oz \
        -o "$OUT/normalized.vcf.gz"

bcftools index -t "$OUT/normalized.vcf.gz"

echo "[7/7] Generating variant statistics"

bcftools stats "$OUT/normalized.vcf.gz" \
    > "$OUT/variant_stats.txt"

echo "Pipeline completed: $(date)"

echo "Variant summary:"
grep '^SN' "$OUT/variant_stats.txt"

# Export tabular variant calls and coverage metrics for reporting.
# FORMAT/AD can be absent for some calls; inspect the resulting TSV before interpretation.
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%DP\t%AD]\n' \
    "$OUT/normalized.vcf.gz" > "$OUT/variants.tsv"
python3 scripts/analyze_sequencing_coverage.py
python3 scripts/create_enhanced_qc_report.py
