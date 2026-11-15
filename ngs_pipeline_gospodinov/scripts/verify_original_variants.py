#!/usr/bin/env python3
"""Compare independently called normalized VCF variants to the committed TSV baseline.

Only CHROM/POS/REF/ALT are compared. QUAL, depth and allele depth are reported
by the original pipeline but are not required to be bit-for-bit identical.
"""
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path


def key(chrom, pos, ref, alt):
    return (chrom, int(pos), ref.upper(), alt.upper())


def read_baseline(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        for line_no, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not row:
                continue
            if len(row) < 4:
                raise ValueError(f"Invalid baseline row {line_no}: {row!r}")
            rows.append(key(*row[:4]))
    return rows


def read_vcf(path):
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError("Malformed VCF record")
            for alt in fields[4].split(","):
                rows.append(key(fields[0], fields[1], fields[3], alt))
    return rows


def kind(record):
    return "SNP" if len(record[2]) == len(record[3]) == 1 else "indel_or_complex"


def main():
    if len(sys.argv) != 5:
        sys.exit("Usage: verify_original_variants.py BASELINE.tsv CALLED.vcf.gz COMPARISON.tsv SUMMARY.md")
    baseline_path, vcf_path, comparison_path, summary_path = map(Path, sys.argv[1:])
    expected = read_baseline(baseline_path)
    observed = read_vcf(vcf_path)
    expected_counts = Counter(expected)
    observed_counts = Counter(observed)
    missing = expected_counts - observed_counts
    extra = observed_counts - expected_counts
    passed = (len(expected) == 19 and len(observed) == 19 and not missing and not extra)
    all_keys = sorted(set(expected_counts) | set(observed_counts))
    comparison_path.parent.mkdir(parents=True, exist_ok=True)
    with comparison_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["chrom", "pos", "ref", "alt", "expected_count", "observed_count", "status"])
        for record in all_keys:
            ec, oc = expected_counts[record], observed_counts[record]
            writer.writerow([*record, ec, oc, "MATCH" if ec == oc else "MISMATCH"])
    summary = [
        "# Phase 1B — Original dataset reproduction",
        "",
        f"**Result:** {'PASS' if passed else 'FAIL'}",
        "",
        f"- Baseline records: {len(expected)}",
        f"- Newly called records: {len(observed)}",
        f"- Baseline SNPs: {sum(kind(r) == 'SNP' for r in expected)}",
        f"- Baseline indels/complex: {sum(kind(r) != 'SNP' for r in expected)}",
        f"- Newly called SNPs: {sum(kind(r) == 'SNP' for r in observed)}",
        f"- Newly called indels/complex: {sum(kind(r) != 'SNP' for r in observed)}",
        f"- Missing variant records: {sum(missing.values())}",
        f"- Unexpected variant records: {sum(extra.values())}",
        "",
        "Comparison uses normalized chromosome, position, reference and alternate allele. ",
        "Allele case is ignored. This is computational reproduction against the committed baseline, not biological truth validation.",
        "",
    ]
    summary_path.write_text("\n".join(summary), encoding="utf-8")
    print("\n".join(summary))
    if not passed:
        sys.exit(1)


if __name__ == "__main__":
    main()
