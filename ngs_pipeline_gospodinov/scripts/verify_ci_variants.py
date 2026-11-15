#!/usr/bin/env python3
"""Validate called variants against deterministic synthetic ground truth."""
from pathlib import Path
import csv
import subprocess

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = ROOT / 'ci_expected_variants.tsv'
VCF = ROOT / 'results/variant_pipeline/normalized.vcf.gz'
BAM = ROOT / 'results/variant_pipeline/mutant.sorted.bam'


def main():
    if not VCF.is_file() or not BAM.is_file():
        raise SystemExit('FAIL: missing VCF or BAM')
    subprocess.run(['samtools', 'quickcheck', '-v', str(BAM)], check=True)
    subprocess.run(['bcftools', 'index', '-n', str(VCF)], check=True)
    with EXPECTED.open(newline='') as handle:
        expected = {(r['chrom'], int(r['pos']), r['ref'], r['alt']) for r in csv.DictReader(handle, delimiter='\t')}
    result = subprocess.run(['bcftools', 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\n', str(VCF)], capture_output=True, text=True, check=True)
    observed = set()
    rows = [line for line in result.stdout.splitlines() if line.strip()]
    for line in rows:
        chrom, pos, ref, alt = line.split('\t')
        observed.add((chrom, int(pos), ref, alt))
    missing = expected - observed
    unexpected = observed - expected
    print(f'Expected: {len(expected)} variants; observed: {len(rows)} records')
    print('Missing:', sorted(missing))
    print('Unexpected:', sorted(unexpected))
    if missing or unexpected or len(rows) != len(expected):
        raise SystemExit('FAIL: synthetic variant calls differ from expected truth')
    print('PASS: exact synthetic variant identity comparison')

if __name__ == '__main__':
    main()
