#!/usr/bin/env python3
"""Generate deterministic paired-end FASTQs for end-to-end CI (not biological truth)."""
from pathlib import Path
import random

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / 'reference/variant_demo/wildtype.fna'
R1 = ROOT / 'data/variant_demo/mutant_R1.fastq'
R2 = ROOT / 'data/variant_demo/mutant_R2.fastq'
SEED = 20261115
LENGTH = 3000
POSITIONS = (601, 1201, 2001)  # 1-based
COMP = str.maketrans('ACGT', 'TGCA')

def main():
    rng = random.Random(SEED)
    reference = ''.join(rng.choices('ACGT', k=LENGTH))
    mutated = list(reference)
    expected = []
    for pos in POSITIONS:
        ref = reference[pos - 1]
        alt = next(base for base in 'ACGT' if base != ref)
        mutated[pos - 1] = alt
        expected.append((pos, ref, alt))
    mutated = ''.join(mutated)
    for path in (REF, R1, R2):
        path.parent.mkdir(parents=True, exist_ok=True)
    with REF.open('w') as handle:
        handle.write('>Wildtype\n')
        for offset in range(0, LENGTH, 80):
            handle.write(reference[offset:offset+80] + '\n')
    read_length, fragment_length, stride = 150, 350, 8
    with R1.open('w') as forward, R2.open('w') as reverse:
        for index, start in enumerate(range(0, LENGTH-fragment_length+1, stride), 1):
            seq1 = mutated[start:start+read_length]
            seq2 = mutated[start+fragment_length-read_length:start+fragment_length]
            seq2 = seq2.translate(COMP)[::-1]
            for handle, number, sequence in ((forward, 1, seq1), (reverse, 2, seq2)):
                handle.write(f'@ci_read_{index}/{number}\n{sequence}\n+\n{"I"*read_length}\n')
    expected_path = ROOT / 'ci_expected_variants.tsv'
    with expected_path.open('w') as handle:
        handle.write('chrom\tpos\tref\talt\n')
        for pos, ref, alt in expected:
            handle.write(f'Wildtype\t{pos}\t{ref}\t{alt}\n')
    print(f'Generated {index} paired reads and {len(expected)} expected SNPs')

if __name__ == '__main__':
    main()
