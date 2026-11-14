import pysam
from collections import Counter

bam_path = "results/variant_pipeline/mutant.sorted.bam"

# Deleted reference bases, using 0-based coordinates
targets = [160547, 160552]

counts = Counter()

with pysam.AlignmentFile(bam_path, "rb") as bam:
    for read in bam.fetch("Wildtype", 160540, 160560):

        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        observations = {}

        for query_pos, ref_pos in read.get_aligned_pairs(
            matches_only=False
        ):
            if ref_pos in targets:
                observations[ref_pos] = (
                    "DEL" if query_pos is None else "BASE"
                )

        if not all(pos in observations for pos in targets):
            continue

        state = tuple(
            observations[pos] == "DEL"
            for pos in targets
        )

        counts[state] += 1

print("\nVALIDATED DELETION PHASING\n")
print("Both deletions:", counts[(True, True)])
print("Only first:", counts[(True, False)])
print("Only second:", counts[(False, True)])
print("Neither:", counts[(False, False)])
print("Total informative alignments:", sum(counts.values()))
