import pysam

bam_path = "data/processed/wt_H3K4me3.sorted.bam"

total_all_n = 0
mapped_all_n = 0
unmapped_all_n = 0

mate_mapped = 0
mate_unmapped = 0

with pysam.AlignmentFile(bam_path, "rb") as bam:

    for read in bam.fetch(until_eof=True):

        if not read.is_read2:
            continue

        if read.is_secondary or read.is_supplementary:
            continue

        sequence = read.query_sequence

        if not sequence:
            continue

        if set(sequence.upper()) == {"N"}:

            total_all_n += 1

            if read.is_unmapped:
                unmapped_all_n += 1
            else:
                mapped_all_n += 1

            if read.mate_is_unmapped:
                mate_unmapped += 1
            else:
                mate_mapped += 1

print("\nALL-N R2 READ INVESTIGATION\n")

print("Total all-N R2 reads:", total_all_n)
print("Mapped all-N R2 reads:", mapped_all_n)
print("Unmapped all-N R2 reads:", unmapped_all_n)

print("\nMate alignment status:")
print("R1 mate mapped:", mate_mapped)
print("R1 mate unmapped:", mate_unmapped)
