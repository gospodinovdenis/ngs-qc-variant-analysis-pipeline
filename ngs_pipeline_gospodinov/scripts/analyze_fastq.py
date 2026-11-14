from pathlib import Path
import gzip
import pandas as pd


def analyze_fastq(filepath):
    total_reads = 0
    all_n_reads = 0
    reads_with_n = 0
    total_bases = 0
    n_bases = 0

    with gzip.open(filepath, "rt") as f:
        while True:
            header = f.readline()

            if not header:
                break

            sequence = f.readline().strip()
            separator = f.readline()
            quality = f.readline()

            if not separator or not quality:
                raise ValueError(
                    f"Incomplete FASTQ record in {filepath}"
                )

            total_reads += 1
            total_bases += len(sequence)

            count_n = sequence.upper().count("N")
            n_bases += count_n

            if count_n > 0:
                reads_with_n += 1

            if sequence and set(sequence.upper()) == {"N"}:
                all_n_reads += 1

    return {
        "file": Path(filepath).name,
        "total_reads": total_reads,
        "reads_with_N": reads_with_n,
        "all_N_reads": all_n_reads,
        "percentage_reads_with_N": (
            100 * reads_with_n / total_reads
            if total_reads else 0
        ),
        "percentage_all_N_reads": (
            100 * all_n_reads / total_reads
            if total_reads else 0
        ),
        "percentage_N_bases": (
            100 * n_bases / total_bases
            if total_bases else 0
        ),
    }


files = [
    "data/raw/wt_H3K4me3_read1.fastq.gz",
    "data/raw/wt_H3K4me3_read2.fastq.gz",
]

results = [analyze_fastq(file) for file in files]

df = pd.DataFrame(results)

print("\nFASTQ QUALITY CONTROL SUMMARY\n")
print(df.to_string(index=False))

Path("results/qc").mkdir(parents=True, exist_ok=True)

df.to_csv(
    "results/qc/fastq_n_content_summary.csv",
    index=False,
)

print("\nQC summary saved successfully.")
