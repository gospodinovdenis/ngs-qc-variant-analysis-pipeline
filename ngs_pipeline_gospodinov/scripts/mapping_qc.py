import subprocess
import pandas as pd

bam = "data/processed/wt_H3K4me3.sorted.bam"

def count_reads(required_flag, excluded_flag):
    result = subprocess.run(
        [
            "samtools",
            "view",
            "-c",
            "-f", str(required_flag),
            "-F", str(excluded_flag),
            bam
        ],
        capture_output=True,
        text=True,
        check=True
    )

    return int(result.stdout.strip())


total_r1 = count_reads(64, 2304)
total_r2 = count_reads(128, 2304)

mapped_r1 = count_reads(64, 2308)
mapped_r2 = count_reads(128, 2308)

results = pd.DataFrame([
    {
        "Read": "R1",
        "Total": total_r1,
        "Mapped": mapped_r1,
        "Mapping_Percentage": (
            100 * mapped_r1 / total_r1
            if total_r1 else 0
        )
    },
    {
        "Read": "R2",
        "Total": total_r2,
        "Mapped": mapped_r2,
        "Mapping_Percentage": (
            100 * mapped_r2 / total_r2
            if total_r2 else 0
        )
    }
])

print("\nPAIRED-END ALIGNMENT QC\n")
print(results.to_string(index=False))

results.to_csv(
    "results/qc/mapping_summary.csv",
    index=False
)
