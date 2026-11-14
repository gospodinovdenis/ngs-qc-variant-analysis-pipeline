import pandas as pd
from pathlib import Path

output = Path("results/qc/alignment_qc_summary.csv")

data = [
    {
        "read": "R1",
        "total_reads": 50000,
        "mapped_reads": 49702,
        "unmapped_reads": 298,
        "mapping_rate": 99.404,
        "mapq30_reads": 46886,
        "mapq30_among_mapped": 94.334,
        "all_n_reads": 0,
    },
    {
        "read": "R2",
        "total_reads": 50000,
        "mapped_reads": 49172,
        "unmapped_reads": 828,
        "mapping_rate": 98.344,
        "mapq30_reads": 46374,
        "mapq30_among_mapped": 94.310,
        "all_n_reads": 318,
    },
]

df = pd.DataFrame(data)

output.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(output, index=False)

print(df.to_string(index=False))
print(f"\nSaved: {output}")
