import pandas as pd
from pathlib import Path

out = Path("results/variant_pipeline")

columns = [
    "chrom", "position", "ref", "alt",
    "qual", "depth", "allele_depth"
]

df = pd.read_csv(
    out / "variants.tsv",
    sep="\t",
    names=columns
)

df["type"] = df.apply(
    lambda row: "SNP"
    if len(row["ref"]) == len(row["alt"]) == 1
    else "INDEL",
    axis=1
)

ad = df["allele_depth"].str.split(",", expand=True)
df["ref_depth"] = pd.to_numeric(ad[0])
df["alt_depth"] = pd.to_numeric(ad[1])

df["vaf"] = (
    df["alt_depth"] /
    (df["ref_depth"] + df["alt_depth"])
)

summary = pd.DataFrame({
    "Metric": [
        "Total variants",
        "SNPs",
        "Indels",
        "Mean read depth",
        "Mean VAF",
        "Minimum VAF"
    ],
    "Value": [
        len(df),
        (df["type"] == "SNP").sum(),
        (df["type"] == "INDEL").sum(),
        round(df["depth"].mean(), 2),
        round(df["vaf"].mean(), 3),
        round(df["vaf"].min(), 3)
    ]
})

html = f"""
<html>
<head>
<title>NGS Variant QC Report</title>
<style>
body {{
    font-family: Arial, sans-serif;
    margin: 40px;
    color: #243447;
}}
h1, h2 {{ color: #174a75; }}
table {{
    border-collapse: collapse;
    width: 100%;
    margin-bottom: 25px;
}}
th, td {{
    border: 1px solid #ddd;
    padding: 8px;
    text-align: left;
}}
th {{ background: #eaf2f8; }}
tr:nth-child(even) {{ background: #f7f9fb; }}
</style>
</head>
<body>
<h1>NGS Variant Quality Control Report</h1>
<p>Educational haploid DNA variant analysis</p>
<h2>Variant Summary</h2>
{summary.to_html(index=False)}
<h2>Normalized Variant Calls</h2>
{df.drop(columns=["allele_depth"]).to_html(index=False)}
<p><b>Interpretation:</b> These are candidate variants
from an educational dataset, not clinically validated calls.</p>
</body>
</html>
"""

(out / "variant_qc_report.html").write_text(html)

print("HTML report created successfully.")
print(out / "variant_qc_report.html")
print(summary.to_string(index=False))
