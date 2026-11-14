from pathlib import Path
from datetime import datetime, timezone
import re
import sys

import pandas as pd
import plotly.express as px
from plotly.offline import plot


# ============================================================
# 1. CONFIGURATION
# ============================================================

OUT = Path("results/variant_pipeline")

VCF_TABLE = (
    Path(sys.argv[1])
    if len(sys.argv) > 1
    else OUT / "variants.tsv"
)

FLAGSTAT = OUT / "flagstat.txt"
COVERAGE_SUMMARY = OUT / "coverage_summary.tsv"
COVERAGE_POSITIONS = OUT / "coverage_by_position.tsv"
LOW_COVERAGE = OUT / "low_coverage_regions.tsv"

REPORT = OUT / "enhanced_qc_report.html"

# Illustrative training thresholds only.
# These are NOT clinical acceptance criteria.
MIN_MAPPING_RATE = 95.0
MIN_PROPER_PAIR_RATE = 90.0
MIN_VARIANT_DEPTH = 10
MIN_VARIANT_QUAL = 30
MIN_VARIANT_VAF = 0.90

MIN_MEAN_COVERAGE = 10.0
MIN_BREADTH_1X = 90.0
MIN_BREADTH_10X = 90.0


# ============================================================
# 2. VALIDATE INPUT FILES
# ============================================================

required_files = [
    VCF_TABLE,
    FLAGSTAT,
    COVERAGE_SUMMARY,
    COVERAGE_POSITIONS,
    LOW_COVERAGE,
]

for file_path in required_files:
    if not file_path.exists():
        raise FileNotFoundError(
            f"Required input file not found: {file_path}"
        )

print("\nLoading NGS analysis results...")


# ============================================================
# 3. LOAD VARIANT DATA
# ============================================================

columns = [
    "chrom",
    "position",
    "ref",
    "alt",
    "qual",
    "depth",
    "allele_depth",
]

df = pd.read_csv(
    VCF_TABLE,
    sep="\t",
    names=columns,
)

if df.empty:
    raise ValueError("Variant table is empty.")

df["type"] = df.apply(
    lambda row: (
        "SNP"
        if len(str(row["ref"])) == 1
        and len(str(row["alt"])) == 1
        else "INDEL"
    ),
    axis=1,
)

ad = df["allele_depth"].astype(str).str.split(
    ",",
    expand=True,
)

if ad.shape[1] < 2:
    raise ValueError(
        "Allele-depth column must contain reference and alternate depths."
    )

df["ref_depth"] = pd.to_numeric(ad[0], errors="coerce")
df["alt_depth"] = pd.to_numeric(ad[1], errors="coerce")
df["depth"] = pd.to_numeric(df["depth"], errors="coerce")
df["qual"] = pd.to_numeric(df["qual"], errors="coerce")

if df[["depth", "qual", "ref_depth", "alt_depth"]].isna().any().any():
    raise ValueError("Variant table contains missing or invalid numeric values.")

total_ad = df["ref_depth"] + df["alt_depth"]

df["vaf"] = (
    df["alt_depth"] / total_ad.where(total_ad > 0)
)

df["vaf_percent"] = df["vaf"] * 100


# ============================================================
# 4. LOAD ALIGNMENT STATISTICS
# ============================================================

flagstat_text = FLAGSTAT.read_text()


def extract_flagstat(pattern):
    match = re.search(
        pattern,
        flagstat_text,
        flags=re.MULTILINE,
    )

    if not match:
        raise ValueError(
            f"Could not find flagstat metric: {pattern}"
        )

    return int(match.group(1))


total_records = extract_flagstat(
    r"^(\d+)\s+\+\s+\d+\s+in total"
)

mapped_records = extract_flagstat(
    r"^(\d+)\s+\+\s+\d+\s+mapped\s+\("
)

properly_paired = extract_flagstat(
    r"^(\d+)\s+\+\s+\d+\s+properly paired\s+\("
)

paired_records = extract_flagstat(
    r"^(\d+)\s+\+\s+\d+\s+paired in sequencing"
)

duplicate_records = extract_flagstat(
    r"^(\d+)\s+\+\s+\d+\s+duplicates"
)

mapping_rate = (
    100 * mapped_records / total_records
    if total_records else 0
)

proper_pair_rate = (
    100 * properly_paired / paired_records
    if paired_records else 0
)


# ============================================================
# 5. LOAD WHOLE-REFERENCE COVERAGE DATA
# ============================================================

coverage_df = pd.read_csv(
    COVERAGE_SUMMARY,
    sep="\t",
)

if not {"Metric", "Value"}.issubset(coverage_df.columns):
    raise ValueError("Coverage summary must contain Metric and Value columns.")

coverage_metrics = dict(
    zip(
        coverage_df["Metric"].astype(str),
        coverage_df["Value"],
    )
)


def coverage_value(metric):
    if metric not in coverage_metrics:
        raise ValueError(
            f"Coverage metric missing: {metric}"
        )

    return float(coverage_metrics[metric])


reference_name = str(coverage_metrics["Reference"])

reference_length = int(
    coverage_value("Reference length (bp)")
)

mean_coverage = coverage_value("Mean depth (x)")
breadth_1x = coverage_value("Breadth >=1x (%)")
breadth_10x = coverage_value("Breadth >=10x (%)")
breadth_20x = coverage_value("Breadth >=20x (%)")

bases_covered_1x = int(
    coverage_value("Bases covered >=1x")
)

bases_covered_10x = int(
    coverage_value("Bases covered >=10x")
)

bases_covered_20x = int(
    coverage_value("Bases covered >=20x")
)

bases_below_10x = int(
    coverage_value("Bases below 10x")
)

low_coverage_regions = pd.read_csv(
    LOW_COVERAGE,
    sep="\t",
)

position_df = pd.read_csv(
    COVERAGE_POSITIONS,
    sep="\t",
)

expected_position_columns = {"chrom", "position", "depth"}

if not expected_position_columns.issubset(position_df.columns):
    raise ValueError(
        "Coverage-by-position file is missing required columns."
    )

if len(position_df) != reference_length:
    raise ValueError(
        "Coverage position count does not match reference length: "
        f"{len(position_df)} observed versus {reference_length} expected."
    )

if position_df["depth"].isna().any():
    raise ValueError("Coverage depth contains missing values.")

if len(low_coverage_regions) > 0:
    expected_region_columns = {
        "chrom", "start", "end", "length",
        "mean_depth", "min_depth", "max_depth",
    }

    if not expected_region_columns.issubset(low_coverage_regions.columns):
        raise ValueError(
            "Low-coverage region file is missing required columns."
        )

number_low_coverage_regions = len(low_coverage_regions)


# ============================================================
# 6. AUTOMATED QC CHECKS
# ============================================================

checks = []


def add_check(metric, observed, threshold, passed):
    checks.append({
        "Metric": metric,
        "Observed": observed,
        "Training threshold": threshold,
        "Status": "PASS" if passed else "WARN",
    })


# Alignment checks
add_check(
    "Mapping rate (includes supplementary records)",
    f"{mapping_rate:.2f}%",
    f">= {MIN_MAPPING_RATE}%",
    mapping_rate >= MIN_MAPPING_RATE,
)

add_check(
    "Properly paired reads",
    f"{proper_pair_rate:.2f}%",
    f">= {MIN_PROPER_PAIR_RATE}%",
    proper_pair_rate >= MIN_PROPER_PAIR_RATE,
)

# Variant checks
add_check(
    "Minimum variant depth",
    f'{df["depth"].min():.0f}',
    f">= {MIN_VARIANT_DEPTH}",
    df["depth"].min() >= MIN_VARIANT_DEPTH,
)

add_check(
    "Minimum variant QUAL",
    f'{df["qual"].min():.2f}',
    f">= {MIN_VARIANT_QUAL}",
    df["qual"].min() >= MIN_VARIANT_QUAL,
)

add_check(
    "Minimum variant VAF",
    f'{df["vaf_percent"].min():.2f}%',
    f">= {MIN_VARIANT_VAF * 100:.0f}%",
    df["vaf"].min() >= MIN_VARIANT_VAF,
)

# Whole-reference coverage checks
add_check(
    "Mean sequencing coverage",
    f"{mean_coverage:.2f}x",
    f">= {MIN_MEAN_COVERAGE:.1f}x",
    mean_coverage >= MIN_MEAN_COVERAGE,
)

add_check(
    "Breadth of coverage >=1x",
    f"{breadth_1x:.2f}%",
    f">= {MIN_BREADTH_1X:.1f}%",
    breadth_1x >= MIN_BREADTH_1X,
)

add_check(
    "Breadth of coverage >=10x",
    f"{breadth_10x:.2f}%",
    f">= {MIN_BREADTH_10X:.1f}%",
    breadth_10x >= MIN_BREADTH_10X,
)

qc = pd.DataFrame(checks)

overall = (
    "PASS"
    if (qc["Status"] == "PASS").all()
    else "REVIEW REQUIRED"
)


# ============================================================
# 7. CREATE VARIANT CHARTS
# ============================================================

fig_types = px.bar(
    df["type"].value_counts()
    .rename_axis("Variant type")
    .reset_index(name="Count"),
    x="Variant type",
    y="Count",
    title="SNP and Indel Distribution",
    text="Count",
)

fig_depth = px.histogram(
    df,
    x="depth",
    color="type",
    nbins=12,
    title="Variant Read Depth Distribution",
    labels={"depth": "Read depth"},
)

fig_vaf = px.scatter(
    df,
    x="depth",
    y="vaf_percent",
    color="type",
    hover_data=["position", "ref", "alt", "qual"],
    title="Variant Depth vs Alternate Allele Fraction",
    labels={
        "depth": "Read depth",
        "vaf_percent": "Alternate allele fraction (%)",
    },
)

fig_vaf.update_yaxes(range=[0, 105])

fig_qual = px.histogram(
    df,
    x="qual",
    color="type",
    nbins=12,
    title="Variant QUAL Distribution",
    labels={"qual": "Variant QUAL"},
)


# ============================================================
# 8. CREATE COVERAGE CHARTS
# ============================================================

fig_coverage_hist = px.histogram(
    position_df,
    x="depth",
    nbins=40,
    title="Whole-Reference Sequencing Depth Distribution",
    labels={
        "depth": "Read depth (x)",
        "count": "Reference positions",
    },
)

# Calculate mean depth in 500-bp windows.
coverage_plot_df = (
    position_df.assign(
        window=(
            (position_df["position"] - 1) // 500
        ) * 500 + 1
    )
    .groupby("window", as_index=False)
    .agg(
        mean_depth=("depth", "mean"),
        min_depth=("depth", "min"),
        max_depth=("depth", "max"),
    )
)

fig_coverage_profile = px.line(
    coverage_plot_df,
    x="window",
    y="mean_depth",
    title="Sequencing Coverage Across Reference",
    labels={
        "window": "Reference position (bp)",
        "mean_depth": "Mean depth (x)",
    },
    hover_data=["min_depth", "max_depth"],
)

fig_coverage_profile.add_hline(
    y=10,
    line_dash="dash",
    annotation_text="10x threshold",
)

fig_coverage_profile.add_hline(
    y=20,
    line_dash="dot",
    annotation_text="20x threshold",
)

fig_coverage_profile.update_yaxes(rangemode="tozero")


# ============================================================
# 9. RENDER EACH CHART EXACTLY ONCE
# ============================================================

# The first chart loads Plotly JavaScript.
# The remaining charts reuse that same JavaScript library.
all_figures = [
    ("coverage_hist", fig_coverage_hist),
    ("coverage_profile", fig_coverage_profile),
    ("variant_types", fig_types),
    ("variant_depth", fig_depth),
    ("variant_vaf", fig_vaf),
    ("variant_qual", fig_qual),
]

chart_html = {}

for index, (chart_name, figure) in enumerate(all_figures):
    chart_html[chart_name] = plot(
        figure,
        include_plotlyjs=(index == 0),
        output_type="div",
    )


# ============================================================
# 10. FORMAT HTML TABLES
# ============================================================

def style_status(value):
    if value == "PASS":
        return "background-color: #d8f3dc"

    return "background-color: #fff3cd"


qc_html = (
    qc.style
    .map(style_status, subset=["Status"])
    .to_html(index=False)
)

variant_html = df[
    [
        "chrom",
        "position",
        "ref",
        "alt",
        "type",
        "qual",
        "depth",
        "ref_depth",
        "alt_depth",
        "vaf_percent",
    ]
].round(2).to_html(
    index=False,
    classes="data-table",
)

if not low_coverage_regions.empty:
    low_coverage_display = (
        low_coverage_regions
        .sort_values(
            ["length", "mean_depth"],
            ascending=[False, True],
        )
        .head(20)
        .copy()
    )

    low_coverage_html = low_coverage_display.round(
        {"mean_depth": 2}
    ).to_html(
        index=False,
        classes="data-table",
    )
else:
    low_coverage_html = (
        "<p>No positions below the configured coverage threshold.</p>"
    )


# ============================================================
# 11. GENERATE HTML REPORT
# ============================================================

generated_at = datetime.now(
    timezone.utc
).strftime("%Y-%m-%d %H:%M UTC")

report_html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Enhanced NGS Quality Control Report</title>

<style>
body {{
    font-family: Arial, sans-serif;
    max-width: 1150px;
    margin: 40px auto;
    padding: 0 20px;
    color: #263238;
    line-height: 1.5;
}}

h1, h2, h3 {{
    color: #164b73;
}}

.summary {{
    background: #eef4f8;
    padding: 20px;
    border-radius: 8px;
    margin-bottom: 25px;
}}

table {{
    border-collapse: collapse;
    width: 100%;
    margin: 20px 0;
    font-size: 14px;
}}

th, td {{
    padding: 9px;
    border: 1px solid #ddd;
    text-align: left;
}}

th {{
    background: #eaf2f8;
}}

.warning {{
    background: #fff3cd;
    padding: 15px;
    margin-top: 25px;
    border-radius: 6px;
}}

.metric-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
    gap: 12px;
    margin: 20px 0;
}}

.metric {{
    background: #f1f5f8;
    padding: 15px;
    border-radius: 7px;
}}

.metric strong {{
    display: block;
    font-size: 24px;
    color: #164b73;
    margin-bottom: 5px;
}}

.metric span {{
    font-size: 13px;
    color: #52616b;
}}

.small-note {{
    color: #52616b;
    font-size: 13px;
}}

.chart {{
    margin: 25px 0;
}}
</style>
</head>

<body>

<h1>Enhanced NGS Quality Control Report</h1>
<p>Generated: {generated_at}</p>

<div class="summary">
    <h2>Pipeline Summary</h2>
    <p><b>Overall training QC:</b> {overall}</p>
    <p><b>Total candidate variants:</b> {len(df)}</p>
    <p><b>SNPs:</b> {(df["type"] == "SNP").sum()}</p>
    <p><b>Indels:</b> {(df["type"] == "INDEL").sum()}</p>
    <p><b>Mapped alignment records:</b> {mapped_records:,} / {total_records:,}</p>
    <p><b>Mapping rate (includes supplementary records):</b> {mapping_rate:.2f}%</p>
    <p><b>Properly paired rate:</b> {proper_pair_rate:.2f}%</p>
    <p><b>Duplicate-marked records:</b> {duplicate_records:,}</p>
    <p class="small-note">
        Zero duplicate-marked records means no records were flagged as duplicates;
        it does not establish that the reads contain no biological or technical duplicates.
    </p>
</div>

<h2>Automated QC Checks</h2>
{qc_html}

<h2>Whole-Reference Sequencing Coverage</h2>

<p>
    Reference: <b>{reference_name}</b><br>
    Reference length: <b>{reference_length:,} bp</b>
</p>

<div class="metric-grid">
    <div class="metric">
        <strong>{mean_coverage:.2f}x</strong>
        <span>Mean sequencing depth</span>
    </div>

    <div class="metric">
        <strong>{breadth_1x:.2f}%</strong>
        <span>Breadth at 1x ({bases_covered_1x:,} bases)</span>
    </div>

    <div class="metric">
        <strong>{breadth_10x:.2f}%</strong>
        <span>Breadth at 10x ({bases_covered_10x:,} bases)</span>
    </div>

    <div class="metric">
        <strong>{breadth_20x:.2f}%</strong>
        <span>Breadth at 20x ({bases_covered_20x:,} bases)</span>
    </div>
</div>

<p>
    <b>Bases below 10x:</b> {bases_below_10x:,}<br>
    <b>Contiguous low-coverage regions below 10x:</b>
    {number_low_coverage_regions:,}
</p>

<h3>Coverage Visualizations</h3>

<div class="chart">
{chart_html["coverage_hist"]}
</div>

<div class="chart">
{chart_html["coverage_profile"]}
</div>

<h3>Longest Low-Coverage Regions</h3>

<p class="small-note">
    Showing up to 20 longest contiguous intervals below 10x.
    Coordinates are 1-based and inclusive. The complete list is
    available in low_coverage_regions.tsv.
</p>

{low_coverage_html}

<h2>Interactive Variant QC Visualizations</h2>

<div class="chart">
{chart_html["variant_types"]}
</div>

<div class="chart">
{chart_html["variant_depth"]}
</div>

<div class="chart">
{chart_html["variant_vaf"]}
</div>

<div class="chart">
{chart_html["variant_qual"]}
</div>

<h2>Detailed Variant Calls</h2>
{variant_html}

<div class="warning">
    <b>Important interpretation notes:</b>
    <ul>
        <li>This is an educational NGS analysis using training data.</li>
        <li>QC thresholds are illustrative and have not been validated
            for clinical diagnostic use.</li>
        <li>Coverage metrics describe the supplied 197,394-bp reference
            fragment, not a complete bacterial genome.</li>
        <li>The overall training QC status does not establish biological
            truth or clinical reportability.</li>
        <li>Variant quality, biological validity, and clinical significance
            require additional evidence and review.</li>
    </ul>
</div>

</body>
</html>
"""


# ============================================================
# 12. SAVE REPORT AND PRINT SUMMARY
# ============================================================

REPORT.write_text(
    report_html,
    encoding="utf-8",
)

print("\nENHANCED NGS QUALITY CONTROL REPORT")
print("-----------------------------------")
print(f"Variants: {len(df)}")
print(f"SNPs: {(df['type'] == 'SNP').sum()}")
print(f"Indels: {(df['type'] == 'INDEL').sum()}")
print(f"Mapping rate: {mapping_rate:.2f}%")
print(f"Properly paired rate: {proper_pair_rate:.2f}%")
print(f"Mean coverage: {mean_coverage:.2f}x")
print(f"Breadth >=1x: {breadth_1x:.2f}%")
print(f"Breadth >=10x: {breadth_10x:.2f}%")
print(f"Breadth >=20x: {breadth_20x:.2f}%")
print(f"Low-coverage regions: {number_low_coverage_regions}")
print(f"Overall training QC: {overall}")

print("\nQC CHECKS")
print(qc.to_string(index=False))

print(f"\nReport saved: {REPORT}")
