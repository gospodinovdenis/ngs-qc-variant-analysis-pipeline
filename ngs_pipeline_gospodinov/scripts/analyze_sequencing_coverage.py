from pathlib import Path
import subprocess
import csv
import sys

# -------------------------------------------------
# Configuration
# -------------------------------------------------

BAM = Path("results/variant_pipeline/mutant.sorted.bam")
OUT = Path("results/variant_pipeline")
OUT.mkdir(parents=True, exist_ok=True)

REFERENCE = "Wildtype"
LOW_COVERAGE_THRESHOLD = 10

POSITION_FILE = OUT / "coverage_by_position.tsv"
SUMMARY_FILE = OUT / "coverage_summary.tsv"
REGIONS_FILE = OUT / "low_coverage_regions.tsv"

# Exclude:
# 0x004 UNMAP
# 0x100 SECONDARY
# 0x200 QCFAIL
# 0x400 DUP
# 0x800 SUPPLEMENTARY
EXCLUDE_FLAGS = "0xF04"

# -------------------------------------------------
# Validate inputs
# -------------------------------------------------

if not BAM.exists():
    sys.exit(f"ERROR: BAM not found: {BAM}")

if not Path(str(BAM) + ".bai").exists():
    sys.exit("ERROR: BAM index not found. Run samtools index first.")

reference_info = subprocess.check_output(
    ["samtools", "idxstats", str(BAM)],
    text=True
)

reference_length = None

for line in reference_info.splitlines():
    fields = line.split("\t")
    if fields[0] == REFERENCE:
        reference_length = int(fields[1])
        break

if not reference_length:
    sys.exit(f"ERROR: Reference {REFERENCE} not found in BAM.")

print("\nSEQUENCING COVERAGE ANALYSIS")
print("----------------------------")
print(f"Reference: {REFERENCE}")
print(f"Reference length: {reference_length:,} bp")
print(f"Excluded SAM flags: {EXCLUDE_FLAGS}")
print(f"Low-coverage threshold: < {LOW_COVERAGE_THRESHOLD}x")

# -------------------------------------------------
# Generate depth at every reference position
# -------------------------------------------------

command = [
    "samtools",
    "depth",
    "-aa",
    "-d", "0",
    "-G", EXCLUDE_FLAGS,
    "-r", REFERENCE,
    str(BAM)
]

print("\nRunning:")
print(" ".join(command))

process = subprocess.Popen(
    command,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True
)

total_depth = 0
covered_1x = 0
covered_10x = 0
covered_20x = 0
positions_processed = 0

depth_values = []

# Low-coverage interval tracking
regions = []
region_start = None
region_depths = []
previous_position = None

with POSITION_FILE.open("w", newline="") as position_handle:

    writer = csv.writer(position_handle, delimiter="\t")
    writer.writerow(["chrom", "position", "depth"])

    for line in process.stdout:
        chrom, position, depth = line.rstrip().split("\t")

        position = int(position)
        depth = int(depth)

        positions_processed += 1
        total_depth += depth
        depth_values.append(depth)

        if depth >= 1:
            covered_1x += 1

        if depth >= 10:
            covered_10x += 1

        if depth >= 20:
            covered_20x += 1

        writer.writerow([chrom, position, depth])

        # Build contiguous regions below threshold
        if depth < LOW_COVERAGE_THRESHOLD:

            if (
                region_start is None
                or previous_position is None
                or position != previous_position + 1
            ):

                if region_start is not None:
                    regions.append({
                        "chrom": REFERENCE,
                        "start": region_start,
                        "end": previous_position,
                        "length": len(region_depths),
                        "mean_depth": (
                            sum(region_depths) / len(region_depths)
                        ),
                        "min_depth": min(region_depths),
                        "max_depth": max(region_depths)
                    })

                region_start = position
                region_depths = []

            region_depths.append(depth)

        else:

            if region_start is not None:
                regions.append({
                    "chrom": REFERENCE,
                    "start": region_start,
                    "end": previous_position,
                    "length": len(region_depths),
                    "mean_depth": (
                        sum(region_depths) / len(region_depths)
                    ),
                    "min_depth": min(region_depths),
                    "max_depth": max(region_depths)
                })

            region_start = None
            region_depths = []

        previous_position = position

    stderr = process.stderr.read()
    return_code = process.wait()

if return_code != 0:
    sys.exit(f"SAMtools depth failed:\n{stderr}")

# Close final interval if reference ends below threshold
if region_start is not None:
    regions.append({
        "chrom": REFERENCE,
        "start": region_start,
        "end": previous_position,
        "length": len(region_depths),
        "mean_depth": sum(region_depths) / len(region_depths),
        "min_depth": min(region_depths),
        "max_depth": max(region_depths)
    })

# -------------------------------------------------
# Validate position count
# -------------------------------------------------

if positions_processed != reference_length:
    sys.exit(
        "ERROR: Number of depth positions does not match "
        f"reference length. Processed {positions_processed:,}; "
        f"expected {reference_length:,}."
    )

# -------------------------------------------------
# Calculate coverage metrics
# -------------------------------------------------

mean_depth = total_depth / reference_length

breadth_1x = 100 * covered_1x / reference_length
breadth_10x = 100 * covered_10x / reference_length
breadth_20x = 100 * covered_20x / reference_length

below_1x = reference_length - covered_1x
below_10x = reference_length - covered_10x
below_20x = reference_length - covered_20x

# -------------------------------------------------
# Write summary
# -------------------------------------------------

summary = [
    ("Reference", REFERENCE),
    ("Reference length (bp)", reference_length),
    ("Mean depth (x)", round(mean_depth, 4)),
    ("Bases covered >=1x", covered_1x),
    ("Breadth >=1x (%)", round(breadth_1x, 4)),
    ("Bases covered >=10x", covered_10x),
    ("Breadth >=10x (%)", round(breadth_10x, 4)),
    ("Bases covered >=20x", covered_20x),
    ("Breadth >=20x (%)", round(breadth_20x, 4)),
    ("Bases below 1x", below_1x),
    ("Bases below 10x", below_10x),
    ("Bases below 20x", below_20x),
    ("Low-coverage threshold (x)", LOW_COVERAGE_THRESHOLD),
    ("Low-coverage regions", len(regions))
]

with SUMMARY_FILE.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["Metric", "Value"])
    writer.writerows(summary)

# -------------------------------------------------
# Write low-coverage regions
# Coordinates are 1-based inclusive.
# -------------------------------------------------

with REGIONS_FILE.open("w", newline="") as handle:

    fields = [
        "chrom",
        "start",
        "end",
        "length",
        "mean_depth",
        "min_depth",
        "max_depth"
    ]

    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()

    for region in regions:
        region["mean_depth"] = round(region["mean_depth"], 4)
        writer.writerow(region)

# -------------------------------------------------
# Display results
# -------------------------------------------------

print("\nCOVERAGE RESULTS")
print("----------------")
print(f"Reference length: {reference_length:,} bp")
print(f"Mean depth: {mean_depth:.4f}x")
print(f"Breadth >=1x: {breadth_1x:.2f}%")
print(f"Breadth >=10x: {breadth_10x:.2f}%")
print(f"Breadth >=20x: {breadth_20x:.2f}%")
print(f"Bases below 10x: {below_10x:,}")
print(f"Low-coverage regions (<10x): {len(regions):,}")

print("\nOUTPUT FILES")
print("------------")
print(POSITION_FILE)
print(SUMMARY_FILE)
print(REGIONS_FILE)

print("\nCoverage analysis completed successfully.")
