import subprocess
import sys
import tempfile
import shutil
from pathlib import Path

# --------------------------------------------------
# Configuration
# --------------------------------------------------

PROJECT = Path(__file__).resolve().parents[1]

SOURCE = PROJECT / "results/variant_pipeline/variants.tsv"

SCRIPT = PROJECT / "scripts/create_enhanced_qc_report.py"

REPORT = (
    PROJECT
    / "results/variant_pipeline/enhanced_qc_report.html"
)

METRICS = [
    "Mapping rate",
    "Properly paired reads",
    "Minimum variant depth",
    "Minimum variant QUAL",
    "Minimum variant VAF",
]


# --------------------------------------------------
# Run enhanced QC report
# --------------------------------------------------

def run_report(input_path):
    """Execute the QC report generator."""

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            str(input_path)
        ],
        cwd=TEST_ROOT,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        raise RuntimeError(
            "QC report generation failed:\n"
            + result.stderr
        )

    return result.stdout


# --------------------------------------------------
# Create modified test datasets
# --------------------------------------------------

def create_modified_file(
    directory,
    name,
    column,
    value
):
    """Modify one value without changing original data."""

    rows = [
        line.split("\t")
        for line in SOURCE.read_text().splitlines()
        if line.strip()
    ]

    if not rows:
        raise ValueError("Source variant table is empty")

    rows[0][column] = str(value)

    path = directory / name

    path.write_text(
        "\n".join(
            "\t".join(row)
            for row in rows
        ) + "\n"
    )

    return path


# --------------------------------------------------
# Extract a specific QC metric
# --------------------------------------------------

def get_metric_line(output, metric_name):
    """Find a QC result row ending in PASS or WARN."""

    for line in output.splitlines():

        stripped = line.strip()

        if not stripped.startswith(metric_name):
            continue

        remainder = stripped[len(metric_name):]

        # Avoid matching a longer, unrelated metric.
        if remainder and not remainder[0].isspace():
            continue

        parts = remainder.split()

        if parts and parts[-1] in ("PASS", "WARN"):
            return stripped

    raise AssertionError(
        f"QC metric not found: {metric_name}\n"
        f"Report output:\n{output}"
    )


def get_metric_status(output, metric_name):
    """Return PASS or WARN for a metric."""

    line = get_metric_line(
        output,
        metric_name
    )

    return line.split()[-1]


def get_metric_value(output, metric_name):
    """Extract the observed metric value."""

    line = get_metric_line(
        output,
        metric_name
    )

    remainder = line[len(metric_name):].strip()

    return remainder.split()[0]


# --------------------------------------------------
# Test 1: Original dataset
# --------------------------------------------------

def test_original_dataset():

    result = run_report(SOURCE)

    assert "Overall training QC: PASS" in result, (
        "Original dataset should pass QC"
    )

    for metric in METRICS:

        status = get_metric_status(
            result,
            metric
        )

        assert status == "PASS", (
            f"Expected PASS for {metric}, got {status}"
        )

    print("PASS: Original dataset accepted")


# --------------------------------------------------
# Test 2: Low-depth variant
# --------------------------------------------------

def test_low_depth(directory):

    low_depth = create_modified_file(
        directory,
        "low_depth.tsv",
        5,
        3
    )

    result = run_report(low_depth)

    assert (
        "Overall training QC: REVIEW REQUIRED"
        in result
    ), "Low-depth dataset should require review"

    status = get_metric_status(
        result,
        "Minimum variant depth"
    )

    value = get_metric_value(
        result,
        "Minimum variant depth"
    )

    assert status == "WARN", (
        f"Expected depth WARN, got {status}"
    )

    assert float(value) == 3.0, (
        f"Expected depth 3, got {value}"
    )

    # All other checks should still pass.
    for metric in METRICS:

        if metric == "Minimum variant depth":
            continue

        assert get_metric_status(
            result,
            metric
        ) == "PASS", (
            f"Unexpected warning for {metric}"
        )

    print("PASS: Low-depth dataset flagged")


# --------------------------------------------------
# Test 3: Low-quality variant
# --------------------------------------------------

def test_low_quality(directory):

    low_quality = create_modified_file(
        directory,
        "low_quality.tsv",
        4,
        5
    )

    result = run_report(low_quality)

    assert (
        "Overall training QC: REVIEW REQUIRED"
        in result
    ), "Low-quality dataset should require review"

    status = get_metric_status(
        result,
        "Minimum variant QUAL"
    )

    value = get_metric_value(
        result,
        "Minimum variant QUAL"
    )

    assert status == "WARN", (
        f"Expected QUAL WARN, got {status}"
    )

    assert float(value) == 5.0, (
        f"Expected QUAL 5, got {value}"
    )

    # All other checks should still pass.
    for metric in METRICS:

        if metric == "Minimum variant QUAL":
            continue

        assert get_metric_status(
            result,
            metric
        ) == "PASS", (
            f"Unexpected warning for {metric}"
        )

    print("PASS: Low-quality dataset flagged")


# --------------------------------------------------
# Main test runner
# --------------------------------------------------

def main():

    print()
    print("NGS QC AUTOMATED REGRESSION TESTS")
    print("---------------------------------")

    if not SOURCE.is_file():
        raise FileNotFoundError(
            f"Missing source file: {SOURCE}"
        )

    if not SCRIPT.is_file():
        raise FileNotFoundError(
            f"Missing report script: {SCRIPT}"
        )

    global TEST_ROOT
    with tempfile.TemporaryDirectory() as tmp:

        directory = Path(tmp)
        TEST_ROOT = directory
        target = directory / "results" / "variant_pipeline"
        target.mkdir(parents=True)
        for name in ("flagstat.txt", "coverage_summary.tsv", "coverage_by_position.tsv", "low_coverage_regions.tsv"):
            shutil.copy2(SOURCE.parent / name, target / name)

        test_original_dataset()

        test_low_depth(directory)

        test_low_quality(directory)

    print()
    print("ALL 3 REGRESSION TESTS PASSED")


# Tests run entirely in a temporary directory and do not overwrite published results.
if __name__ == "__main__":
    main()
