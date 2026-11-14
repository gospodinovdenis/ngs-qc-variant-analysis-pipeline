import pandas as pd
import pysam
from pathlib import Path

VCF = "results/variant_pipeline/normalized.vcf.gz"
OUTPUT = Path("results/variant_qc/variant_qc_summary.csv")

records = []

with pysam.VariantFile(VCF) as vcf:
    samples = list(vcf.header.samples)
    if len(samples) != 1:
        raise ValueError("Expected exactly one sample")

    sample_name = samples[0]

    for variant in vcf:
        sample = variant.samples[sample_name]

        ref = variant.ref
        alt = variant.alts[0]

        if len(ref) == 1 and len(alt) == 1:
            variant_type = "SNP"
        elif len(ref) != len(alt):
            variant_type = "INDEL"
        else:
            variant_type = "OTHER"

        dp = sample.get("DP")
        ad = sample.get("AD")

        ref_depth = ad[0] if ad and len(ad) > 0 else None
        alt_depth = ad[1] if ad and len(ad) > 1 else None

        if (
            ref_depth is not None
            and alt_depth is not None
            and ref_depth + alt_depth > 0
        ):
            vaf = alt_depth / (ref_depth + alt_depth)
        else:
            vaf = None

        records.append({
            "chrom": variant.chrom,
            "position": variant.pos,
            "ref": ref,
            "alt": alt,
            "type": variant_type,
            "qual": variant.qual,
            "depth": dp,
            "ref_depth": ref_depth,
            "alt_depth": alt_depth,
            "vaf": vaf,
        })

df = pd.DataFrame(records)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUTPUT, index=False)

print("\nVARIANT QC SUMMARY\n")
print(df.to_string(index=False))

print("\nVariant counts:")
print(df["type"].value_counts())

print(f"\nSaved: {OUTPUT}")
