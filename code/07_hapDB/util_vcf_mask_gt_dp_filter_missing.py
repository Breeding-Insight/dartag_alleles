#!/usr/bin/env python3
import pandas as pd
import numpy as np
from datetime import datetime
import argparse



def read_vcf(vcf_file):
    # Find header line
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#CHROM"):
                header = line.strip().lstrip("#").split("\t")
                break
    df = pd.read_csv(vcf_file, comment="#", sep="\t", header=None, names=header)
    fixed_cols = header[:9]
    sample_cols = header[9:]
    return df, fixed_cols, sample_cols


def get_dp_matrix(df, sample_cols):
    # Looks at the FORMAT column for the first variant (iloc[0])
    format_keys = df["FORMAT"].iloc[0].split(":")
    if "DP" not in format_keys:
        raise ValueError("No DP field found in FORMAT column.")
    
    # Finds the index position of "DP" in the format keys list ["GT", "AD", "DP"]
    dp_index = format_keys.index("DP")
    
    # Creates an empty pandas DataFrame to hold DP values
    dp_matrix = pd.DataFrame(index=df.index, columns=sample_cols, dtype=float)
    
    # Empty dictionary to hold ploidy information
    ploidy_dict = {}
    for sample in sample_cols:
        dp_matrix[sample] = df[sample].str.split(":").str[dp_index].replace(".", pd.NA).astype(float)
        # detect ploidy from first non-missing GT string
        non_missing = df[sample].dropna()
        if len(non_missing) > 0:
            ploidy_dict[sample] = len(non_missing.iloc[0].split(":")[0].split("/"))
        else:
            ploidy_dict[sample] = None  # no data to detect ploidy
    return dp_matrix, ploidy_dict


def calculate_missingness(dp_matrix, minDP, maxDP):
    # Create a boolean DataFrame the same shape as dp_matrix
    # True: this genotype should be considered missing
    # False: this genotype passes depth thresholds
    # missing_mask:
    #          S1     S2     S3
    # Var1   False  False   True  # S3 is NaN, missing
    # Var2   True   False   False # S1=5 (<10) missing
    # Var3   True   True    True  # all NaN
    missing_mask = (dp_matrix < minDP) | (dp_matrix > maxDP) | (dp_matrix.isna())
    
    # Missingness rate per variant: proportion of samples missing for that variant
    # axis=1: for each row (variant), take the mean of the boolean values across samples
    # True = 1, False = 0
    # Var1 missing_mask row: [False, False, True] → mean = 1/3 ≈ 0.333  
    # Var2: [True, False, False] → mean = 1/3 ≈ 0.333  
    # Var3: [True, True, True] → mean = 3/3 = 1.0
    # variant_missing_rate: pandas Series (one value per variant)
    # values are proportion of missing samples for that variant
    variant_missing_rate = missing_mask.mean(axis=1)
    
    # Missingness rate per sample: proportion of variants missing for that sample
    # axis=0: for each column (sample), take the mean of the boolean values across variants
    sample_missing_rate = missing_mask.mean(axis=0)
    return missing_mask, variant_missing_rate, sample_missing_rate


def filter_variants_samples(df, missing_mask, variant_missing_rate, sample_missing_rate, marker_missing_cutoff, sample_missing_cutoff):
    # Create a boolean Series indicating which variants pass the missingness threshold
    keep_variants_mask = variant_missing_rate <= marker_missing_cutoff
    
    # Keep only variants where keep_variants_mask is True
    # Create a separate copy, so we don't modify the original df accidentally
    df_filtered = df.loc[keep_variants_mask].copy()
    
    # Apply the same variant filter to the missing mask DataFrame so it matches df_filtered
    missing_mask_filtered = missing_mask.loc[keep_variants_mask]

    # Create a boolean Series for samples
    keep_samples_mask = sample_missing_rate <= sample_missing_cutoff
    kept_sample_cols = [s for s, keep in keep_samples_mask.items() if keep]
    dropped_samples = [s for s, keep in keep_samples_mask.items() if not keep]

    return df_filtered, missing_mask_filtered, kept_sample_cols, keep_variants_mask, dropped_samples


def mask_gt_keep_rest(field_str, ploidy):
    if pd.isna(field_str) or ploidy is None:
        return field_str
    parts = field_str.split(":")
    missing_gt = "/".join(["."] * ploidy)
    parts[0] = missing_gt
    return ":".join(parts)


def apply_masking(df_filtered, kept_sample_cols, missing_mask_filtered, ploidy_dict):
    for sample in kept_sample_cols:
        ploidy = ploidy_dict[sample]
        mask = missing_mask_filtered[sample]
        df_filtered.loc[mask, sample] = df_filtered.loc[mask, sample].apply(
            lambda x: mask_gt_keep_rest(x, ploidy)
        )
    return df_filtered


def write_vcf(df_filtered, fixed_cols, kept_sample_cols, vcf_file, output_vcf):
    final_cols = fixed_cols + kept_sample_cols
    df_out = df_filtered[final_cols]
    with open(output_vcf, "w") as out:
        # write meta lines
        with open(vcf_file) as src:
            for line in src:
                if line.startswith("##"):
                    out.write(line)
                elif line.startswith("#CHROM"):
                    break
        out.write("#" + "\t".join(final_cols) + "\n")
        df_out.to_csv(out, sep="\t", header=False, index=False)


def write_log(log_file, sample_miss_out, marker_miss_out, args, dropped_variants, dropped_samples, sample_missing_rate, missing_mask_filtered, df_filtered):
    with open(log_file, "w") as log:
        log.write(f"VCF DP filtering log - {datetime.now()}\n")
        log.write("Dropped variants: {}\n".format(len(dropped_variants)))
        log.write("Dropped samples: {}\n".format(len(dropped_samples)))
        print(f"\n# VCF DP filtering log - {datetime.now()}")
        print("# Dropped variants: {}".format(len(dropped_variants)))
        print("# Dropped samples: {}".format(len(dropped_samples)))
        
        for k, v in vars(args).items():
            log.write(f"{k}: {v}\n")
        log.write("\nDropped variants: {}\n".format(len(dropped_variants)))
        for chrom, pos in dropped_variants:
            log.write(f"  {chrom}:{pos}\n")
        log.write("\nDropped samples: {}\n".format(len(dropped_samples)))
        for s in dropped_samples:
            log.write(f"  {s} (missing={sample_missing_rate[s]:.3f})\n")
        
        log.write("\nPer-sample missingness (after variant filtering):\n")
        sample_outp = open(sample_miss_out, "w")
        for s in missing_mask_filtered.columns:
            rate = missing_mask_filtered[s].mean()
            sample_outp.write(f"{s},{rate * 100:.2f}\n")
            log.write(f"  {s}: {rate * 100:.2f}\n")
        sample_outp.close()
        
        log.write("\nPer-variant missingness (after filtering):\n")
        marker_outp = open(marker_miss_out, "w")
        for idx, rate in zip(df_filtered.index, missing_mask_filtered.mean(axis=1)):
            log.write(f"  {df_filtered.at[idx, 'CHROM']}:{df_filtered.at[idx, 'POS']}: {rate * 100:.2f}\n")
            markerID = f"{df_filtered.at[idx, 'CHROM']}_{str(df_filtered.at[idx, 'POS']).zfill(9)}"
            marker_outp.write(f"{markerID},{rate * 100:.2f}\n")
        marker_outp.close()
        

# ----------------- Helper Functions -----------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Mask GT based on DP thresholds, keep AD and DP, remove markers/samples by missingness."
    )
    parser.add_argument("--vcf", required=True, help="Input VCF file (uncompressed or .gz)")
    parser.add_argument("--out", required=True, help="Output filtered VCF file")
    parser.add_argument("--minDP", type=int, default=10, help="Minimum DP to keep genotype [default=10]")
    parser.add_argument("--maxDP", type=int, default=200, help="Maximum DP to keep genotype [default=200]")
    parser.add_argument("--marker-missing", type=float, default=0.2,
                        help="Max missingness fraction allowed per variant [default=0.2]")
    parser.add_argument("--sample-missing", type=float, default=0.2,
                        help="Max missingness fraction allowed per sample [default=0.2]")
    return parser.parse_args()


# ----------------- Main -----------------
if __name__ == "__main__":
    args = parse_args()

    # Read VCF
    df, fixed_cols, sample_cols = read_vcf(args.vcf)

    # Get DP matrix and ploidy info
    dp_matrix, ploidy_dict = get_dp_matrix(df, sample_cols)

    # Get missingness masks
    missing_mask, variant_missing_rate, sample_missing_rate = calculate_missingness(
        dp_matrix, args.minDP, args.maxDP)

    # Filter variants and samples
    df_filtered, missing_mask_filtered, kept_sample_cols, keep_variants_mask, dropped_samples = filter_variants_samples(
        df, missing_mask, variant_missing_rate, sample_missing_rate,
        args.marker_missing, args.sample_missing
    )

    # Mask GTs but keep AD, DP
    df_filtered = apply_masking(df_filtered, kept_sample_cols, missing_mask_filtered, ploidy_dict)

    # Dropped variant list
    dropped_variants = df.loc[~keep_variants_mask, ["CHROM", "POS"]].values.tolist()

    # Write VCF
    write_vcf(df_filtered, fixed_cols, kept_sample_cols, args.vcf, args.out)

    
    # Log details
    sample_miss_out = args.out.replace(".vcf", "_sample_miss.csv")
    marker_miss_out = args.out.replace(".vcf", "_marker_miss.csv")
    log = args.out.replace(".vcf", ".log")
    write_log(log, sample_miss_out, marker_miss_out, args, dropped_variants, dropped_samples, sample_missing_rate, missing_mask_filtered, df_filtered)

    print(f"Filtered VCF written to {args.out}")
    print(f"Log written to {log}")
