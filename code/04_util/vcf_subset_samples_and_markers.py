#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd

def get_id_list_from_csv(file, id_column=None):
    """
    Load IDs from a CSV file's specified column.
    If id_column is None, use the first column in the CSV.
    """
    df = pd.read_csv(file)  # header assumed
    if id_column is None:
        id_column = df.columns[0]
        print(f"# Using first column '{id_column}' from {file}")
    elif id_column not in df.columns:
        raise ValueError(f"Column '{id_column}' not found in {file}. "
                         f"Available columns: {list(df.columns)}")
    ids = df[id_column].dropna().astype(str).tolist()
    return ids

def get_vcf_data_lines(vcf, outf):
    """Parse VCF into DataFrame keyed by variant ID."""
    gt_dict = {}
    samples = []
    with open(vcf) as inp, open(outf, 'w') as outp:
        for line in inp:
            if line.startswith('##'):
                outp.write(line)
            elif line.startswith('#CHROM'):
                samples = line.strip().split('\t')
                outp.write(line)  # keep header in output
            else:
                arr = line.strip().split('\t')
                markerID = arr[0].replace('"', '') + '_' + arr[1].zfill(9)
                gt_dict[markerID] = arr
    return pd.DataFrame.from_dict(gt_dict, orient='index', columns=samples)

def subset_samples(df, keep_samples):
    if keep_samples is None:
        return df
    meta_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    drop_cols = [c for c in df.columns if c not in meta_cols + keep_samples]
    return df.drop(columns=drop_cols)

def subset_markers(df, keep_markers):
    if keep_markers is None:
        return df
    marker_set = set(keep_markers)
    return df.loc[df.index.intersection(marker_set)]


def update_INFO_read_counts(df):
    """
    Vectorized replacement for update_read_counts() using correct field mapping:
    Ref = index 1 after split
    Alt = index 2 after split
    DP = last field
    """
    vcf_samples = df.columns[9:]  # sample columns after FORMAT

    # Flatten all sample genotype strings into a Series
    flat = pd.Series(df[vcf_samples].to_numpy().ravel())

    # Split into columns
    parts = flat.str.split(':|,', expand=True)

    # Convert to numeric (NaN→0), then int
    parts = parts.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

    n_markers = len(df)
    n_samples = len(vcf_samples)

    # Reshape each needed part back to (n_markers, n_samples)
    dp_vals  = parts.iloc[:, -1].to_numpy().reshape(n_markers, n_samples)
    ref_vals = parts.iloc[:, 1].to_numpy().reshape(n_markers, n_samples)
    alt_vals = parts.iloc[:, 2].to_numpy().reshape(n_markers, n_samples)

    # Sum across samples for each marker
    dp_sum  = dp_vals.sum(axis=1)
    ref_sum = ref_vals.sum(axis=1)
    alt_sum = alt_vals.sum(axis=1)

    # Assign INFO
    df = df.copy()
    df["INFO"] = [f"DP={dp};ADS={ref},{alt}" for dp, ref, alt in zip(dp_sum, ref_sum, alt_sum)]
    return df

def write_vcf(df, outf):
    df.to_csv(outf, sep="\t", mode='a', index=False)


# ----------------
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Subset VCF by samples and/or markers")
    parser.add_argument('--keep_samples', help='CSV file with samples to keep')
    parser.add_argument('--sample_col', help='Column name in sample CSV (defaults to first column)')
    parser.add_argument('--keep_markers', help='CSV file with markers to keep')
    parser.add_argument('--marker_col', help='Column name in marker CSV (defaults to first column)')
    parser.add_argument('vcf', help='VCF file with read counts')

    args = parser.parse_args()

    if not args.keep_samples and not args.keep_markers:
        parser.error("At least one of --keep_samples or --keep_markers must be provided")

    # Load sample and marker lists as needed
    if args.keep_samples:
        sample_list = get_id_list_from_csv(args.keep_samples, args.sample_col)  
    else:
        sample_list = None
    if sample_list is not None:
        print(f'# Number of samples to keep: {len(sample_list)}')

    if args.keep_markers:
        marker_list = get_id_list_from_csv(args.keep_markers, args.marker_col)
    else:
        marker_list = None
    if marker_list is not None:
        print(f'# Number of markers to keep: {len(marker_list)}')

    # Read VCF
    outf = args.vcf.replace('.vcf', '_sub.vcf')
    df = get_vcf_data_lines(args.vcf, outf)

    # Apply subsetting
    df = subset_samples(df, sample_list)
    df = subset_markers(df, marker_list)

    # Update INFO field
    if args.keep_samples:
        df = update_INFO_read_counts(df)
        print('# Updated read counts\n', df)
    else:
        pass

    # Save
    write_vcf(df, outf)
