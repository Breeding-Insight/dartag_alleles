#!/usr/bin/env python3
# Merge multi‑allelic and polyploid VCFs (ploidy > 2 supported).
# Accept different sample sets from each file (union of all samples).
# If the same variant (CHROM,POS,REF) occurs in multiple files with different ALT lists:
# Merge all ALT alleles into one list.
# Safely remap GT allele indices to match the merged ALT ordering.
# Remap AD arrays to match the merged ALT ordering (zero-filled for absent alleles).
# Preserve original GT separators (/ or |).
# Fill missing genotypes with correct ploidy ././… and AD placeholders.
# Sort output by (CHROM, POS) numeric order.
# CSV input supported: pass --vcf-list vcf_list.csv or direct file list.
# Handles empty allele fields (''), partial missing data, and phased genotypes.

import argparse
import gzip
import csv
from collections import defaultdict

def open_vcf(path):
    """Open VCF or VCF.GZ"""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def parse_vcf(vcf_path):
    """Parse VCF into meta lines, sample list, and variant dict."""
    meta = []
    samples = []
    records = {}
    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("##"):
                if not line.startswith("##SAMPLE"):
                    # Preserve non-SAMPLE meta
                    meta.append(line.rstrip())
            elif line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
            else:
                cols = line.strip().split("\t")
                chrom, pos, vid, ref, alt, qual, flt, info, fmt = cols[:9]
                alts = alt.split(",") if alt != "." else []
                gdata = cols[9:]
                records[(chrom, pos, ref)] = {
                    "id": vid, "alt": alts, "qual": qual, "filter": flt,
                    "info": info, "format": fmt, "samples": dict(zip(samples, gdata))
                }
    return meta, samples, records


def safe_allele_to_index(a, mapping):
    """Convert allele string to merged ALT index safely."""
    if a == "." or a == "":
        return "."
    try:
        return mapping.get(int(a), ".")
    except ValueError:
        return "."



def merge_vcfs(file_paths, ploidy):
    all_meta = []
    all_samples = []
    variants = defaultdict(list)

    # Collect meta, samples, and variants
    for idx, path in enumerate(file_paths):
        meta, samples, recs = parse_vcf(path)
        if idx == 0:
            all_meta = meta
        for s in samples:
            if s not in all_samples:
                all_samples.append(s)
        for key, rec in recs.items():
            variants[key].append(rec)

    merged_records = {}
    for key, rec_list in variants.items():
        chrom, pos, ref = key
        all_alts = []
        for rec in rec_list:
            for a in rec["alt"]:
                if a not in all_alts:
                    all_alts.append(a)

        fmt = "GT:AD:DP"
        sample_data = {}
        rec_by_sample = {}
        for rec in rec_list:
            for sname, call in rec["samples"].items():
                rec_by_sample[sname] = (call, rec["alt"])

        for sample in all_samples:
            call_alt = rec_by_sample.get(sample)
            if call_alt:
                call, old_alts = call_alt
                gt, ad, dp = call.split(":")
                sep = "/" if "/" in gt else "|"

                if gt.startswith("."):
                    new_gt = sep.join(["."] * ploidy)
                    new_ad = ",".join(["0"] * (1 + len(all_alts)))  # zeros not .
                    new_dp = "0"
                else:
                    if all_alts == old_alts:
                        new_gt = gt
                        new_ad = ad if ad != "." else ",".join(["0"] * (1 + len(all_alts)))
                    else:
                        mapping = {0: 0}
                        for idx_alt, alt_allele in enumerate(old_alts, start=1):
                            mapping[idx_alt] = all_alts.index(alt_allele) + 1
                        gt_alleles = gt.split(sep)
                        new_gt = sep.join(str(safe_allele_to_index(a, mapping)) for a in gt_alleles)

                        if ad != ".":
                            ad_vals = ad.split(",")
                            new_ad_list = [0] * (1 + len(all_alts))
                            for ai, val in enumerate(ad_vals):
                                if ai == 0:
                                    new_ad_list[0] = val
                                else:
                                    if ai - 1 < len(old_alts):
                                        old_alt_seq = old_alts[ai - 1]
                                        new_index = all_alts.index(old_alt_seq) + 1
                                        new_ad_list[new_index] = val
                            new_ad = ",".join(str(v) for v in new_ad_list)
                        else:
                            new_ad = ",".join(["0"] * (1 + len(all_alts)))
                    new_dp = dp if dp != "." else "0"
            else:
                new_gt = "/".join(["."] * ploidy)
                new_ad = ",".join(["0"] * (1 + len(all_alts)))
                new_dp = "0"
            sample_data[sample] = f"{new_gt}:{new_ad}:{new_dp}"

        merged_records[key] = {
            "chrom": chrom, "pos": pos, "id": ".", "ref": ref,
            "alt": all_alts, "qual": ".", "filter": ".", "info": ".",
            "format": fmt, "samples": sample_data
        }

    return all_meta, all_samples, merged_records


def write_vcf(meta, samples, records, output_path):
    """Write merged VCF with skipped ##SAMPLE lines."""
    out_fh = gzip.open(output_path, "wt") if output_path.endswith(".gz") else open(output_path, "w")
    for m in meta:
        if m.startswith("##SAMPLE="):
            continue
        out_fh.write(m + "\n")
    out_fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
    for key in sorted(records.keys(), key=lambda x: (x[0], int(x[1]))):
        rec = records[key]
        alt_field = ",".join(rec["alt"]) if rec["alt"] else "."
        row = [
            rec["chrom"], rec["pos"], rec["id"], rec["ref"], alt_field,
            rec["qual"], rec["filter"], rec["info"], rec["format"]
        ] + [rec["samples"][s] for s in samples]
        out_fh.write("\t".join(row) + "\n")
    out_fh.close()


def read_vcf_list(csv_path):
    """Read CSV with one path per row."""
    file_paths = []
    with open(csv_path, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            file_paths.append(row[0].strip())
    return file_paths



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge multi-allelic/polyploid VCFs with ALT and sample union (fast generic Python version).")
    parser.add_argument("-o", "--output", required=True, help="Output merged VCF file (.vcf or .vcf.gz)")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--vcf-files", nargs="+", help="List of input VCF files (.vcf or .vcf.gz)")
    group.add_argument("--vcf-list", help="CSV file with one VCF file path per row")
    parser.add_argument("--ploidy", type=int, default=2, help="Ploidy level (default: 2)")
    args = parser.parse_args()

    if args.vcf_files:
        vcf_files = args.vcf_files
    else:
        vcf_files = read_vcf_list(args.vcf_list)

    meta, samples, merged_records = merge_vcfs(vcf_files, args.ploidy)
    write_vcf(meta, samples, merged_records, args.output)

    print(f"[INFO] Merged {len(merged_records)} variants across {len(samples)} samples into {args.output}")
