#!/usr/bin/env python3

import sys
import math
import csv
from collections import defaultdict
from scipy.stats import chi2

# ---------------------------
# Extract genotype and allele counts for a site
# ---------------------------
def extract_genotypes(var, ploidy, minDP, sample_stats):
    """
    Extract genotype info for one variant, update per-sample stats, and 
    return allele counts, heterozygote count, missing count, and total samples.
    """
    gts = var.genotypes  # list of allele indices, plus phased flag at end
    try:
        dp_vals = var.format("DP")
    except KeyError:
        dp_vals = None

    allele_counts = {}
    het_count = 0
    missing_count = 0
    total_samples = len(gts)

    for i, gt in enumerate(gts):
        alleles = gt[:ploidy]
        sample_stats[i]['N_VARIANTS'] += 1

        # Missing if all alleles < 0 or fails DP threshold
        if all(a < 0 for a in alleles) or (minDP and dp_vals is not None and dp_vals[i][0] < minDP):
            sample_stats[i]['MISS'] += 1
            continue

        sample_stats[i]['N_CALLED'] += 1
        if dp_vals is not None:
            sample_stats[i]['DP_SUM'] += dp_vals[i][0]

        unique_alleles = set(a for a in alleles if a >= 0)

        # Count alleles
        for a in alleles:
            if a >= 0:
                allele_counts[a] = allele_counts.get(a, 0) + 1

        # Heterozygote if >1 allele present
        if len(unique_alleles) > 1:
            het_count += 1
            sample_stats[i]['HET'] += 1
        else:
            if list(unique_alleles) == [0]:
                sample_stats[i]['HOM_REF'] += 1
            elif list(unique_alleles):
                sample_stats[i]['HOM_ALT'] += 1

    return allele_counts, het_count, missing_count, total_samples

# ---------------------------
# Calculate per-site metrics
# ---------------------------
def calculate_site_stats(var, allele_counts, het_count, missing_count, total_samples, ploidy):
    """
    Compute per-site metrics: MAF_strict, AAF_combined, AAF_max, pi, HOBS, etc.
    """

    total_called = total_samples - missing_count
    total_alleles = sum(allele_counts.values())

    if total_called == 0 or total_alleles == 0:
        return None, "all_missing"

    ref_count = allele_counts.get(0, 0)
    alt_counts_list = [allele_counts.get(i, 0) for i in range(1, len(var.ALT) + 1)]
    alt_freqs_list = [ac / total_alleles for ac in alt_counts_list]

    aaf_combined = sum(alt_counts_list) / total_alleles            # All ALT alleles combined
    aaf_max = max(alt_freqs_list) if alt_freqs_list else 0.0        # Most common ALT frequency

    freqs = [count / total_alleles for count in allele_counts.values()]
    maf_strict = min(freqs) if len(freqs) > 1 else 0.0              # Minimum across *all* alleles
    if maf_strict == 0.0:
        return None, "monomorphic"

    hobs = het_count / total_called
    miss_pct = missing_count / total_samples
    n_alleles_obs = len(allele_counts)

    pi_val = 1 - sum(f**2 for f in freqs)

    # Hardy–Weinberg for diploid, bi-allelic
    hwe_p = None
    if ploidy == 2 and len(var.ALT) == 1:
        n_hom_ref = allele_counts.get(0, 0) // 2
        n_hom_alt = allele_counts.get(1, 0) // 2
        n_het = total_called - n_hom_ref - n_hom_alt
        p = ref_count / total_alleles
        q = 1 - p
        exp = [(p**2)*total_called, (2*p*q)*total_called, (q**2)*total_called]
        obs = [n_hom_ref, n_het, n_hom_alt]
        chi2_stat = sum((o-e)**2/e for o, e in zip(obs, exp) if e > 0)
        hwe_p = 1 - chi2.cdf(chi2_stat, df=1)

    # Store stats
    stats = {
        "CALLED": total_called,
        "MISS_PCT": miss_pct,
        "MAF_strict": maf_strict,
        "AAF_combined": aaf_combined,
        "AAF_max": aaf_max,
        "HOBS": hobs,
        "N_ALLELES": n_alleles_obs,
        "PI": pi_val,
        "MULTI_ALLELIC": len(var.ALT) > 1
    }
    for i, (ac, af) in enumerate(zip(alt_counts_list, alt_freqs_list), start=1):
        stats[f"ALT{i}_AC"] = ac
        stats[f"ALT{i}_AF"] = af
    if hwe_p is not None:
        stats["HWE_P"] = hwe_p
    stats["ALLELE_COUNTS"] = allele_counts

    return stats, None

# ---------------------------
# Write per-site row
# ---------------------------
def write_site_stats_csv(csvwriter, var, total_samples, stats):
    row = {
        "CHROM": var.CHROM,
        "POS": var.POS,
        "ID": var.ID if var.ID else ".",
        "REF": var.REF,
        "ALT": ",".join(var.ALT),
        "N_SAMPLES": total_samples,
        "CALLED": stats["CALLED"],
        "MISS_PCT": f"{stats['MISS_PCT']:.3f}",
        "MAF_strict": f"{stats['MAF_strict']:.4f}",
        "AAF_combined": f"{stats['AAF_combined']:.4f}",
        "AAF_max": f"{stats['AAF_max']:.4f}",
        "HOBS": f"{stats['HOBS']:.4f}",
        "N_ALLELES": stats["N_ALLELES"],
        "PI": f"{stats['PI']:.4f}",
        "MULTI_ALLELIC": stats["MULTI_ALLELIC"]
    }
    for k, v in stats.items():
        if k.startswith("ALT") and (k.endswith("_AC") or k.endswith("_AF")):
            row[k] = v if "_AC" in k else f"{v:.4f}"
    if "HWE_P" in stats:
        row["HWE_P"] = f"{stats['HWE_P']:.4g}"
    row["ALLELE_COUNTS"] = stats["ALLELE_COUNTS"]
    csvwriter.writerow(row)

# ---------------------------
# Write per-sample QC
# ---------------------------
def write_sample_stats_csv(csvwriter, sample_names, sample_stats):
    for sname, sstat in zip(sample_names, sample_stats):
        missing_rate = sstat['MISS'] / sstat['N_VARIANTS'] if sstat['N_VARIANTS'] > 0 else 0
        hobs = sstat['HET'] / sstat['N_CALLED'] if sstat['N_CALLED'] > 0 else 0
        mean_dp = sstat['DP_SUM'] / sstat['N_CALLED'] if sstat['N_CALLED'] > 0 else 0
        csvwriter.writerow({
            "Sample_ID": sname,
            "N_VARIANTS": sstat['N_VARIANTS'],
            "N_CALLED": sstat['N_CALLED'],
            "MISSING_RATE": f"{missing_rate:.4f}",
            "HOBS": f"{hobs:.4f}",
            "HOM_REF": sstat['HOM_REF'],
            "HET": sstat['HET'],
            "HOM_ALT": sstat['HOM_ALT'],
            "MEAN_DP": f"{mean_dp:.2f}"
        })

# ---------------------------
# Main driver
# ---------------------------
def process_vcf(vcf_file, ploidy, minDP=None, minAAF_combined=None, minAAF_max=None):
    from cyvcf2 import VCF, Writer

    vcf = VCF(vcf_file)
    sample_names = list(vcf.samples)
    sample_stats = [defaultdict(int) for _ in sample_names]
    # Count variants
    variant_counter = 0
    for _ in vcf:
        variant_counter += 1
    initial_variants = variant_counter

    # Re-open for further processing
    vcf.close()
    vcf = VCF(vcf_file)

    # Output file paths
    aaf_suffix = str(int(minAAF_combined * 100))
    out_vcf = vcf_file.replace(".vcf", "_aaf" + aaf_suffix + ".vcf").replace(".vcf.gz", "_aaf.vcf.gz")
    out_site_csv = vcf_file.replace(".vcf", "_site_stats.csv").replace(".vcf.gz", "_site_stats.csv")
    out_sample_csv = vcf_file.replace(".vcf", "_sample_stats.csv").replace(".vcf.gz", "_sample_stats.csv")

    w = Writer(out_vcf, vcf)

    # Prepare CSV for site stats
    max_alts = max(len(var.ALT) for var in vcf)
    vcf.close(); vcf = VCF(vcf_file)
    site_fields = ["CHROM","POS","ID","REF","ALT","N_SAMPLES","CALLED","MISS_PCT",
                   "MAF_strict","AAF_combined","AAF_max","HOBS","N_ALLELES","PI","MULTI_ALLELIC"]
    for i in range(1, max_alts + 1):
        site_fields.extend([f"ALT{i}_AC", f"ALT{i}_AF"])
    site_fields.append("HWE_P")
    site_fields.append("ALLELE_COUNTS")
    fout_site = open(out_site_csv, "w", newline='')
    site_writer = csv.DictWriter(fout_site, fieldnames=site_fields)
    site_writer.writeheader()

    # Counters for summary
    kept, drop_monomorphic, drop_aaf_combined, drop_aaf_max = 0, 0, 0, 0

    for var in vcf:
        allele_counts, het_count, missing_count, total_samples = extract_genotypes(var, ploidy, minDP, sample_stats)
        stats, drop_reason = calculate_site_stats(var, allele_counts, het_count, missing_count, total_samples, ploidy)
        if stats is None:
            if drop_reason == "monomorphic": drop_monomorphic += 1
            continue
        if minAAF_combined is not None and stats["AAF_combined"] < minAAF_combined:
            drop_aaf_combined += 1
            continue
        if minAAF_max is not None and stats["AAF_max"] < minAAF_max:
            drop_aaf_max += 1
            continue

        write_site_stats_csv(site_writer, var, total_samples, stats)
        w.write_record(var)
        kept += 1

    fout_site.close(); w.close(); vcf.close()

    # Write per-sample CSV
    fout_sample = open(out_sample_csv, "w", newline='')
    sample_fields = ["Sample_ID","N_VARIANTS","N_CALLED","MISSING_RATE","HOBS","HOM_REF","HET","HOM_ALT","MEAN_DP"]
    sample_writer = csv.DictWriter(fout_sample, fieldnames=sample_fields)
    sample_writer.writeheader()
    write_sample_stats_csv(sample_writer, sample_names, sample_stats)
    fout_sample.close()

    # SUMMARY
    sys.stderr.write("\n=== SUMMARY ===\n")
    sys.stderr.write(f"Initial variants: {variant_counter}\n")
    sys.stderr.write(f"Dropped monomorphic: {drop_monomorphic}\n")
    sys.stderr.write(f"Dropped by minAAF_combined {minAAF_combined}: {drop_aaf_combined}\n")
    sys.stderr.write(f"Dropped by minAAF_max: {drop_aaf_max}\n")
    sys.stderr.write(f"Kept after filtering: {kept}\n")
    sys.stderr.write(f"Filtered VCF: {out_vcf}\n")
    sys.stderr.write(f"Per-site CSV: {out_site_csv}\n")
    sys.stderr.write(f"Per-sample CSV: {out_sample_csv}\n")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="VCF QC: per-site/per-sample metrics, multiallelic aware, optional AAF filters.")
    parser.add_argument('vcf', help="Input VCF (gzipped or uncompressed)")
    parser.add_argument('-p','--ploidy', type=int, required=True, help="Alleles per genotype (e.g., 2=diploid, 6=hexaploid)")
    parser.add_argument('--minDP', type=int, default=None, help="Minimum DP to count genotype as call")
    parser.add_argument('--minAAF_combined', type=float, default=None, help="Filter by minimum combined AAF")
    parser.add_argument('--minAAF_max', type=float, default=None, help="Filter by minimum max single ALT freq")
    args = parser.parse_args()

    process_vcf(args.vcf, args.ploidy, minDP=args.minDP,
                minAAF_combined=args.minAAF_combined,
                minAAF_max=args.minAAF_max)
