#!/usr/bin/python3

def get_target_snpID(lut):
    inp = open(lut)
    line = inp.readline()
    target_snpID = []
    while line:
        line_array = line.strip().split(',')
        target_snpID.append(line_array[1])
        line = inp.readline()
    inp.close()
    return(target_snpID)


def get_missing_and_maf(marker):
    marker_array = marker.strip().split('\t')
    # chr1.1	194324	chr1.1_000194324	A	T	.	.	DP=99863;ADS=15261,84602	DP:RA:AD	320:0:0,320
    total_samples = len(marker_array) - 9
    missing_samples = 0
    no_ref_samples = 0
    no_alt_samples = 0
    index = 9
    import re
    while index < len(marker_array):
        dp_array = re.split(':|,', marker_array[index])
        if int(dp_array[0]) == 0:
            missing_samples += 1
        else:
            if int(dp_array[2]) != 0:
                no_ref_samples += 1
            if int(dp_array[3]) != 0:
                no_alt_samples += 1
        index += 1

    if missing_samples == total_samples:
        f_miss = 100.00
        maf = "na"
    else:
        f_miss = round(float(missing_samples) / total_samples, 2)
        alt_af = round(float(no_alt_samples) / (total_samples - missing_samples), 2)
        ref_af = round(float(no_ref_samples) / (total_samples - missing_samples), 2)
        if ref_af > alt_af:
            maf = alt_af
        else:
            maf = ref_af
    #print(marker_array[2], f_miss, ref_af, alt_af, "maf", maf)
    return(f_miss, maf)


def filter_vcf(report, target_snpID, f_miss_threshold, maf_threshold):
    inp = open(report)
    line = inp.readline()
    outp_vcf = open(report.replace('.vcf', '_filter.vcf'), 'w')
    while line:
        # chr1.1	194324	chr1.1_000194324	A	T	.	.	DP=99863;ADS=15261,84602	DP:RA:AD	320:0:0,320
        line_array = line.strip().split('\t')
        if line.startswith('#'):
            outp_vcf.write(line)
        else:
            f_miss, maf = get_missing_and_maf(line)
            # Keep all on-target SNPs
            if line_array[2] in target_snpID:
                outp_vcf.write('\t'.join(line_array[:8]) + ':TARGET=1:MAF=' + str(maf) + '\t' + '\t'.join(line_array[8:]) + '\n')
            else:
            # Only retain off-target SNPs that meet missing data and maf
                if f_miss <= float(f_miss_threshold):
                    if maf >= float(maf_threshold):
                        outp_vcf.write('\t'.join(line_array[:8]) + ':TARGET=0:MAF=' + str(maf) + '\t' + '\t'.join(line_array[8:]) + '\n')
                    else:
                        pass
        line = inp.readline()
    inp.close()
    outp_vcf.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('lut',
                        help='SNP ID lookup table')
    
    parser.add_argument('hap2snp_vcf',
                        help='SNPs extracted from MADC')

    parser.add_argument('f_miss_threshold', help='Threshold of f_miss')

    parser.add_argument('maf_threshold', help='Threshold of maf')

    args=parser.parse_args()
    
    target_snpID = get_target_snpID(args.lut)

    filter_vcf(args.hap2snp_vcf, target_snpID, args.f_miss_threshold, args.maf_threshold)
