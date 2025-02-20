#!/usr/bin/python3

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
    

def parse_counts(report, vcf_header):
    headers = open(vcf_header)
    lines = headers.readlines()
    outp = open(report.replace('.csv', '_snps.vcf'), 'w')
    for line in lines:
        outp.write(line)
    outp.write('\n')

    inp = open(report)
    line = inp.readline()
    while line:
        line_array = line.strip().split(',')
        if line_array[0] != '' and line_array[0] != '*':
            if line_array[0] == 'AlleleID':
                print(line_array)
                outp.write('\t'.join(['#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT']))
                # For duplicate sample names, rename them
                samples = []
                for i in line_array[3:]:
                    if i in samples:
                        suffix = samples.count(i) + 1
                        i_updated = i + '_' + str(suffix)
                        outp.write('\t' + i_updated)
                    else:
                        outp.write('\t' + i)
                    samples.append(i)
                outp.write('\n')
            else:
                # snp_position = ['alfalfaRep2vsXJDY1_shared_2402219': ['chr8.1', '20667639', '[G/A]'],...]
                # Write out snp position information
                cloneID = line_array[0].split('|')
                snp_chr = snp_position[cloneID[0]][0]
                snp_bp = snp_position[cloneID[0]][1]
                snp_ref = snp_position[cloneID[0]][2][1]
                snp_alt = snp_position[cloneID[0]][2][3]
                outp.write('\t'.join([snp_chr, snp_bp, cloneID[0], snp_ref, snp_alt, '.', '.']) + '\t')
                # Convert string to integer
                alt_line = inp.readline()
                alt_line_array = alt_line.strip().split(',')
                ref_dp = sum(list(map(int, line_array[15:])))
                alt_dp = sum(list(map(int, alt_line_array[15:])))
                dp = ref_dp + alt_dp
                outp.write('DP=' + str(dp) + ';ADS=' + str(ref_dp) + ',' + str(alt_dp) + '\tDP:RA:AD')
                idx = 3
                while idx < len(line_array):
                    sample_ref_ad = line_array[idx]
                    sample_alt_ad = alt_line_array[idx]
                    sample_dp = str(int(sample_ref_ad) + int(sample_alt_ad))
                    sample_ra = line_array[idx]
                    outp.write('\t' + sample_dp + ':' + sample_ra + ':' + sample_ref_ad + ',' + sample_alt_ad)
                    idx += 1
                outp.write('\n') 
        else:
            pass
        line = inp.readline()


def parse_probe_info(probe_info):
    probe = open(probe_info)
    line = probe.readline() # Header
    # MarkerName	TargetSequence	ReferenceGenome	Chrom	Pos	VariantAllelesDef	Required
    line = probe.readline()
    snp_position = {}
    while line:
        line_array = line.strip().split('\t')
        snp_name = line_array[3] + '_' + line_array[4].zfill(9)
        snp_position[snp_name] = line_array[3:6]
        # snp_position = ['chr8.1_020667639': ['chr8.1', '20667639', '[G/A]'],...]
        line = probe.readline()
    probe.close()
    return(snp_position)


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate VCF file of Ref and Alt read count")

    parser.add_argument('probe_info', help='Probe information in txt format')

    parser.add_argument('vcf_header',
                        help='VCF header')
    
    parser.add_argument('input',
                        help='DArTag read count CSV input file - Original file with header lines')

    args=parser.parse_args()

    snp_position = parse_probe_info(args.probe_info)

    parse_counts(args.input, args.vcf_header, snp_position)
