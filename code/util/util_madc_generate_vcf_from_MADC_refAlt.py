#!/usr/bin/python3

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

def parse_probe_info(probe_info):
    probe = open(probe_info)
    line = probe.readline() # Header
    # MarkerName	TargetSequence	ReferenceGenome	Chrom	Pos	VariantAllelesDef	Required
    line = probe.readline()
    snp_position = {}
    while line:
        line_array = line.strip().split('\t')
        alleles = line_array[5].replace('[', '').replace(']', '').split('/')
        # If there are indels in the markers
        if alleles[0] == '-':
            snp_position[line_array[3] + '_' + line_array[4].zfill(9)] = line_array[3:5] + ['.', alleles[1]]
        elif alleles[1] == '-':
            snp_position[line_array[3] + '_' + line_array[4].zfill(9)] = line_array[3:5] + [alleles[0], '.']
        else:
            snp_position[line_array[3] + '_' + line_array[4].zfill(9)] = line_array[3:5] + alleles
        # snp_position = {'alfalfaRep2vsXJDY1_shared_2402219': ['chr8.1', '20667639', 'G', 'A'],...}
        line = probe.readline()
    probe.close()
    return(snp_position)


def get_bottom_loci(botloci):
    inp = open(botloci)
    lines = inp.readlines()
    bottom_loci = [line.strip() for line in lines]
    return(bottom_loci)


def rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-':'-'}
    rev_seq = ''
    for base in seq:
        base_com = complement.get(base)
        rev_seq = base_com + rev_seq
    return(rev_seq)


def collect_ref_alt_read_counts(report):
    inp = open(report)
    line = inp.readline()
    ref_alt_read_counts = {}
    while line:
        line_array = line.strip().split(',')
        # chr1.1_000194324|Alt_0002,chr1.1_000194324,CGAAATAATAACCCAAGTTCTGCCAGTTTATGTTAAAACTTTTCTTACAAGGTACAAGTTCGGTGACAACTTAACAAGTAA,16
        if line_array[0] != '' and line_array[0] != '*':
            if line_array[0] == 'AlleleID':
                print('Header line:', line_array[0:6])
                sample_line = ['#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT'] + line_array[3:]
            else:
                if line_array[0].endswith("Ref_0001") or line_array[0].endswith("Alt_0002"):
                    if line_array[1] not in ref_alt_read_counts:
                        ref_alt_read_counts[line_array[1]] = [line_array]
                    else:
                        ref_alt_read_counts[line_array[1]].append(line_array)
        else:
            pass
        line = inp.readline()
    inp.close()
    return(ref_alt_read_counts, sample_line)


def generate_vcf(report, sample_line, ref_alt_read_counts, vcf_header, snp_position, bottom_loci):
    headers = open(vcf_header)
    lines = headers.readlines()
    outp = open(report.replace('.csv', '_snps_target.vcf'), 'w')
    for line in lines:
        outp.write(line)
    outp.write('\t'.join(sample_line) + '\n')   
    # snp_position = ['alfalfaRep2vsXJDY1_shared_2402219': ['chr8.1', '20667639', '[G/A]'],...]
    # Write out snp position information
    for key, value in ref_alt_read_counts.items():
        if value[0][0].endswith("Ref_0001"):
            ref_count = value[0]
            alt_count = value[1]
        elif value[0][0].endswith("Alt_0002"):
            ref_count = value[1]
            alt_count = value[0]
        else:  
            print('Check this marker locus:', value)
            break
    
        cloneID = ref_count[1]
        snp_chr = snp_position[cloneID][0]
        snp_bp = snp_position[cloneID][1]
        if cloneID in bottom_loci:
            snp_ref = rev_complement(snp_position[cloneID][2])
            snp_alt = rev_complement(snp_position[cloneID][3])
        else:
            snp_ref = snp_position[cloneID][2]
            snp_alt = snp_position[cloneID][3]
        outp.write('\t'.join([snp_chr, snp_bp, cloneID, snp_ref, snp_alt, '.', '.']) + '\t')
        ref_dp = sum(list(map(int, ref_count[3:])))
        alt_dp = sum(list(map(int, alt_count[3:])))
        dp = ref_dp + alt_dp
        outp.write('DP=' + str(dp) + ';ADS=' + str(ref_dp) + ',' + str(alt_dp) + '\tDP:RA:AD')
        idx = 3
        while idx < len(ref_count):
            sample_ref_ad = ref_count[idx]
            sample_alt_ad = alt_count[idx]
            sample_dp = str(int(sample_ref_ad) + int(sample_alt_ad))
            sample_ra = ref_count[idx]
            outp.write('\t' + sample_dp + ':' + sample_ra + ':' + sample_ref_ad + ',' + sample_alt_ad)
            idx += 1
        outp.write('\n')



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate VCF file of Ref and Alt read count")

    parser.add_argument('probe_info',
                        help='Probe information in txt format')

    parser.add_argument('botloci', 
                        help='DArTag probes on the bottom strand')

    parser.add_argument('vcf_header',
                        help='VCF header')
    
    parser.add_argument('report',
                        help='DArTag read count CSV input file - after assigning fixed alleleIDs')

    args=parser.parse_args()

    snp_position = parse_probe_info(args.probe_info)
    print('Number of marker loci from probe design:', len(snp_position))

    ref_alt_read_counts, sample_line = collect_ref_alt_read_counts(args.report)
    print('Number of marker loci from MADC:', len(ref_alt_read_counts))
    
    bottom_loci = get_bottom_loci(args.botloci)

    generate_vcf(args.report, sample_line, ref_alt_read_counts, args.vcf_header, snp_position, bottom_loci)
