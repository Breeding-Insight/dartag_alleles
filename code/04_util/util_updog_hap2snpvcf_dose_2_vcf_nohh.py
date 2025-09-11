#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def get_updog_parameters(par):
    inp = open(par)
    header = inp.readline()
    line = inp.readline()
    par_dict = {}
    while line:
        line_array = line.strip().split(',')
        par_dict[line_array[0].strip('"')] = str(round(float(line_array[1]),3))
        line = inp.readline()
    return(par_dict)


def determine_genotype(dose, ploidy):
    # For updog, the genotype is the estimated reference allele dosage for a given individual at a given SNP.
    # Calculate number of reference alleles (0's) and alternative alleles (1's)
    if dose.isdigit():
        ref_alleles = int(dose)  # number of 0's
        alt_alleles = int(ploidy) - int(dose)  # number of 1's
        # Create genotype string
        gt = '/'.join(['0'] * ref_alleles + ['1'] * alt_alleles)
    else:
        gt = '/'.join(['.'] * int(ploidy))  # missing data
        print('Invalid dose:', dose)
        print('Mark it as missing data:', gt)
    return(gt)


def convert_updog_dose_to_gt(updog_dose, ploidy):
    inp = open(updog_dose)
    header = inp.readline()
    line = inp.readline()
    updog_gt = {}
    while line:
        line_array = line.strip().split(',')
        updog_gt[line_array[0].strip('"')] = []
        index = 1
        while index < len(line_array):
            gt = determine_genotype(line_array[index], ploidy)
            updog_gt[line_array[0].strip('"')].append(gt)
            index += 1
        line = inp.readline()
    return(updog_gt)


def convert_dosage2vcf(updog_dose, new_vcf_header, vcf_readCount, bias_dict, od_dict, pmc_dict, updog_gt, ploidy, miss_cutoff):
    outp = open(updog_dose.replace('.csv', '.vcf'), 'w')
    inp = open(new_vcf_header)
    lines = inp.readlines()
    for line in lines:
        outp.write(line)
    inp.close()
    
    inp = open(vcf_readCount)
    line = inp.readline()
    cnt = 0
    while line:
        if line.startswith("#CHROM"):
            outp.write(line)
        elif not line.startswith('#'):
            # Chr01	85423	Chr01_000085423	A	G	.	.	DP=68459;ADS=65959,2500	DP:RA:AD
            line_array = line.strip().split()
            if line_array[2] in updog_gt:
                outp.write('\t'.join(line_array[:8]) + ';BIAS=' + bias_dict[line_array[2]] + ';OD=' + od_dict[line_array[2]] + ';PMC=' + pmc_dict[line_array[2]] + '\tGT:DP:RA:AD')
                vcf_index = 9
                while vcf_index < len(line_array):
                    dp = line_array[vcf_index].split(':')[0]
                    if int(dp) < int(miss_cutoff):
                        outp.write('\t' + '/'.join(['.'] * int(ploidy)) + ':' + line_array[vcf_index])
                        cnt += 1
                    else:
                        outp.write('\t' + updog_gt[line_array[2]][vcf_index - 9] + ':' + line_array[vcf_index])
                    vcf_index += 1
                outp.write('\n')
            else:
                pass
        else:
            pass
        line = inp.readline()
    print('Number of data points/genotypes converted to missing data if a genotype has less than', miss_cutoff + ' reads:', cnt)
    inp.close()
    outp.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Convert GT calls to dosage")
    
    parser.add_argument('bias', help='Estimated bias for the SNP from updog')

    parser.add_argument('od', help='Estimated overdispersion for the SNP from updog')

    parser.add_argument('pmc', help='Proportion of mis-classified individuals from updog')

    parser.add_argument('updog_dose', help='Dosage calls from updog')

    parser.add_argument('ploidy', help='Ploidy of the samples')
    
    parser.add_argument('vcf_readCount',
                        help='vcf containing read counts from running hap2snpvcf')

    parser.add_argument('new_vcf_header',
                        help='New vcf header with additional features added to the INFO field')
    
    parser.add_argument('miss_cutoff',
                        help='Minimum number of reads required to call a genotype')

    args=parser.parse_args()
    
    bias_dict = get_updog_parameters(args.bias)
    print('bias', str(dict(list(bias_dict.items())[0: 5])), len(bias_dict))
    
    od_dict = get_updog_parameters(args.od)
    print('od', str(dict(list(od_dict.items())[0: 5])), len(od_dict))
    
    pmc_dict = get_updog_parameters(args.pmc)
    print('pmc', str(dict(list(pmc_dict.items())[0: 5])), len(pmc_dict))
    
    updog_gt = convert_updog_dose_to_gt(args.updog_dose, args.ploidy)
    #print('updog_gt', str(dict(list(updog_gt.items())[0: 2])), len(updog_gt))
    
    convert_dosage2vcf(args.updog_dose, args.new_vcf_header, args.vcf_readCount, bias_dict, od_dict, pmc_dict, updog_gt, args.ploidy, args.miss_cutoff)
