#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
    
def get_updog_parameters(par):
    inp = open(par)
    header = inp.readline()
    line = inp.readline()
    par_dict = {}
    while line:
        line_array = line.strip().split(',')
        par_dict[line_array[0].strip('"')] = str(round(float(line_array[1]), 3))
        line = inp.readline()
    inp.close()
    return(par_dict)


def convert_updog_dose_to_gt(updog_dose, ploidy):
    df_gt = pd.read_csv(updog_dose, index_col=[0])
    all_markers_wDose = df_gt.index.tolist()
    # For updog, the genotype is the estimated reference allele dosage for a given individual at a given SNP.
    if ploidy == '2':
        df_gt = df_gt.replace([0.0, 1.0, 2.0], ['1/1', '0/1', '0/0'])
        df_gt = df_gt.fillna('./.')
    elif ploidy == '4':
        df_gt = df_gt.replace([0.0, 1.0, 2.0, 3.0, 4.0], ['1/1/1/1', '0/1/1/1', '0/0/1/1', '0/0/0/1', '0/0/0/0'])
        df_gt = df_gt.fillna('./.')
    elif ploidy == '6':
        df_gt = df_gt.replace([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], ['1/1/1/1/1/1', '0/1/1/1/1/1', '0/0/1/1/1/1', '0/0/0/1/1/1', '0/0/0/0/1/1', '0/0/0/0/0/1', '0/0/0/0/0/0'])
        df_gt = df_gt.fillna('./.')
    return(df_gt, all_markers_wDose)


def add_missing_markers_2_par(all_markers_wDose, par_dict):
    for i in all_markers_wDose:
        if i not in par_dict:
            par_dict[i] = 'NA'

    import pandas as pd
    par_df = pd.DataFrame.from_dict(par_dict, orient='index', columns=['par'])
    return(par_df)


def parse_probe_info(probe_info):
    probe = open(probe_info)
    header = probe.readline() # Header
    # MarkerName	TargetSequence	ReferenceGenome	Chrom	Pos	VariantAllelesDef	Required
    line = probe.readline()
    snp_position = {}
    while line:
        line_array = line.strip().split('\t')
        alleles = line_array[5].replace('[', '').replace(']', '').split('/')
        markerID = line_array[3] + '_' + line_array[4].zfill(9)
        # If there are indels in the markers
        if alleles[0] == '-':
            snp_position[line_array[0]] = line_array[3:5] + [markerID, '.', alleles[1]]
        elif alleles[1] == '-':
            snp_position[line_array[0]] = line_array[3:5] + [markerID, alleles[0], '.']
        else:
            snp_position[line_array[0]] = line_array[3:5] + [markerID] + alleles
        # snp_position = {'alfalfaRep2vsXJDY1_shared_2402219': ['chr8.1', '20667639', 'G', 'A'],...}
        line = probe.readline()
    probe.close()

    snp_df = pd.DataFrame.from_dict(snp_position, orient='index')
    snp_df.columns = ["#CHROM", "POS", "ID", "REF", "ALT"]
    snp_df['QUAL'] = '.'
    snp_df['FILTER'] = '.'
    print(snp_df)
    return(snp_df)


def convert_dosage2vcf(snp_df, updog_dose, new_vcf_header, df_gt, bias_df, od_df, pmc_df, ref_mat_df, size_mat_df):
    outp = open(updog_dose.replace('.csv', '.vcf'), 'w')
    inp = open(new_vcf_header)
    lines = inp.readlines()
    for line in lines:
        outp.write(line)
    inp.close()

    for marker in snp_df.index.tolist():
        if marker not in df_gt.index.tolist():
            snp_df = snp_df.drop(marker, axis=0)

    # DP=14675;ADS=11691,2984;
    size_df = size_mat_df.sum(axis=1).to_frame('DP')
    ref_df = ref_mat_df.sum(axis=1).to_frame('ref')
    alt = size_df['DP'] - ref_df['ref']
    alt_df = alt.to_frame('alt')
    snp_df['INFO'] = 'DP=' + size_df['DP'].astype(str) + ';ADS=' + ref_df['ref'].astype(str) + ',' + alt_df['alt'].astype(str) + ';BIAS=' + bias_df['par'].astype(str) + ';OD=' + od_df['par'].astype(str) + ';PMC=' + pmc_df['par'].astype(str)
    snp_df['FORMAT'] = 'GT:DP:RA:AD'
    
    alt_mat_df = size_mat_df - ref_mat_df
    for col in df_gt.columns.tolist():
        snp_df[col] = df_gt[col] + ':' + size_mat_df[col].astype(str) + ':' + ref_mat_df[col].astype(str) + ':' + ref_mat_df[col].astype(str) + ',' + alt_mat_df[col].astype(str)
    snp_df.to_csv(outp, sep="\t", mode='a', index=False)
    outp.close()




if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Convert GT calls to dosage")
    
    parser.add_argument('bias', help='Estimated bias for the SNP from updog')

    parser.add_argument('od', help='Estimated overdispersion for the SNP from updog')

    parser.add_argument('pmc', help='Proportion of mis-classified individuals from updog')

    parser.add_argument('updog_dose', help='Dosage calls from updog')

    parser.add_argument('ploidy', help='Ploidy of the samples')
    
    parser.add_argument('ref_mat',
                        help='Reference read depth in csv format after updog filtering')

    parser.add_argument('size_mat', help='Total read depth in csv format after updog filtering')

    parser.add_argument('probe_info', help='Probe information in txt format')

    parser.add_argument('new_vcf_header',
                        help='New vcf header with additional features added to the INFO field')

    args=parser.parse_args()
    
    bias_dict = get_updog_parameters(args.bias)
    print('bias', str(dict(list(bias_dict.items())[0: 5])), len(bias_dict))
    
    od_dict = get_updog_parameters(args.od)
    print('od', str(dict(list(od_dict.items())[0: 5])), len(od_dict))
    
    pmc_dict = get_updog_parameters(args.pmc)
    print('pmc', str(dict(list(pmc_dict.items())[0: 5])), len(pmc_dict))
    
    df_gt, all_markers_wDose = convert_updog_dose_to_gt(args.updog_dose, args.ploidy)

    # Add missing markers in parameter files
    bias_df = add_missing_markers_2_par(all_markers_wDose, bias_dict)
    print('# Number of markers with bias: ', len(bias_df.index.tolist()))

    od_df = add_missing_markers_2_par(all_markers_wDose, od_dict)
    print('# Number of markers with od: ', len(od_df.index.tolist()))

    pmc_df = add_missing_markers_2_par(all_markers_wDose, pmc_dict)
    print('# Number of markers with pmc: ', len(pmc_df.index.tolist()))

    ref_mat_df = pd.read_csv(args.ref_mat, index_col=0)
    print('# Number of markers with ref read depth: ', len(ref_mat_df.index.tolist()))

    size_mat_df = pd.read_csv(args.size_mat, index_col=0)
    print('# Number of markers with total read depth: ', len(size_mat_df.index.tolist()))

    snp_df = parse_probe_info(args.probe_info)
    
    convert_dosage2vcf(snp_df, args.updog_dose, args.new_vcf_header, df_gt, bias_df, od_df, pmc_df, ref_mat_df, size_mat_df)
