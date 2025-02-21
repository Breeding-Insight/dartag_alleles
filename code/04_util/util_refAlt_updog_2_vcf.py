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


def convert_dosage2vcf(updog_dose, new_vcf_header, df_gt, vcf_readCount, hh_df, bias_df, od_df, pmc_df):
    gt_markers = df_gt.index.tolist()
    outp = open(updog_dose.replace('.csv', '.vcf'), 'w')
    inp = open(new_vcf_header)
    lines = inp.readlines()
    for line in lines:
        outp.write(line)
    inp.close()
    
    inp = open(vcf_readCount)
    line = inp.readline()
    read_count = {}
    while line:
        if line.startswith('#CHROM'):
            colnames = line.strip().split('\t')
        elif not line.startswith('##'):
            # Chr01	85423	Chr01_000085423	A	G	.	.	DP=68459;ADS=65959,2500	DP:RA:AD
            line_array = line.strip().split('\t')
            markerID = line_array[0].replace('"', '') + '_' + line_array[1].zfill(9)
            if markerID in gt_markers:
                read_count[markerID] = line_array
            else:
                pass
        line = inp.readline()
    inp.close()
    
    df_read_count = pd.DataFrame.from_dict(read_count, orient='index', columns=colnames)
    print('Read count vcf:', df_read_count[df_read_count.index.duplicated()])
    print('Updog GT: ', df_gt[df_gt.index.duplicated()])
    outp_col = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    outp_df = df_read_count.copy()
    outp_df = outp_df[outp_col]
    outp_df['INFO'] = outp_df['INFO'] + ';HH=' + hh_df['par'].astype(str) + ';BIAS=' + bias_df['par'].astype(str) + ';OD=' + od_df['par'].astype(str) + ';PMC=' + pmc_df['par'].astype(str)
    print(outp_df['INFO'])

    df_gt_col = df_gt.columns.tolist()
    for col in df_read_count.columns.tolist():
        if col in df_gt_col:
            outp_df[col] = df_gt[col] + ':' + df_read_count[col]
    outp_df.to_csv(outp, sep="\t", mode='a', index=False)
    outp.close()




if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Convert GT calls to dosage")

    parser.add_argument('hh', help='Hind/He for markers from polyRAD')
    
    parser.add_argument('bias', help='Estimated bias for the SNP from updog')

    parser.add_argument('od', help='Estimated overdispersion for the SNP from updog')

    parser.add_argument('pmc', help='Proportion of mis-classified individuals from updog')

    parser.add_argument('updog_dose', help='Dosage calls from updog')

    parser.add_argument('ploidy', help='Ploidy of the samples')
    
    parser.add_argument('vcf_readCount',
                        help='vcf containing read counts from running hap2snpvcf')

    parser.add_argument('new_vcf_header',
                        help='New vcf header with additional features added to the INFO field')

    args=parser.parse_args()

    hh_dict = get_updog_parameters(args.hh)
    print('hh', str(dict(list(hh_dict.items())[0: 5])), len(hh_dict))
    
    bias_dict = get_updog_parameters(args.bias)
    print('bias', str(dict(list(bias_dict.items())[0: 5])), len(bias_dict))
    
    od_dict = get_updog_parameters(args.od)
    print('od', str(dict(list(od_dict.items())[0: 5])), len(od_dict))
    
    pmc_dict = get_updog_parameters(args.pmc)
    print('pmc', str(dict(list(pmc_dict.items())[0: 5])), len(pmc_dict))

    df_gt, all_markers_wDose = convert_updog_dose_to_gt(args.updog_dose, args.ploidy)

    # Add missing markers in parameter files
    hh_df = add_missing_markers_2_par(all_markers_wDose, hh_dict)
    print('# Number of markers with hh: ', len(hh_df.index.tolist()))

    bias_df = add_missing_markers_2_par(all_markers_wDose, bias_dict)
    print('# Number of markers with bias: ', len(bias_df.index.tolist()))

    od_df = add_missing_markers_2_par(all_markers_wDose, od_dict)
    print('# Number of markers with od: ', len(od_df.index.tolist()))

    pmc_df = add_missing_markers_2_par(all_markers_wDose, pmc_dict)
    print('# Number of markers with pmc: ', len(pmc_df.index.tolist()))
    
    convert_dosage2vcf(args.updog_dose, args.new_vcf_header, df_gt, args.vcf_readCount, hh_df, bias_df, od_df, pmc_df)
