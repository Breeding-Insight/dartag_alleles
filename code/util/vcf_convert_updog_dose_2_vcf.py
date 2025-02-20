#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def get_updog_parameters(par):
    inp = open(par)
    header = inp.readline()
    line = inp.readline()
    par_dict = {}
    while line:
        line_array = line.strip().split(',')
        par_dict[line_array[0].strip('"')] = str(round(float(line_array[1]), 3))
        line = inp.readline()
    return (par_dict)    
    


def convert_updog_dose_to_gt(updog_dose, ploidy):
    import pandas as pd
    df_gt = pd.read_csv(updog_dose, index_col=[0])
    print(df_gt)
    # For updog, the genotype is the estimated reference allele dosage for a given individual at a given SNP.
    if ploidy == '2':
        df_gt = df_gt.replace([0.0, 1.0, 2.0], ['1/1', '0/1', '0/0'])
        df_gt = df_gt.fillna('./.')
    elif ploidy == '4':
        df_gt = df_gt.replace([0.0, 1.0, 2.0, 3.0, 4.0], ['1/1/1/1', '0/1/1/1', '0/0/1/1', '0/0/0/1', '0/0/0/0'])
        df_gt = df_gt.fillna('./.')
    print(df_gt)
    return (df_gt)


def convert_dosage2vcf(updog_dose, df_gt, vcf_readCount):
    import pandas as pd
    gt_markers = df_gt.index.tolist()
    outp = open(updog_dose.replace('.csv', '.vcf'), 'w')
    inp = open(vcf_readCount)
    line = inp.readline()
    read_count = {}
    while line:
        if line.startswith("##"):
            outp.write(line)
        elif line.startswith('#'):
            colnames = line.strip().split('\t')
        else:
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
    outp_col = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    outp_df = df_read_count.copy()
    outp_df = outp_df[outp_col]
    
    df_gt_col = df_gt.columns.tolist()
    for col in df_read_count.columns.tolist():
        if col in df_gt_col:
            outp_df[col] = df_gt[col] + ':' + df_read_count[col]
    print(outp_df)
    outp_df.to_csv(outp, sep="\t", mode='a', index=False)
    outp.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Convert GT calls to dosage")

    parser.add_argument('updog_dose', help='Dosage calls from updog')

    parser.add_argument('ploidy', help='Ploidy of the samples')

    parser.add_argument('vcf_readCount', help='vcf containing read counts from running hap2snpvcf')
    
    args=parser.parse_args()

    df_gt = convert_updog_dose_to_gt(args.updog_dose, args.ploidy)
    #print('updog_gt', str(dict(list(updog_gt.items())[0: 2])), len(updog_gt))

    convert_dosage2vcf(args.updog_dose, df_gt, args.vcf_readCount)
