#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def get_vcf_data_lines(vcf, first_vcf):
    inp = open(vcf)
    line = inp.readline()
    gt_dict = {}
    if first_vcf == 'true': # include snp info in df
        while line:
            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')
            elif not line.startswith('##'):
                # Chr01	85423	Chr01_000085423	A	G	.	.	DP=68459;ADS=65959,2500	DP:RA:AD
                line_array = line.strip().split('\t')
                markerID = line_array[0].replace('"', '') + '_' + line_array[1].zfill(9)
                gt_dict[markerID] = line_array
            else:
                pass
            line = inp.readline()
    else:
        while line:
            # Just get sample data
            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
            elif not line.startswith('##'):
                # Chr01	85423	Chr01_000085423	A	G	.	.	DP=68459;ADS=65959,2500	DP:RA:AD
                line_array = line.strip().split('\t')
                if len(samples) != len(line_array[9:]):
                    print('check', len(samples), len(line_array[9:]), line_array)
                markerID = line_array[0].replace('"', '') + '_' + line_array[1].zfill(9)
                gt_dict[markerID] = line_array[9:]
            else:
                pass
            line = inp.readline()
    inp.close()
    import pandas as pd
    df_sample_gt = pd.DataFrame.from_dict(gt_dict, orient='index', columns=samples)
    return(df_sample_gt)
    


def merge_vcf(vcf_list, out_vcf):
    import pandas as pd
    outp = open(out_vcf, 'w')
    vcf_list = vcf_list.split(',')
    index = 0
    while index < len(vcf_list):
        if index == 0:
            sample_df = get_vcf_data_lines(vcf_list[index], 'true')
            merge_df = sample_df.copy()
            inp = open(vcf_list[index])
            line = inp.readline()
            while line:
                if line.startswith('##'):
                    outp.write(line)
                else:
                    pass
                line = inp.readline()
        else:
            sample_df = get_vcf_data_lines(vcf_list[index], 'false')
            merge_df = pd.concat([merge_df, sample_df], join='inner', axis=1)
        index += 1

    merge_df.to_csv(outp, sep="\t", mode='a', index=False)
    outp.close()




if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Merge multiple VCF files into one")

    parser.add_argument('vcf_list', help='VCF files separated by comma')
    
    parser.add_argument('outf', help='Output VCF file name')

    args=parser.parse_args()

    merge_vcf(args.vcf_list, args.outf)
