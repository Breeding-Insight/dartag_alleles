#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def get_subset_samples_and_markers(file):
    import pandas as pd
    df = pd.read_csv(file, index_col=0)
    IDs = df.index.tolist()
    return(IDs)


def get_vcf_data_lines(vcf, outf):
    inp = open(vcf)
    line = inp.readline()
    outp = open(outf, 'w')
    gt_dict = {}
    while line:
        if line.startswith('##'):
            outp.write(line)
        elif line.startswith('#CHROM'):
            samples = line.strip().split('\t')
        else:
            # Chr01	85423	Chr01_000085423	A	G	.	.	DP=68459;ADS=65959,2500	DP:RA:AD
            line_array = line.strip().split('\t')
            markerID = line_array[0].replace('"', '') + '_' + line_array[1].zfill(9)
            gt_dict[markerID] = line_array
        line = inp.readline()
    inp.close()

    import pandas as pd
    df_sample_gt = pd.DataFrame.from_dict(gt_dict, orient='index', columns=samples)
    return(df_sample_gt)


def subset_samples_and_markers(df_sample_gt, samples, markers, outp):
    import pandas as pd
    out_col = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples
    vcf_samples = df_sample_gt.columns.tolist()
    
    # Calculate sum of read depth for reference and alternative alleles
    df_sample_gt_reads = df_sample_gt.copy(deep=True)
    for i in vcf_samples:
        if i not in samples:
            df_sample_gt_reads = df_sample_gt_reads.drop(i, axis=1) # drop sample/column
    df_sample_gt_reads_t = df_sample_gt_reads.transpose()
    df_sample_gt_reads_t_split = df_sample_gt_reads_t.apply(lambda x: x.str.split(':|,'))
    #                    Chr01_000084077     Chr01_000084078    Chr01_000084094    Chr01_000084095   
    # S_RL_84__3_A1   [271, 217, 217, 54]  [271, 217, 217, 54]  [271, 271, 271, 0]  [271, 271, 271, 0] 
    size = {}
    ref = {}
    alt = {}
    for column in df_sample_gt_reads_t_split:
        column_size = df_sample_gt_reads_t_split[column].str[0].astype(int).sum()
        size[column] = str(column_size)
        
        column_ref = df_sample_gt_reads_t_split[column].str[2].astype(int).sum()
        ref[column] = str(column_ref)

        column_alt = df_sample_gt_reads_t_split[column].str[-1].astype(int).sum()
        alt[column] = str(column_alt)
    
    size_df = pd.DataFrame.from_dict(size, orient='index', columns=['Size'])
    ref_df = pd.DataFrame.from_dict(ref, orient='index', columns=['Ref'])
    alt_df = pd.DataFrame.from_dict(alt, orient='index', columns=['Alt'])
    # Drop samples not present in sample list
    for i in vcf_samples:
        if i not in out_col:
            df_sample_gt = df_sample_gt.drop(i, axis=1) # drop sample/column
    
    # Drop markers not present in marker list
    vcf_markers = df_sample_gt.index.tolist()
    for j in vcf_markers:
        if j not in markers:
            df_sample_gt = df_sample_gt.drop(j, axis=0)

    # DP=39418;ADS=27713,11705
    df_sample_gt['INFO'] = 'DP=' + size_df['Size'] + ';ADS=' + ref_df['Ref'] + ',' + alt_df['Alt']
    print('# Updated read counts\n', df_sample_gt)
    df_sample_gt.to_csv(outp, sep="\t", mode='a', index=False)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Merge multiple VCF files into one")

    parser.add_argument('keep_samples', help='Samples to keep')
    
    parser.add_argument('keep_markers', help='Markers to keep')

    parser.add_argument('vcf', help='VCF file with read counts')

    args=parser.parse_args()

    samples = get_subset_samples_and_markers(args.keep_samples)
    print('# Number of samples to keep:', len(samples))
    
    markers = get_subset_samples_and_markers(args.keep_markers)
    print('# Number of markers to keep:', len(markers))

    outf = args.vcf.replace('.vcf', '_polyRADcleaned.vcf')
    df_sample_gt = get_vcf_data_lines(args.vcf, outf)
    
    subset_samples_and_markers(df_sample_gt, samples, markers, outf)
