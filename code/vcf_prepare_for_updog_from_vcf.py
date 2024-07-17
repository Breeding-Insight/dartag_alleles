#!/usr/bin/python3

import pandas as pd
pd.options.mode.chained_assignment = None

def prepare_for_updog(vcf):
    inp = open(vcf)
    # An output file with read depth for the reference allele
    out_refmat = open(vcf.replace('.vcf', '_refmat.csv'), 'w')
    # An output file with total read depth 
    out_sizemat = open(vcf.replace('.vcf', '_sizemat.csv'), 'w')
    line = inp.readline()
    while line:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                line_array = line.strip().split('\t')
                out_refmat.write('AlleleID,' + ','.join(line_array[9:]) + '\n')
                out_sizemat.write('AlleleID,' + ','.join(line_array[9:]) + '\n')
        else:
            # DP:RA:AD	137:83:83,54
            line_array = line.strip().split('\t')
            out_refmat.write(line_array[0] + '_' + line_array[1].zfill(9))
            out_sizemat.write(line_array[0] + '_' + line_array[1].zfill(9))
            index = 9
            while index < len(line_array):
                ref = line_array[index].split(':')[1]
                size = line_array[index].split(':')[0]
                out_refmat.write(',' + ref)
                out_sizemat.write(',' + size)
                index += 1
            out_refmat.write('\n')
            out_sizemat.write('\n')
        line = inp.readline()
    inp.close()
    out_refmat.close()
    out_sizemat.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('report',
                        help='Reformatted DArTag report')
    
    args=parser.parse_args()

    prepare_for_updog(args.report)
