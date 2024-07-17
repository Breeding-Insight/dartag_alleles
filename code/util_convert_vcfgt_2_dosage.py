#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def convert_gf2dosage(vcf):
    inp = open(vcf)
    line = inp.readline()
    outp = open(vcf.replace('.vcf', '_dose.txt'), 'w')
    while line:
        if line.startswith('#CHROM'):
            line_array = line.strip().split('\t')
            outp.write('\t'.join(line_array[0:2]) + '\t' + '\t'.join(line_array[3:5]) + '\t' + '\t'.join(line_array[9:]) + '\n')
        elif line.startswith('##'):
            pass
        else:
            line_array = line.strip().split('\t')
            outp.write('\t'.join(line_array[0:2]) + '\t' + '\t'.join(line_array[3:5]))
            index = 9
            while index < len(line_array):
                gt_array = line_array[index].split(':')
                ref_cnt = gt_array[0].count('0')
                outp.write('\t' + str(ref_cnt))
                index += 1
            outp.write('\n')
        line = inp.readline()
    inp.close()
    outp.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Convert GT calls to dosage")

    parser.add_argument('vcf',
                        help='vcf from running genotype call using polyRAD')


    args=parser.parse_args()

    convert_gf2dosage(args.vcf)
