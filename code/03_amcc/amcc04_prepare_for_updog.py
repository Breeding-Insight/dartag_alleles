#!/usr/bin/python3

import pandas as pd
pd.options.mode.chained_assignment = None

def prepare_for_updog(file):
    inp = open(file)
    outp = open(file.replace('.csv', '_refmat.csv'), 'w')
    outp2 = open(file.replace('.csv', '_sizemat.csv'), 'w')
    header = inp.readline()
    outp.write(header)
    outp2.write(header)
    line = inp.readline()
    while line:
        line_array = line.strip().split(',')
        next = inp.readline()
        next_array = next.strip().split(',')
        if line_array[0][:-4] == next_array[0][:-4]:
            outp.write(line_array[0][:-4])
            outp2.write(line_array[0][:-4])
            index = 1
            while index < len(line_array):
                outp.write(',' + str(int(line_array[index])))
                totalDepth = int(line_array[index]) + int(next_array[index])
                outp2.write(',' + str(totalDepth))
                index += 1
            outp.write('\n')
            outp2.write('\n')
        line = inp.readline()
    inp.close()
    outp.close()
    outp2.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('report',
                        help='Reformatted DArTag report')
    
    args=parser.parse_args()

    prepare_for_updog(args.report)
