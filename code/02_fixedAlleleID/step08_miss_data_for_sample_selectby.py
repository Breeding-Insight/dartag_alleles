#!/usr/bin/python3
# Parse intermediate files generated duing processing DArTag report (Samples with >=95% missing data were removed)
# e.g., DBlue21-6097_Counts_missing_allele_discovery_rename_miss_sample.csv


import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

def calculate_mean_missing_data_by_group(miss_data, miss_marker_rate, selectby, column):
    inp = open(selectby, 'r', encoding="ISO-8859-1")
    line = inp.readline() # header
    line = inp.readline()
    miss_data_by_group = {}
    while line:
        line_array = line.strip().split('\t')
        if line_array[column - 1] != '':
            if line_array[column - 1] in miss_data_by_group:
                if line_array[1] in miss_marker_rate:
                    miss_data_by_group[line_array[column - 1]].append(int(miss_marker_rate[line_array[1]]))
                else:
                    print('Sample not present in missing data file: ', line_array[1])
            else:
                if line_array[1] in miss_marker_rate:
                    miss_data_by_group[line_array[column - 1]] = [int(miss_marker_rate[line_array[1]])]
                else:
                    print('Sample not present in missing data file: ', line_array[1])
        line = inp.readline()
    inp.close()

    import os
    import datetime
    outf = miss_data.replace('.csv', '_groups.csv')
    if os.path.exists(outf):
        outp = open(outf, 'a')
        outp.write('\n######\n\n')
    else:
        outp = open(outf, 'w')
        timenow = datetime.datetime.now()
        outp.write('Processed on ' + timenow.strftime('%Y-%m-%d %H:%M:%S') + '\n\n')

    import statistics
    for k in sorted(miss_data_by_group, key=lambda k: len(miss_data_by_group[k]), reverse=True):
        average = round(float(sum(miss_data_by_group[k])/len(miss_data_by_group[k])), 2)
        percentage = round(average/3000*100, 2)
        std = statistics.pstdev(miss_data_by_group[k])
        outp.write(k + ',' + str(len(miss_data_by_group[k])) + ',' + str(average) + ',' + str(percentage) + ',' + str(std) + '\n')
    outp.close()


def put_miss_data_in_dict(miss_data):
    inp = open(miss_data)
    line = inp.readline() # header
    line = inp.readline()
    miss_marker_rate = {}
    while line:
        line_array = line.strip().split(',')
        if line_array[0] not in miss_marker_rate:
            miss_marker_rate[line_array[0]] = line_array[1]
        else:
            print('This sample has duplicate name with other samples: ', line_array[0])
        line = inp.readline()
    return(miss_marker_rate)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Subsetting samples from missing allele count in CSV format")

    parser.add_argument('missing_data',
                        help='Missing marker rate in CSV format')

    parser.add_argument('selectby',
                        help='A CSV file containing sample passport information')

    parser.add_argument('column',
                        help='The column number of the trait to be classified from the sample passport info')

    args=parser.parse_args()

    miss_marker_rate = put_miss_data_in_dict(args.missing_data)

    calculate_mean_missing_data_by_group(args.missing_data, miss_marker_rate, args.selectby, int(args.column))

