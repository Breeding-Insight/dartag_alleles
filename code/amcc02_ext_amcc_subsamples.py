#!/usr/bin/python3

import pandas as pd
pd.options.mode.chained_assignment = None

def get_sample_names(file):
    inp = open(file)
    header = inp.readline()
    line = inp.readline()
    sample_names = []
    dup_names = []
    while line:
        line_array = line.strip().split(',')
        sample = line_array[0].strip()
        if sample not in sample_names:
            sample_names.append(sample)
        else:
            #print('Sample name occurred before:', line_array)
            dup_names.append(line)
        line = inp.readline()
    inp.close()
    if len(dup_names) != 0:
        outp = open(file.replace('.csv', '_dupSamples.csv'), 'w')
        for i in dup_names:
            outp.write(i)
        outp.close()
        print("Number of samples with duplicate names in 'Sample list': ", len(dup_names))
        print("Samples with duplicate names in 'Sample list': ", dup_names)
        print("Check the output file '_dupSamples.csv' for those duplicate sample names.")
    else:
        print("No sample with duplicate names")
    return(sample_names)


def get_subsamples(file, sample_names, report):
    df = pd.read_csv(report)
    df = df.set_index('AlleleID')
    columns = []
    noData = []
    for sample in sample_names:
        if sample in df.columns:
            columns.append(sample)
        else:
            noData.append(sample)
            
    if len(noData) != 0:
        outp = open(file.replace('.csv', '_noData.csv'), 'w')
        print("Number of samples not present in DArTag report: ", len(noData))
        print("Check the output file 'noData.csv' for samples not present in DArTag report.")
        for i in noData:
            outp.write(i + '\n')
        outp.close()
    else:
        print("Number of samples NOT present in DArTag reports: 0")
        
    df_sub = df[columns]
    df_sub.fillna(0, inplace=True)
    df_sub.reset_index(inplace=True)
    df_sub = df_sub.rename(columns={'index': 'AlleleID'})
    
    outf = file.replace('.csv', '_amcc.csv')
    df_sub.to_csv(outf, index=False)




if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('samples',
                        help='A file containing sample names')

    parser.add_argument('report',
                        help='Reformatted DArTag report')
    
    args=parser.parse_args()

    sample_names = get_sample_names(args.samples)
    
    get_subsamples(args.samples, sample_names, args.report)
