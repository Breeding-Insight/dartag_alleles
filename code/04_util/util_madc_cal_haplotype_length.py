#!/usr/bin/python3
# The raw MADC file

def cal_haplotype_length(report):
    import pandas as pd
    pd.options.mode.chained_assignment = None  # default='warn'
    df = pd.read_csv(report, index_col='AlleleID', skiprows=7)
    remove_cols = ['ClusterConsensusSequence', 'CallRate', 'OneRatioRef', 'OneRatioSnp', 'FreqHomRef', 'FreqHomSnp', 'FreqHets', 'PICRef', 'PICSnp', 'AvgPIC', 'AvgCountRef', 'AvgCountSnp', 'RatioAvgCountRefAvgCountSnp', 'readCountSum']
    for col in remove_cols:
        if col in df.columns:
            df = df.drop(columns=col)
        else:
            pass

    outp = open(report.replace(".csv", "_hapLen.csv"), 'w')
    import datetime
    now = datetime.datetime.now()
    outp.write('#Processed on ' + str(now) + '\n\n')
    outp.write('haplotype_length' + ',' + 'count' + '\n')
    df.fillna(0, inplace=True)
    
    # Get the frequency of each length
    length_freq = df['AlleleSequence'].str.len().value_counts()
    length_freq = length_freq.to_dict()
    print('haplotype_length', 'count')
    for key, value in length_freq.items():
        outp.write(str(key) + ',' + str(value) + '\n')
        print(str(key), str(value))



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Calculate the length of haplotypes")

    parser.add_argument('report',
                        help='Missing allele report with allele name reformatted and unique sample names')

    args=parser.parse_args()

    cal_haplotype_length(args.report)
