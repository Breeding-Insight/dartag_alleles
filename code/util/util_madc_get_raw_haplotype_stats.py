#!/usr/bin/python3
# The raw MADC file

def get_haplotype_stats(report):
    import pandas as pd
    pd.options.mode.chained_assignment = None  # default='warn'
    df = pd.read_csv(report, index_col='AlleleID', skiprows=7)
    remove_cols = ['AlleleSequence', 'ClusterConsensusSequence', 'CallRate', 'OneRatioRef', 'OneRatioSnp', 'FreqHomRef', 'FreqHomSnp', 'FreqHets', 'PICRef', 'PICSnp', 'AvgPIC', 'AvgCountRef', 'AvgCountSnp', 'RatioAvgCountRefAvgCountSnp', 'readCountSum']
    for col in remove_cols:
        if col in df.columns:
            df = df.drop(columns=col)
        else:
            pass

    outp_summary = open(report.replace(".csv", "_hapCnt_summary.csv"), 'w')
    import datetime
    now = datetime.datetime.now()
    outp_summary.write('#Processed on ' + str(now) + '\n')
    df.fillna(0, inplace=True)
    grouped_df = df.groupby('CloneID')
    alleleCnt = {}
    for marker, group in grouped_df:
        group = group.drop(['CloneID'], axis=1)
        group = group.astype(int)
    
        if len(group.index) not in alleleCnt:
            alleleCnt[len(group.index)] = 1
        else:
            alleleCnt[len(group.index)] += 1

    for i in range(0, 10):
        if i not in alleleCnt:
            alleleCnt[i] = 0
    ten_plus = 0
    outp_summary.write('#Microhaplotypes,#Marker_loci\n')

    for key in sorted(alleleCnt.keys()):
        print(key, alleleCnt[key])
        if key < 10:
            outp_summary.write(str(key) + ',' + str(alleleCnt[key]) + '\n')
        else:
            ten_plus += alleleCnt[key]
    outp_summary.write(str('>=10') + ',' + str(ten_plus) + '\n')



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('report',
                        help='Missing allele report with allele name reformatted and unique sample names')

    args=parser.parse_args()

    get_haplotype_stats(args.report)
