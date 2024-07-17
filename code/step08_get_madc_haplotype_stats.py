#!/usr/bin/python3
# The modified MADC file with unique allele and sample names


def get_allele_counts(report, df_filterSamples, hap_threshold):
    import pandas as pd
    grouped_df = df_filterSamples.groupby('CloneID')
    df_mono = pd.DataFrame()
    df_10plus = pd.DataFrame()
    alleleCnt = {}
    total_sample = len(df_filterSamples.columns) - 1
    outp_summary = open(report.replace(".csv", "_hapCnt_summary.csv"), 'w')
    import datetime
    now = datetime.datetime.now()
    outp_summary.write('#Processed on ' + str(now) + '\n')

    threshold_10samples = round(10.0 / total_sample * 100, 2)
    if threshold_10samples < 5:
        threshold = 100 - threshold_10samples
        outp_summary.write('# Use threshold of at least 10 samples for which each sample having a minimum of 2 reads per microhaplotype: ' + str(threshold) + '\n')
    else:
        threshold = 95
        outp_summary.write('# Total number of samples:' + str(total_sample) + '\n# Use threshold of at least 5% samples for which each sample having a minimum of 2 reads per microhaplotype\n')

    for marker, group in grouped_df:
        group = group.drop(['CloneID'], axis=1)
        group = group.astype(int)
        # Transpose group and dropping microhaps that has >=95% missing data (<1 read per microhap)
        group_t = group.T # marker names as column labels
        # AlleleID  Chr04_007050770|Ref_0001  Chr04_007050770|Alt_0002
        # 95-145                         285                         0 
        for column in group_t:
            columnSeriesObj = len(group_t[group_t[column] < 2])
            missing = columnSeriesObj/total_sample * 100
            if missing >= threshold:
                # Note that this will drop Ref and Alt too
                group.drop(column, inplace=True)
            else:
                pass

        if len(group.index) == 1:
            df_mono = pd.concat([df_mono, group])
        elif len(group.index) >= hap_threshold:
            df_10plus = pd.concat([df_10plus, group])
        else:
            pass

        if len(group.index) not in alleleCnt:
            alleleCnt[len(group.index)] = 1
        else:
            alleleCnt[len(group.index)] += 1

    outp_mono = report.replace('.csv', '_monoLoci.csv')
    df_mono.to_csv(outp_mono)
    suffix = "_" + str(hap_threshold) + "plusHaps.csv"
    outp_10plus = report.replace(".csv", suffix)
    df_10plus.to_csv(outp_10plus)

    for i in range(0, 10):
        if i not in alleleCnt:
            alleleCnt[i] = 0
    ten_plus = 0
    outp_summary.write('#Microhaplotypes,#Marker_loci\n')
    
    for key in sorted(alleleCnt.keys()):
        print(key, alleleCnt[key])
        if key <= 10:
            outp_summary.write(str(key) + ',' + str(alleleCnt[key]) + '\n')
        else:
            ten_plus += alleleCnt[key]
    outp_summary.write(str('>10') + ',' + str(ten_plus) + '\n')


def get_marker_missing_rate_and_filter(df, df_sum):
    ## Missing markers per sample
    sample_miss_marker_rate = {}
    marker_count = len(df_sum.index)
    for column in df_sum:
        columnSeriesObj = len(df_sum[df_sum[column] < 10])
        missing_percent = round(columnSeriesObj / marker_count * 100, 2)
        sample_with_marker = 100.00 - missing_percent
        sample_miss_marker_rate[column] = [columnSeriesObj, missing_percent, sample_with_marker]

    # For a single sample, if less than i% of markers with missing data
    print('# Retaining a sample by requesting more than 5% of the markers with data (read count per marker locus >=10)')
    print('With_%marker\t#Samples')
    for i in range(0, 100, 5):
        sample_count = len([key for key, value in sample_miss_marker_rate.items() if float(value[2]) >= i])
        stdout = '>=' + str(i) + '%\t' + str(sample_count)
        print(stdout)
    
    # 100% data
    sample_count = len([key for key, value in sample_miss_marker_rate.items() if float(value[2]) == 100])
    stdout = '100%\t' + str(sample_count)
    print(stdout)
    
    # Retaining a sample by requesting more than 5% of the markers with data (read count per marker ≥10)
    remove_samples = [key for key, value in sample_miss_marker_rate.items() if float(value[1]) >= 95]
    keep_samples = [key for key, value in sample_miss_marker_rate.items() if float(value[1]) < 95]

    print('# Total samples: ', len(sample_miss_marker_rate), '\n# Number of samples with ≥95% missing data (remove from the report):', len(remove_samples))
    print('# Remove samples with high missing data:\n', remove_samples, '\n')
    df_filterSamples = df[['CloneID'] + keep_samples]
    return (df_filterSamples)


def convert_read_count_to_df_and_preprocess_df(report):
    import pandas as pd
    pd.options.mode.chained_assignment = None  # default='warn'
    df = pd.read_csv(report, index_col='AlleleID')
    remove_cols = ['ClusterConsensusSequence', 'CallRate', 'OneRatioRef', 'OneRatioSnp', 'FreqHomRef', 'FreqHomSnp', 'FreqHets', 'PICRef', 'PICSnp', 'AvgPIC', 'AvgCountRef', 'AvgCountSnp', 'RatioAvgCountRefAvgCountSnp', 'readCountSum']
    for col in remove_cols:
        if col in df.columns:
            df = df.drop(columns=col)
        else:
            pass

    df.fillna(0, inplace=True)
    grouped_df = df.groupby('CloneID')
    # Reading through the df and get the read count sum of each markers of all possible haplotypes
    df_sum = pd.DataFrame()
    for clone, group in grouped_df:
        # sample_sum is a Series
        sample_sum = group.sum(axis=0, skipna=True, numeric_only=True)
        sample_df = sample_sum.to_frame().transpose()
        sample_df['CloneID'] = clone
        df_sum = pd.concat([df_sum, sample_df])
    df_sum = df_sum.set_index('CloneID')
    return (df, df_sum)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('report',
                        help='Missing allele report with allele name reformatted and unique sample names')
    
    parser.add_argument('hap_threshold',
                        help='Minimum number of microhaplotypes per marker loci to be written to an output file for later filtering')

    args=parser.parse_args()

    df, df_sum = convert_read_count_to_df_and_preprocess_df(args.report)

    df_filterSamples = get_marker_missing_rate_and_filter(df, df_sum)
    
    get_allele_counts(args.report, df_filterSamples, int(args.hap_threshold))
