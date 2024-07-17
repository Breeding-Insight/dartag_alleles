#!/usr/bin/python3
# Parse the original file that DArT provided


def update_madc_report(report, df, remove_markers, keep_samples):
    outf = report.replace('.csv', '_filter_miss.csv')
    haplotypes_removed = 0
    for index in df.index.tolist():
        index_array = index.split('|')
        if index_array[0] in remove_markers:
            df.drop(index, inplace=True)
            haplotypes_removed += 1
        else:
            pass
    print('Number of haplotypes/alleles removed (>=95% missing data):', haplotypes_removed)

    df.reset_index(inplace=True)
    df = df.rename(columns={'index': 'AlleleID'})
    use_cols = ['AlleleID'] + keep_samples
    print('Before removing samples with missing data: ', len(df.columns.tolist())-2)
    df_keep_samples = df[use_cols]
    print('After removing samples with missing data: ', len(df_keep_samples.columns.tolist())-1)
    df_keep_samples.to_csv(outf, index=False)


    
def get_marker_missing_sample_rate_and_filter(report, df_sum_filter_samples):
    ## Missing samples per marker
    df_sum_filter_samples_t = df_sum_filter_samples.T
    marker_miss_sample_rate = {}
    sample_count = len(df_sum_filter_samples_t.index)
    outp_marker = open(report.replace('.csv', '_miss_marker.csv'), 'w')
    outp_marker.write('MarkerID,#Samples_missing,#Sample/Total\n')
    for column in df_sum_filter_samples_t:
        columnSeriesObj = len(df_sum_filter_samples_t[df_sum_filter_samples_t[column] < 10])
        marker_missing_rate = round(float(columnSeriesObj)/sample_count * 100, 2)
        marker_miss_sample_rate[column] = [columnSeriesObj, marker_missing_rate]
        outp_marker.write(column + ',' + str(columnSeriesObj) + ',' + str(marker_missing_rate) + '\n')
    outp_marker.close()

    # Output table
    outp_marker_summary = open(report.replace('.csv', '_miss_marker_summary.csv'), 'w')
    outp_marker_summary.write('# For each marker, if less than or equal to i% of samples with missing data (<10 reads)\n\n')
    outp_marker_summary.write('Missing_samples/marker,#Markers\n')
    # For a single marker, if less than i% of samples with missing data
    print('Missing_samples/marker,#Markers')
    for i in range(0, 101, 5):
        marker_count = len([key for key, value in marker_miss_sample_rate.items() if int(value[1]) >= i])
        outp_marker_summary.write(str(i) + '%,' + str(marker_count) + '\n')
        stdout = '>=' + str(i) + '%\t' + str(marker_count)
        print(stdout)

    # Remove markers with >=95% missing data
    remove_markers = [key for key, value in marker_miss_sample_rate.items() if int(value[1]) >= 95]
    keep_markers = [key for key, value in marker_miss_sample_rate.items() if int(value[1]) < 95]
    outp_marker_summary.write('==> Total markers: ' + str(len(marker_miss_sample_rate)) + '\n')
    outp_marker_summary.write('==> Number of samples with >=95% missing data (remove from the report): ' + str(len(remove_markers)) + '\n\n')
    outp_marker_summary.write('# Markers with >=95% missing data:\n' + '\n'.join(remove_markers) + '\n')
    print('==> Total markers: ', len(marker_miss_sample_rate), '\nNumber of markers with >=95% missing data (remove from the report):', len(remove_markers), '\n')
    df_sum_filter_markers = df_sum_filter_samples_t[keep_markers]
    df_sum_filter_markers = df_sum_filter_markers.T
    outp_marker_filter = open(report.replace('.csv', '_miss_sample_marker_filter.csv'), 'w')
    df_sum_filter_markers.to_csv(outp_marker_filter)
    return(df_sum_filter_markers, remove_markers)


def get_sample_missing_marker_rate_and_filter(report, df_sum):
    ## Missing markers per sample
    sample_miss_marker_rate = {}
    outp_sample = open(report.replace('.csv', '_miss_sample.csv'), 'w')
    outp_sample.write('SampleID,#Markers_missing,#Markers/Total\n')
    for column in df_sum:
        columnSeriesObj = len(df_sum[df_sum[column] < 10])
        missing_percent = round(columnSeriesObj/3000*100, 2)
        sample_miss_marker_rate[column] = [columnSeriesObj, missing_percent]
        outp_sample.write(column + ',' + str(columnSeriesObj) + ',' + str(missing_percent) + '\n')
    outp_sample.close()

    outp_sample_summary = open(report.replace('.csv', '_miss_sample_summary.csv'), 'w')
    outp_sample_summary.write('For each sample, determine the number of markers with missing data (<10 reads)\n\n')
    outp_sample_summary.write('Missing_markers/sample,#Samples\n')
    # For a single sample, if less than i% of markers with missing data
    print('Missing_markers/sample\t#Samples')
    for i in range(0, 101, 5):
        sample_count = len([key for key, value in sample_miss_marker_rate.items() if float(value[1]) >= i])
        outp_sample_summary.write('>=' + str(i) + '%,' + str(sample_count) + '\n')
        stdout = '>=' + str(i) + '%\t' + str(sample_count)
        print(stdout)

    # Remove samples with >=95% missing data
    remove_samples = [key for key, value in sample_miss_marker_rate.items() if float(value[1]) >= 95]
    keep_samples = [key for key, value in sample_miss_marker_rate.items() if float(value[1]) < 95]
    outp_sample_summary.write('==> Total samples: ' + str(len(sample_miss_marker_rate)) + '\n')
    outp_sample_summary.write('==> Number of samples with >=95% missing data (remove from the report): ' + str(len(remove_samples)) + '\n\n')
    outp_sample_summary.write('\n# Samples with >=95% missing data:\n' + '\n'.join(remove_samples) + '\n')

    print('==> Total samples: ', len(sample_miss_marker_rate), '\nNumber of samples with >=95% missing data (remove from the report):', len(remove_samples), '\n')
    df_sum_samples_keep = df_sum[keep_samples]
    outp_sample_filter = report.replace('.csv', '_miss_sample_filter.csv')
    df_sum_samples_keep.to_csv(outp_sample_filter)
    return(df_sum_samples_keep, keep_samples)


def get_haplotype_read_count_sum(df):
    import pandas as pd
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
    return(df_sum)


def convert_read_count_to_df_and_preprocess_df(report):
    with open(report) as inp:
        header = inp.readline()
        
    import pandas as pd
    pd.options.mode.chained_assignment = None  # default='warn'
    if 'AlleleID' in header:
        df = pd.read_csv(report)
    else:
        df = pd.read_csv(report, skiprows=7)
        
    df['CloneID'] = df['AlleleID'].str.split('|', expand=True)[0]
    df = df.set_index('AlleleID')
    remove_cols = ['ClusterConsensusSequence', 'CallRate', 'OneRatioRef', 'OneRatioSnp', 'FreqHomRef', 'FreqHomSnp', 'FreqHets', 'PICRef', 'PICSnp', 'AvgPIC', 'AvgCountRef', 'AvgCountSnp', 'RatioAvgCountRefAvgCountSnp', 'readCountSum']
    for col in remove_cols:
        if col in df.columns:
            df = df.drop(columns=col)
        else:
            pass
    df.fillna(0, inplace=True)
    return(df)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Get missing data stats")

    parser.add_argument('report',
                        help='Missing allele count report in CSV format')

    args=parser.parse_args()

    df = convert_read_count_to_df_and_preprocess_df(args.report)
    # Columns: Index(['USDAMSP1_A01', 'USDAMSP1_C01', 'USDAMSP1_D01', 'USDAMSP1_E01',......'USDAMSP5_F12', 'CloneID'], dtype='object', length=470)
    # Rows: Index(['VaccDscaff11_000042737|Ref_0001', ......, 'VaccDscaff7_041659267|Alt_0002'],dtype='object', name='AlleleID', length=12456)

    df_sum = get_haplotype_read_count_sum(df)

    '''
    Filter the data at SAMPLE level
    1. Generate stats for SAMPLEs' missing data rate and output samples passing filter (<=95% missing rate)
    2. Generate stats for MARKERs' missing data rate and output markers passing filter (<=95% missing rate)
    3. Generate correlation matrix for the filtered data
    '''

    df_sum_samples_keep, keep_samples = get_sample_missing_marker_rate_and_filter(args.report, df_sum)

    # From the output from above, generate stats for marker missing rate
    df_sum_filter_markers, remove_markers = get_marker_missing_sample_rate_and_filter(args.report, df_sum_samples_keep)

    update_madc_report(args.report, df, remove_markers, keep_samples)
