#!/usr/bin/python3
# Parse MADC with haplotypes assigned fixed ID, i.e., after populating microhaplotype db
def update_madc_report(report, df, remove_markers, keep_samples):
    import pandas as pd
    outf = report.replace('.csv', '_filter_miss.csv')
    haplotypes_removed = 0
    for index in df.index.tolist():
        index_array = index.split('|')
        if index_array[0] in remove_markers:
            df.drop(index, inplace=True)
            haplotypes_removed += 1
        else:
            pass
    print('\nNumber of marker loci removed:', len(remove_markers))
    print('Number of microhaplotypes removed:', haplotypes_removed)

    df.reset_index(inplace=True)
    df = df.rename(columns={'index': 'AlleleID'})
    use_cols = ['AlleleID', 'CloneID', 'AlleleSequence'] + keep_samples
    df_keep_samples = df[use_cols]
    df_keep_samples[keep_samples] = df_keep_samples[keep_samples].astype('Int64')
    df_keep_samples.to_csv(outf, index=False)


def generate_filtered_data(report, df_sum_filterSamples):
    ## #Samples with data
    df_sum_filterSamples_t = df_sum_filterSamples.T
    print('# Filtering at marker level')
    print('# Total markers: ', len(df_sum_filterSamples_t.columns))
    marker_miss_sample_rate = {}
    samples_with_data = {}
    sample_count = len(df_sum_filterSamples_t.index)
    for column in df_sum_filterSamples_t:
        columnSeriesObj = len(df_sum_filterSamples_t[df_sum_filterSamples_t[column] < 10])
        marker_missing_rate = round(float(columnSeriesObj)/sample_count * 100, 2)
        marker_with_sample_rate = 100.00 - marker_missing_rate
        # Add sample names with_data into a dict
        if marker_with_sample_rate in samples_with_data:
            samples_with_data[marker_with_sample_rate] = list(set(samples_with_data[marker_with_sample_rate] + df_sum_filterSamples_t[df_sum_filterSamples_t[column] >= 10].index.tolist()))
        else:
            samples_with_data[marker_with_sample_rate] = df_sum_filterSamples_t[df_sum_filterSamples_t[column] >= 10].index.tolist()
        marker_miss_sample_rate[column] = [columnSeriesObj, marker_missing_rate, marker_with_sample_rate]

    # Retaining a marker locus based on whichever is smaller below:
    #       * At least 10 samples, each having a minimum of 2 reads
    #       * At least 5% samples, each having a minimum of 2 reads
    #outp_marker_summary.write('Total markers,' + str(len(marker_miss_sample_rate)) + '\n')
    samples = len(df_sum_filterSamples_t.index)  # Total number of samples
    threshold_10samples = round(10.0/samples * 100, 2)
    if threshold_10samples < 5:
        threshold = 100 - threshold_10samples
    else:
        threshold = 95
    remove_markers = [key for key, value in marker_miss_sample_rate.items() if int(value[1]) >= threshold]
    keep_markers = [key for key, value in marker_miss_sample_rate.items() if int(value[1]) < threshold]
    print('==> Total markers: ', len(marker_miss_sample_rate))
    print('10 samples in %total samples: ', threshold_10samples)
    print('Threshold used in filtering markers with missing data: ', threshold)
    print('Number of markers with missing data (remove from the report):', len(remove_markers))
    print('Check marker missing data summary for markers removed.')
    df_sum_filterMarkers = df_sum_filterSamples_t[keep_markers]
    df_sum_filterMarkers = df_sum_filterMarkers.T
    outp_marker_filter = open(report.replace('.csv', '_miss_sample_marker_filter.csv'), 'w')
    df_sum_filterMarkers.to_csv(outp_marker_filter)
    return(remove_markers)

def get_marker_missing_data_rate(report, df_sum):
    ## All samples
    df_sum_t = df_sum.T
    print('# Get missing data at marker level')
    print('# Total markers: ', len(df_sum_t.columns))
    marker_miss_sample_rate = {}
    samples_with_data = {}
    sample_count = len(df_sum_t.index)
    outp_marker = open(report.replace('.csv', '_miss_marker.csv'), 'w')
    outp_marker.write('MarkerID,#Samples_missing,#Sample/Total, #Samples_w_data\n')
    for column in df_sum_t:
        columnSeriesObj = len(df_sum_t[df_sum_t[column] < 10])
        marker_missing_rate = round(float(columnSeriesObj) / sample_count * 100, 2)
        marker_with_sample_rate = 100.00 - marker_missing_rate
        # Add sample names with_data into a dict
        if marker_with_sample_rate in samples_with_data:
            samples_with_data[marker_with_sample_rate] = list(set(
                samples_with_data[marker_with_sample_rate] + df_sum_t[
                    df_sum_t[column] >= 10].index.tolist()))
        else:
            samples_with_data[marker_with_sample_rate] = df_sum_t[
                df_sum_t[column] >= 10].index.tolist()
        marker_miss_sample_rate[column] = [columnSeriesObj, marker_missing_rate, marker_with_sample_rate]
        
        outp_marker.write(str(column) + ',' + str(columnSeriesObj) + ',' + str(marker_missing_rate) + '\n')
    outp_marker.close()
    
    # Output table
    outp_marker_summary = open(report.replace('.csv', '_miss_marker_summary.csv'), 'w')
    # Get sample names within each with_data range 
    outp_rate = open(report.replace('.csv', '_miss_rate_withSamples.csv'), 'w')
    import datetime
    now = datetime.datetime.now()
    nowf = now.strftime("%Y-%m-%d %H:%M:%S")
    outp_marker_summary.write('## ' + nowf + '\n')
    outp_marker_summary.write('# For each marker locus, if >=i% of samples with data (>=10 reads)\n\n')
    outp_marker_summary.write('Present_in_%sample,#Markers,#Samples\n')
    # For a single marker, if less than i% of samples with missing data
    print('Present_in_%sample,#Markers,#Samples')
    for i in range(0, 100, 5):
        # Added this part on 2023.9.19
        if i == 0:
            unique_samples = df_sum.columns.tolist()
        else:
            samples = []
            for key in sorted(samples_with_data.keys()):
                if key <= i:
                    samples = samples + samples_with_data[key]
            unique_samples = list(set(samples))
            outp_rate.write('>=' + str(i) + '%,' + ','.join(unique_samples) + '\n')
        
        marker_count = len([key for key, value in marker_miss_sample_rate.items() if int(value[2]) >= i])
        outp_marker_summary.write('>=' + str(i) + '%,' + str(marker_count) + ',' + str(len(unique_samples)) + '\n')
        stdout = '>=' + str(i) + '%\t' + str(marker_count) + '\t' + str(len(unique_samples))
        print(stdout)
    # 100% data
    samples = []
    for key in sorted(samples_with_data.keys()):
        if key == 100:
            samples = samples + samples_with_data[key]
    unique_samples = list(set(samples))
    outp_rate.write('>=' + str(100) + '%,' + ','.join(unique_samples) + '\n')
    
    marker_count = len([key for key, value in marker_miss_sample_rate.items() if int(value[2]) == 100])
    outp_marker_summary.write('100%,' + str(marker_count) + ',' + str(len(unique_samples)) + '\n')
    stdout = '100%\t' + str(marker_count) + '\t' + str(len(unique_samples))
    print(stdout)


def get_sample_missing_data_rate_and_filter(report, df_sum):
    ## Missing markers per sample
    sample_miss_marker_rate = {}
    outp_sample = open(report.replace('.csv', '_miss_sample.csv'), 'w')
    outp_sample.write('Sample_ID,#Markers_missing,#Markers/Total\n')
    marker_count = len(df_sum.index)
    print('# Filtering at sample level')
    print('# Total samples:', len(df_sum.columns))
    for column in df_sum:
        columnSeriesObj = len(df_sum[df_sum[column] < 10])
        missing_percent = round(columnSeriesObj/marker_count*100, 2)
        sample_with_marker = 100.00 - missing_percent
        sample_miss_marker_rate[column] = [columnSeriesObj, missing_percent, sample_with_marker]
        outp_sample.write(column + ',' + str(columnSeriesObj) + ',' + str(missing_percent) + '\n')
    outp_sample.close()

    outp_sample_summary = open(report.replace('.csv', '_miss_sample_summary.csv'), 'w')
    import datetime
    now = datetime.datetime.now()
    nowf = now.strftime("%Y-%m-%d, %H:%M:%S")
    outp_sample_summary.write('## ' + nowf + '\n')
    outp_sample_summary.write('# Retaining a sample by requesting more than 5% of the markers with data (read count per marker >=10)\n\n')
    outp_sample_summary.write('With_%marker,#Samples\n')
    # For a single sample, if less than i% of markers with missing data
    print('# Retaining a sample by requesting more than 5% of the markers with data (read count per marker locus >=10)')
    print('With_%marker\t#Samples')
    for i in range(0, 100, 5):
        sample_count = len([key for key, value in sample_miss_marker_rate.items() if float(value[2]) >= i])
        outp_sample_summary.write('>=' + str(i) + '%,' + str(sample_count) + '\n')
        stdout = '>=' + str(i) + '%\t' + str(sample_count)
        print(stdout)

    # 100% data
    sample_count = len([key for key, value in sample_miss_marker_rate.items() if float(value[2]) == 100])
    outp_sample_summary.write('100%,' + str(sample_count) + '\n')
    stdout = '100%\t' + str(sample_count)
    print(stdout)
    
    # Retaining a sample by requesting more than 5% of the markers with data (read count per marker ≥10)
    remove_samples = [key for key, value in sample_miss_marker_rate.items() if float(value[1]) >= 95]
    keep_samples = [key for key, value in sample_miss_marker_rate.items() if float(value[1]) < 95]
    #outp_sample_summary.write('Total samples,' + str(len(sample_miss_marker_rate)) + '\n')
    outp_sample_summary.write('#Sample removed,' + str(len(remove_samples)) + '\n\n')
    outp_sample_summary.write('\nSamples removed\n' + '\n'.join(remove_samples) + '\n')

    print('==> Total samples: ', len(sample_miss_marker_rate), '\nNumber of samples with ≥95% missing data (remove from the report):', len(remove_samples))
    print(remove_samples, '\n')
    df_sum_filterSamples = df_sum[keep_samples]
    outp_sample_filter = report.replace('.csv', '_miss_sample_filter.csv')
    df_sum_filterSamples.to_csv(outp_sample_filter)
    return(df_sum_filterSamples, keep_samples)


def get_haplotype_read_count_sum(df):
    import pandas as pd
    # ['CloneID', 'AlleleSequence', 'Sample_01']
    total_sample = len(df.columns) - 2
    threshold_10samples = round(10.0 / total_sample * 100, 2)
    if threshold_10samples < 5:
        threshold = 100 - threshold_10samples
    else:
        threshold = 95

    # Transpose df and dropping microhaps that has missing data higher than the threshold (<2 read per microhap)
    df_t = df.T
    df_t.drop(['CloneID', 'AlleleSequence'], inplace=True)
    df_t = df_t.astype(int)
    # marker names as column labels
    # AlleleID  Chr04_007050770|Ref_0001  Chr04_007050770|Alt_0002
    # 95-145                         285                         0 
    drop_hap = []
    for column in df_t:
        columnSeriesObj = len(df_t[df_t[column] < 2])
        missing = columnSeriesObj / total_sample * 100
        if missing >= threshold:
            if 'RefMatch' in column or 'AltMatch' in column:
                df.drop(column, axis=0, inplace=True)
                drop_hap.append(column)
            else:
                pass  # Don't remove Ref and Alt
        else:
            pass
    print('# Number of microhaplotypes dropped due to missing data:', len(drop_hap))

    # Reading through the df and get the read count sum of each markers of all possible haplotypes
    grouped_df = df.groupby('CloneID')
    #print(grouped_df.get_group((list(grouped_df.groups)[0])))
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

    # Generate marker missing data rate
    get_marker_missing_data_rate(args.report, df_sum)

    '''
    Filter the data at SAMPLE level
    1. Generate stats for SAMPLEs' missing data rate and output samples by requesting more than 5% of the markers with data (read count per marker ≥10))
    2. Generate stats for MARKERs' missing data rate and output markers by requesting at least 10 samples or 5% of total samples, each sample having a minimum of 10 reads per marker locus
    3. Generate an MADC file with filtered samples and markers
    '''

    # Generate sample missing data rate and filter sample with >=95% missing markers
    df_sum_filterSamples, keep_samples = get_sample_missing_data_rate_and_filter(args.report, df_sum)

    # From the output from above, generate stats for marker missing rate
    remove_markers = generate_filtered_data(args.report, df_sum_filterSamples)

    update_madc_report(args.report, df, remove_markers, keep_samples)
