#!/usr/bin/python3

def collect_passport_data(passport):
    inp = open(passport, encoding='utf-8-sig')
    line = inp.readline() # header
    # ['Sample_ID', 'Breeder', ...]
    print('Header:', line)
    line = inp.readline()
    sample_lut = {}
    sample_count = 0
    while line:
        line_array = line.strip().split(',')
        sample_count += 1
        if line_array[1] not in sample_lut:
            sample_lut[line_array[1]] = ['AlleleID', 'CloneID', 'AlleleSequence', line_array[0].strip()]
        else:
            sample_lut[line_array[1]].append(line_array[0].strip())
        line = inp.readline()
    inp.close()
    print('  # Number of samples in the PASSPORT: ', sample_count)
    print('  # Breeders: ', list(sample_lut.keys()), '\n')
    return(sample_lut)



def split_samples_to_breeders(report, sample_lut):
    import pandas as pd
    df = pd.read_csv(report, header=None, low_memory=False)
    df.columns = df.iloc[7]
    data_col = df.iloc[0].isna().sum().sum()
    print(df)
    outf_summary = open(report.replace('.csv', '_split_summary.csv'), 'w')
    outf_summary.write('# Number of samples in the REPORT: ' + str(len(df.columns.tolist()) - data_col) + '\n')
    print('  # Number of samples in the report: ', len(df.columns.tolist()) - data_col)

    for breeder, sample_list in sample_lut.items():
        outf_summary.write(breeder + ':\n  original sample list: ' + str(len(sample_list)-3) + '\n')
        samples_not_in_report = sample_list.copy()
        df_sub = df
        for col in df.columns.tolist():
            if col not in sample_list:
                df_sub = df_sub.drop(columns=col)
            else:
                samples_not_in_report.remove(col)
        outf_summary.write('    # Number of samples not present in the report: ' + str(len(samples_not_in_report)) + '\n    ' + '\n    '.join(samples_not_in_report) + '\n')
        outf_summary.write('    # Number of samples with data extracted: ' + str(len(df_sub.columns.tolist()) - data_col) + '\n\n')
        
        print('  #', breeder)
        print('    # original sample list:', len(sample_list) - data_col)
        print('    #', len(samples_not_in_report), 'samples not present in report: ', samples_not_in_report)
        print('    # After cleaning:', len(df_sub.columns.tolist()) - data_col)
        outf = report.replace('.csv', '_' + breeder + '.csv')
        df_sub.to_csv(outf, index=False, header=False)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('passport',
                        help='A file containing the passport data of the samples in the report')

    parser.add_argument('report',
                        help='Reformated DArTag report, including fixed alleleIDs')

    args=parser.parse_args()

    sample_lut = collect_passport_data(args.passport)

    split_samples_to_breeders(args.report, sample_lut)
