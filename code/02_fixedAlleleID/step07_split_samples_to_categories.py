#!/usr/bin/python3

def collect_passport_data(passport, category_col_name, sample_number_threshold):
    inp = open(passport, encoding='utf-8-sig')
    line = inp.readline() # header
    # ['Sample_ID', 'Breeder', 'Taxon', 'Source', 'Study_ID', 'Plate_ID', 'Well_ID', 'Plant_ID', 'IVP', 'IVNO', 'ng/ul ', 'G12/H12 Blanks', 'SampleName length QC']
    line_array = line.strip().split(',')
    if category_col_name not in line_array:
        print('ERROR: category_col_name "', category_col_name, '" not found in passport header:', line_array)
        exit()
    else:
        category_col_index = line_array.index(category_col_name)

    line = inp.readline()
    sample_lut = {}
    sample_count = 0
    while line:
        line_array = line.strip().split(',')
        sample_count += 1
        if line_array[category_col_index] not in sample_lut:
            sample_lut[line_array[category_col_index]] = ['AlleleID', 'CloneID', 'AlleleSequence', line_array[0].strip()]
        else:
            sample_lut[line_array[category_col_index]].append(line_array[0].strip())
        line = inp.readline()
    
    # Check if each category contains at least 10 samples
    category_to_remove = []
    sample_lut['Others'] = ['AlleleID', 'CloneID', 'AlleleSequence']
    for category, sample_list in sample_lut.items():
        if len(sample_list) < int(sample_number_threshold) + 3: # 3 fixed columns
            print('Warning: category', category, 'contains less than', sample_number_threshold, 'samples:', len(sample_list) - 3)
            sample_lut['Others'] = sample_lut['Others'] + sample_list[3:]
            category_to_remove.append(category)
    inp.close()
    print('  # Number of samples in the PASSPORT: ', sample_count)
    print('  # ', category_col_name, list(sample_lut.keys()), '\n')
    return(sample_lut, category_to_remove)
    

def split_samples_to_breeders(report, sample_lut, category_to_remove):
    import pandas as pd
    df = pd.read_csv(report)
    outf_summary = open(report.replace('.csv', '_split_summary.csv'), 'w')
    outf_summary.write('# Number of samples in the REPORT: ' + str(len(df.columns.tolist()) - 3) + '\n')
    print('  # Number of samples in the report: ', len(df.columns.tolist()) - 3)
    for category, sample_list in sample_lut.items():
        outf_summary.write(category + ':\n  original sample list: ' + str(len(sample_list)-3) + '\n')
        samples_not_in_report = sample_list.copy()
        if category not in category_to_remove:
            df_sub = df
            for col in df.columns.tolist():
                if col not in sample_list:
                    df_sub = df_sub.drop(columns=col)
                else:
                    samples_not_in_report.remove(col)
            
            # Only retain haplotypes with more than 10 reads
           # df_sub_remove_miss = df_sub[df_sub.sum(axis=1) > 10]
            outf_summary.write('    # Number of samples not present in the report: ' + str(len(samples_not_in_report)) + '\n    ' + '\n    '.join(samples_not_in_report) + '\n')
            outf_summary.write('    # Number of samples with data extracted: ' + str(len(df_sub.columns.tolist()) - 3) + '\n\n')
            
            print('  #', category)
            print('    # original sample list:', len(sample_list) - 3)
            print('    #', len(samples_not_in_report), 'samples not present in report: ', samples_not_in_report)
            print('    # After cleaning:', len(df_sub.columns.tolist()) - 3)
            outf = report.replace('.csv', '_' + category.replace(' ', '_') + '.csv')
            df_sub.to_csv(outf, index=False)
        else:
            outf_summary.write('    # Category merged into "Others" due to insufficient sample size.\n\n')
            print('  #', category, 'merged into "Others" due to insufficient sample size.')



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Split MADC report to breeders and only retain haplotypes with more than 10 reads")

    parser.add_argument('passport',
                        help='A file containing the passport data of the samples in the report')

    parser.add_argument('report',
                        help='Reformated DArTag report, including fixed alleleIDs')

    parser.add_argument('-category_col_name',
                        help='Column containing the category info in the passport file', default='Breeder')

    parser.add_argument('-sample_number_threshold',
                        help='Minimum number of samples required in each category to be retained as a separate group. Categories with fewer samples will be merged into "Others". Default is 10.',
                        default=10)
    args=parser.parse_args()

    sample_lut, category_to_remove = collect_passport_data(args.passport, args.category_col_name, args.sample_number_threshold)

    split_samples_to_breeders(args.report, sample_lut, category_to_remove)
