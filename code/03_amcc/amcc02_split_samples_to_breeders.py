#!/usr/bin/python3

def collect_passport_data(passport):
    inp = open(passport)
    line = inp.readline() # header
    # ['Source', 'Breeder', 'Study_ID', 'Plate_ID', 'Well_ID', 'SampleName', 'Plant_ID', 'IVP', 'IVNO', 'Taxon', 'ng/ul ', 'G12/H12 Blanks', 'SampleName length QC']
    line_array = line.strip().split('\t')
    print('Header:', line_array)
    line = inp.readline()
    sample_lut = {}
    sample_count = 0
    while line:
        line_array = line.strip().split('\t')
        sample_count += 1
        if line_array[1] not in sample_lut:
            sample_lut[line_array[1]] = ['AlleleID',line_array[5].strip()]
        else:
            sample_lut[line_array[1]].append(line_array[5].strip())
        line = inp.readline()
    inp.close()
    print('Number of samples in the passport data: ', sample_count)
    print(sample_lut.keys())
    return(sample_lut)



def split_samples_to_breeders(report, sample_lut):
    import pandas as pd
    df = pd.read_csv(report)
    cols = df.columns.tolist()
    print(cols)
    samples_not_in_report = []
    print('Number of samples in the report: ', len(cols) - 1)
    for breeder, sample_list in sample_lut.items():
        print('original sample list:', len(sample_list))
        for sample in sample_list:
            if sample not in cols:
                sample_list.remove(sample)
                samples_not_in_report.append(sample)
            else:
                pass

        print(breeder, ':', len(samples_not_in_report), 'samples not present in report: ', samples_not_in_report)
        print('After cleaning:', len(sample_list))
        df_subset = df[sample_list]
        outf = report.replace('.csv', '_' + breeder + '.csv')
        df_subset.to_csv(outf, index=False)
        samples_not_in_report = []


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