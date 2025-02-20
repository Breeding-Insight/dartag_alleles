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
        print("Number of sample with duplicate names: ", len(dup_names))
        print("Check the output file for those duplicate sample names.")
    else:
        print("No sample with duplicate names")
    return(sample_names)


def get_madc_for_samples(file, sample_names, report):
    df = pd.read_csv(report)
    df = df.set_index('AlleleID')
    columns = ['CloneID', 'AlleleSequence']
    noData = []
    for sample in sample_names:
        if sample in df.columns:
            columns.append(sample)
        else:
            noData.append(sample)
            
    if len(noData) != 0:
        outp = open(file.replace('.csv', '_noData.csv'), 'w')
        print("Number of samples not present in DArTag reports: ", len(noData))
        print(noData)
        for i in noData:
            outp.write(i + '\n')
        outp.close()
    else:
        print("Number of samples NOT present in DArTag reports: 0")
        
    df_sub = df[columns]
    df_sub.fillna(0, inplace=True)
    df_sub.sort_index(axis=0, inplace=True)
    print('Before removing haplotypes with no reads: ', len(df_sub.index.tolist()))
    for index, sum in df_sub.sum(axis=1, numeric_only=True).items():
        if sum == 0.0 and not index.endswith(('Ref_0001', 'Alt_0002')):
            df_sub.drop(index, inplace=True)

    # Loop through each marker locus, and make the alleles in order of "Ref", "Alt", "RefMatch", "AltMatch"
    # !!!!!! This re-ordering will take a few minutes !!!!
    df_groupby = df_sub.groupby('CloneID')
    df_ordered = pd.DataFrame()
    for cloneID, clone_df in df_groupby:
        ref = cloneID + '|Ref_0001'
        alt = cloneID + '|Alt_0002'
        reindex = [ref, alt]
        idx_sorted = sorted(clone_df.index.to_list())
        for i in idx_sorted:
            if 'RefMatch' in i:
                reindex.append(i)
            else:
                pass

        for i in idx_sorted:
            if 'AltMatch' in i:
                reindex.append(i)
            else:
                pass
        clone_df = clone_df.reindex(reindex)
        df_ordered = pd.concat([df_ordered, clone_df], axis=0)
    print('After removing haplotypes with no reads: ', len(df_ordered.index.tolist()))
    df_ordered.reset_index(inplace=True)
    df_ordered = df_ordered.rename(columns={'index': 'AlleleID'})
    
    outf = file.replace('.csv', '_madc.csv')
    df_ordered.to_csv(outf, index=False)




if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('samples',
                        help='A file containing sample names')

    parser.add_argument('report',
                        help='Reformatted DArTag report')
    
    
    args=parser.parse_args()

    sample_names = get_sample_names(args.samples)
    
    get_madc_for_samples(args.samples, sample_names, args.report)
