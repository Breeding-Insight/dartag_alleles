#!/usr/bin/python3
# Generate a key file for sfetch, then run the following
# esl-sfetch -Cf [options] seqfile subseq-coord-file
# (retrieve multiple subsequences using file of keys and coords)
# /programs/hmmer-3.1b1-linux-intel-x86_64/binaries/esl-sfetch -Cf input.fasta sfetch_keyfile.txt


def parse_dartag(plateIDs, report):
    import pandas as pd
    plateID_list = plateIDs.split(',')
    df = pd.read_csv(report, header=[5,7])
    df = df.T
    df_gb_plates = df.groupby(level=0)
    df_subset = pd.DataFrame()
    for index, row in df_gb_plates:
        if index not in plateID_list:
            df_subset = pd.concat([df_subset, row], axis=0)
        else:
            pass
    df_subset_t = df_subset.T
    outp = report.replace('.csv', '_subset.csv')
    df_subset_t.to_csv(outp, index=False)
            
    

if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Update allele sequences in report with alleles assigned temporary names")

    parser.add_argument('plateIDs',
                        help='plate IDs separated by comma')

    parser.add_argument('report',
                        help='DArTag report')

    args=parser.parse_args()

    parse_dartag(args.plateIDs, args.report)