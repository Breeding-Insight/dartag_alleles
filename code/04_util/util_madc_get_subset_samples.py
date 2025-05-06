#!/usr/bin/python3
# Generate a key file for sfetch, then run the following
# esl-sfetch -Cf [options] seqfile subseq-coord-file
# (retrieve multiple subsequences using file of keys and coords)
# /programs/hmmer-3.1b1-linux-intel-x86_64/binaries/esl-sfetch -Cf input.fasta sfetch_keyfile.txt


import pandas as pd

def get_sample_IDs(samples):
    """
    Get sample IDs from the samples file
    """
    df = pd.read_csv(samples)
    sample_IDs = ["CloneID", "AlleleSequence"] + df['Sample_ID'].tolist()
    print('# Number of samples to be extracted: ', len(sample_IDs)-2)
    return(sample_IDs)


def get_subset_samples(sample_IDs, report, outf):
    df = pd.read_csv(report, index_col='AlleleID')
    df_sub = df[sample_IDs]
    df_sub.to_csv(outf)
            
    

if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Get subset of samples from DArTag MADC report")

    parser.add_argument('samples',
                        help='Samples to be extracted from the report')

    parser.add_argument('report',
                        help='DArTag report')
    
    parser.add_argument('outf',
                        help='Output file name')

    args=parser.parse_args()

    sample_IDs = get_sample_IDs(args.samples)

    get_subset_samples(sample_IDs, args.report, args.outf)
