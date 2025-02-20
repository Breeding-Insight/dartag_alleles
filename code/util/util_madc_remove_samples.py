#!/usr/bin/python3
# Parse the original file that DArT provided



def remove_samples_from_madc(madc, samples):
    import pandas as pd
    df = pd.read_csv(madc, sep=',')
    remove = samples.strip().split(',')
    df_sub = df.drop(columns=remove)
    outp = madc.replace('.csv', '_rmSamples.csv')
    df_sub.to_csv(outp, index=False)
    print('Number of samples to be removed: ', len(remove))
    print('Input file: ', df.shape)
    print('After removing samples: ', df_sub.shape)
    


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Remove samples from MADC report")

    parser.add_argument('madc',
                        help='DArTag MADC after running filter missing data (AlleleID	CloneID	AlleleSequence	1_33642_W6_P1)')

    parser.add_argument('samples',
                        help='A comma-delimited sample names to be removed from MADC report')
                        
    args=parser.parse_args()

    remove_samples_from_madc(args.madc, args.samples)
