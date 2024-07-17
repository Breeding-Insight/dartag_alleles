#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def join_2_dict(file1, file2, outf):
    # The join() function performs a left join by default, so each of the indexes in the first DataFrame are kept.
    import pandas as pd
    df1 = pd.read_csv(file1, sep='\t')
    df1.columns = ['chromosome', 'position', 'ID', 'Ref', 'Alt']
    # Add leading zeros to positions to the total length of 9 digit
    df1['CloneID'] = df1['chromosome'] + '_' + df1['position'].apply(lambda x: '{0:0>9}'.format(x))
    df1 = df1.set_index('CloneID')

    df2 = pd.read_csv(file2, sep='\t')
    df2.columns = ['polyRAD_chromosome', 'polyRAD_position', 'polyRAD_ID', 'polyRAD_Ref', 'polyRAD_Alt']
    # Add leading zeros to positions to the total length of 9 digit
    df2['CloneID'] = df2['polyRAD_chromosome'] + '_' + df2['polyRAD_position'].apply(lambda x: '{0:0>9}'.format(x))
    df2 = df2.set_index('CloneID')
    df_join = df1.join(df2, how='outer', rsuffix='_02')
    df_join.reset_index(inplace=True)
    df_join = df_join.rename(columns={'index': 'CloneID'})
    df_join.to_csv(outf, index=False)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Add sample ID to DArTag report")

    parser.add_argument('file1',
                        help='file1')

    parser.add_argument('file2',
                        help='file2')

    parser.add_argument('outf', help='output file name')
    
    args=parser.parse_args()

    join_2_dict(args.file1, args.file2, args.outf)
