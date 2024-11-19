#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def add_info_to_pp(info, info_key, passport):
    # The join() function performs a left join by default, so each of the indexes in the first DataFrame are kept.
    import pandas as pd
    df_info = pd.read_csv(info, index_col='Sample_ID')
    # If sample IDs are pure number, when reading the data, they may be of type of int or str
    df_info.index = df_info.index.map(str)
    
    df_passport = pd.read_csv(passport, index_col='Sample_ID')
    df_passport = df_passport.dropna(axis=1, how='all')
    df_passport.index = df_passport.index.map(str)
    
    df_join = df_passport.join(df_info) # lef join
    df_join.fillna(101, inplace=True)
    print(df_join)
    outf = passport.replace('.csv', '_' + info_key + '.csv')
    df_join.to_csv(outf)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Add sample ID to DArTag report")

    parser.add_argument('info',
                        help='Info to be added to PCA')

    parser.add_argument('info_key', help='Info key to be added to the output file name')

    parser.add_argument('passport',
                        help='Comm-delimited DArTag reports with sample ID')

    args=parser.parse_args()

    add_info_to_pp(args.info, args.info_key, args.passport)
