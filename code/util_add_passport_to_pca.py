#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def add_pp_info_to_pca(pca, passport):
    # The join() function performs a left join by default, so each of the indexes in the first DataFrame are kept.
    import pandas as pd
    df_pca = pd.read_csv(pca, index_col='Sample_ID')
    # If sample IDs are pure number, when reading the data, they may be of type of int or str
    df_pca.index = df_pca.index.map(str)
    
    df_passport = pd.read_csv(passport, index_col='Sample_ID')
    df_passport = df_passport.dropna(axis=1, how='all')
    df_passport.index = df_passport.index.map(str)
    
    df_join = df_pca.join(df_passport) # lef join
    print(df_join)
    outf = pca.replace('.csv', '_pp.csv')
    df_join.to_csv(outf)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Add sample ID to DArTag report")

    parser.add_argument('pca',
                        help='pca')

    parser.add_argument('passport',
                        help='Comm-delimited DArTag reports with sample ID')

    args=parser.parse_args()

    add_pp_info_to_pca(args.pca, args.passport)
