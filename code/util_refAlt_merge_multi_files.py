#!/usr/bin/python3

def concat_madc(reports, outfile):
    import pandas as pd
    reports_list = reports.strip().split(',')
    index = 0
    df_concat = pd.DataFrame()
    while index < len(reports_list):
        df1 = pd.read_csv(reports_list[index], index_col='AlleleID')
        df_concat = pd.concat([df_concat, df1], axis=1)
        index += 1
    print(df_concat.columns)
    df_concat.to_csv(outfile)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('reports',
                        help='Comma-separated missing allele reports with allele name reformatted and unique sample names')

    parser.add_argument('output',
                        help='Output file name')

    args=parser.parse_args()

    concat_madc(args.reports, args.output)
