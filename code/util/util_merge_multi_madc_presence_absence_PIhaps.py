#!/usr/bin/python3

def concat_madc(reports, outfile):
    import pandas as pd
    reports_list = reports.strip().split(',')
    #cols = list(pd.read_csv(reports_list[0], nrows=1))
    df = pd.read_csv(reports_list[0], index_col='AlleleID')
    index = 1
    while index < len(reports_list):
        df1 = pd.read_csv(reports_list[index], index_col='AlleleID')
        df1 = df1.drop(columns=['AlleleSequence'])
        # The join() function performs a left join by default, so each of the indexes in the first DataFrame are kept.
        df = df.join(df1, how='outer')
        index += 1
        
    # Change row index into a column
    df = df.reset_index()
    df = df.rename(columns={'index': 'AlleleID'})
    df.sort_values(by=['AlleleID'], inplace=True)
    print(df.columns)
    df.to_csv(outfile, index=False)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('reports',
                        help='Comma-separated missing allele reports with allele name reformatted and unique sample names')

    parser.add_argument('output',
                        help='Output file name')

    args=parser.parse_args()

    concat_madc(args.reports, args.output)
