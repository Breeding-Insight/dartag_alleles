#!/usr/bin/python3

def concat_amcc(reports, outfile):
    import pandas as pd
    reports_list = reports.strip().split(',')
    df = pd.read_csv(reports_list[0], skiprows=7)
    df =df.set_index('AlleleID')
    df = df.drop(columns=['CloneID', 'AlleleSequence', 'AvgCountRef', 'AvgCountSnp'])
    print('Are there duplicate row names in', reports_list[0], '?', df[df.index.duplicated()])
    index = 1
    while index < len(reports_list):
        df1 = pd.read_csv(reports_list[index], skiprows=7)
        df1 =df1.set_index('AlleleID')
        df1 = df1.drop(columns=['CloneID', 'AlleleSequence', 'AvgCountRef', 'AvgCountSnp'])
        print('Are there duplicate row names in', reports_list[index], '?', df1[df1.index.duplicated()])
        df_concat = pd.concat([df, df1], axis=1)
        print(df_concat.index)
        df_concat.reset_index(inplace=True)
        df_concat = df_concat.rename(columns = {'index': 'AlleleID'})
        index += 1
    print(df_concat.columns)
    df_concat.to_csv(outfile, index=False)





if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('reports',
                        help='Comma-separated missing allele reports with allele name reformatted and unique sample names')

    parser.add_argument('output',
                        help='Output file name')

    args=parser.parse_args()

    concat_amcc(args.reports, args.output)