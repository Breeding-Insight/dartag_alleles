#!/usr/bin/python3

def concat_madc(reports, outfile):
    import pandas as pd
    inp = open(reports)
    line = inp.readline()
    #cols = list(pd.read_csv(reports_list[0], nrows=1))
    df = pd.read_csv(line.strip(), index_col='AlleleID')
    line = inp.readline()
    while line:
        df1 = pd.read_csv(line.strip(), index_col='AlleleID')
        # The join() function performs a left join by default, so each of the indexes in the first DataFrame are kept.
        df = df.join(df1, how='outer')
        line = inp.readline()
    inp.close()

    # Sort the alleles from Ref, Alt, RefMatch, AltMatch
    df.fillna(0, inplace=True)
    df['CloneID'] = df.index.str.split('|', n=1).str[0]
    # Move to the first column
    first_col = df.pop('CloneID')
    df.insert(0, 'CloneID', first_col)

    df_groupby = df.groupby('CloneID')
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
    df_ordered.fillna(0, inplace=True)
    df_ordered.to_csv(outfile, index=True)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('reports',
                        help='A file contains MADC reports with allele name reformatted and unique sample names, one file per line')

    parser.add_argument('output',
                        help='Output file name')

    args=parser.parse_args()

    concat_madc(args.reports, args.output)
