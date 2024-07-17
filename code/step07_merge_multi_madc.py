#!/usr/bin/python3

def concat_madc(reports, outfile):
    import pandas as pd
    reports_list = reports.strip().split(',')
    index = 0
    df_concat = pd.DataFrame()
    while index < len(reports_list):
        df1 = pd.read_csv(reports_list[index], index_col='AlleleID')
        print('\n### Are there duplicate row names in', reports_list[index], '?', df1[df1.index.duplicated()])
        if 'readCountSum' in df1.columns:
            df1 = df1.drop(columns=['readCountSum'])
        else:
            pass
        
        if 'CloneID' in df_concat.columns:
            df1 = df1.rename(columns={'CloneID': 'CloneID_01', 'AlleleSequence': 'AlleleSequence_01'})
            df_concat = pd.concat([df_concat, df1], axis=1)
            # Replace missing values for CloneID and AlleleSequence and drop CloneID_01 and AlleleSEquence_01
            df_concat['CloneID'] = df_concat['CloneID'].fillna(df_concat.pop('CloneID_01'))
            df_concat['AlleleSequence'] = df_concat['AlleleSequence'].fillna(df_concat.pop('AlleleSequence_01'))
            df_concat.fillna(0, inplace=True)

            # Drop 'CloneID_01', 'AlleleSequence_01'
            #df_concat = df_concat.drop(columns=['CloneID_01', 'AlleleSequence_01'])
            df_concat.reset_index(inplace=True)
            df_concat = df_concat.rename(columns={'index': 'AlleleID'})
        else:
            df_concat = df1
        index += 1
 
    # Sort the alleles from Ref, Alt, RefMatch, AltMatch
    df_concat = df_concat.set_index('AlleleID')
    df_groupby = df_concat.groupby('CloneID')
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
    df_ordered.to_csv(outfile, index=True)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('reports',
                        help='Comma-separated missing allele reports with allele name reformatted and unique sample names')

    parser.add_argument('output',
                        help='Output file name')

    args=parser.parse_args()

    concat_madc(args.reports, args.output)
