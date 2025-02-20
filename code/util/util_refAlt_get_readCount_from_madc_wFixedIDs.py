#!/usr/bin/python3


def ext_RefAlt_readCount(report):
    inp = open(report)
    line = inp.readline()
    ref_alt_dict = {}
    while line:
        line_array = line.strip().split(',')
        # AlleleID,CloneID,AlleleSequence,sample1,
        if 'AlleleID' == line_array[0]:
            columns = line_array[1:]
        else:
            if line_array[0].endswith('|Ref_0001') or line_array[0].endswith('|Alt_0002'):
                ref_alt_dict[line_array[0]] = line_array[1:]
            else:
                pass
        line = inp.readline()
    inp.close()
    
    # Remove negative controls
    import pandas as pd
    df = pd.DataFrame.from_dict(ref_alt_dict, orient='index', columns=columns)
    print('## All samples:\n', df)
    for i in df.columns.tolist():
        if 'CONTROL' in i or 'control' in i:
            df = df.drop(i, axis=1)
    print('## After removing negative controls:\n', df)

    # Loop through each marker locus, and make the alleles in order of "Ref", "Alt", "RefMatch", "AltMatch"
    # !!!!!! This re-ordering will take a few minutes !!!!
    df_groupby = df.groupby('CloneID')
    index_ordered = []
    for cloneID, clone_df in df_groupby:
        ref = cloneID + '|Ref_0001'
        alt = cloneID + '|Alt_0002'
        index_ordered.append(ref)
        index_ordered.append(alt)
    print('# Order index: ', index_ordered[:10])
    df_ordered = df.reindex(index_ordered)
    outp_report = report.replace('.csv', '_RefAlt.csv')
    df_ordered.to_csv(outp_report, index_label='AlleleID')



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract read counts from Reference and Alternative alleles only")
    
    parser.add_argument('report', help='Raw MADC')

    args=parser.parse_args()
    
    ext_RefAlt_readCount(args.report)
