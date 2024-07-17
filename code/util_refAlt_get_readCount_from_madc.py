#!/usr/bin/python3


def ext_RefAlt_readCount(report):
    inp = open(report)
    line = inp.readline()
    data_col = 0
    ref_alt_dict = {}
    while line:
        line_array = line.strip().split(',')
        if line_array[0] == '*':
            data_col = line_array.count('*')
        elif 'AlleleID' == line_array[0]:
            columns = line_array[data_col:]
        else:
            if line_array[0].endswith('|Ref') or line_array[0].endswith('|Alt'):
                ref_alt_dict[line_array[0]] = line_array[data_col:]
            else:
                pass
        line = inp.readline()
    inp.close()
    
    # Remove negative controls
    import pandas as pd
    df = pd.DataFrame.from_dict(ref_alt_dict, orient='index', columns=columns)
    print('#All samples:\n', df)
    for i in df.columns.tolist():
        if 'CONTROL' in i or 'control' in i:
            df = df.drop(i, axis=1)
    print('After removing negative controls:\n', df)
    outp_report = report.replace('.csv', '_RefAlt.csv')
    df.to_csv(outp_report, index_label='AlleleID')



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract read counts from Reference and Alternative alleles only")
    
    parser.add_argument('report', help='Raw MADC')

    args=parser.parse_args()
    
    ext_RefAlt_readCount(args.report)
