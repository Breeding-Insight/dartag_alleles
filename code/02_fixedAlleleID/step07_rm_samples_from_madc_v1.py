#!/usr/bin/python3

def remove_samples_from_madc(rm_samples, report, outfile):
    import pandas as pd
    sample_list = []
    inp = open(rm_samples)
    line = inp.readline()
    while line:
        line_array = line.strip().split(',')
        sample_list.append(line_array[0].strip())
        line = inp.readline()
    inp.close()
    
    df = pd.read_csv(report)
    cnt = 0
    for sample in sample_list:
        if sample in df.columns:
            df = df.drop(columns=[sample])
            cnt += 1
        else:
            print(f"Sample {sample} not found in the report columns.")
    print(f"Removed {cnt} samples from the report.")
    df.to_csv(outfile, index=False)
    




if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")
    
    parser.add_argument('rm_samples',
                        help='A file containing sample names to be removed from the MADC report')

    parser.add_argument('report',
                        help='MADC report with allele name reformatted and unique sample names')

    parser.add_argument('output',
                        help='Output file name')

    args=parser.parse_args()

    remove_samples_from_madc(args.rm_samples, args.report, args.output)
