#!/usr/bin/python3

def plot_corr(corr):
    import matplotlib.pyplot as plt
    plt.plot(corr)
    plt.title('Correlation between ')
    plt.show()
    
    
def cal_rep_corr(dosage):
    # updog dosage calls with missing data denoted as NA
    import pandas as pd
    df = pd.read_csv(dosage)
    # Set 1st column as row index and give it a new name 'AlleleID'
    df.set_index('Unnamed: 0', inplace=True)
    idx = df.index
    idx.rename('AlleleID', inplace=True)
    samples = df.columns
    # All replicates end with '_1', '_2', '_3'
    uniq_accessions = {}
    for i in samples:
        if i.count('_') == 1:
            i_array = i.split('_')
            if i_array[0] not in uniq_accessions:
                uniq_accessions[i_array[0]] = [i]
            else:
                uniq_accessions[i_array[0]].append(i)
        else:
            print('This name has more than one underscore', i)
    outp = open(dosage.replace('.csv', '_corr.csv'), 'w')
    outp.write('Rep1,Rep2,Corr\n')
    plot = 0
    for j, h in uniq_accessions.items():
        plot += 1
        # Write correlation matrix of 4 accessions to files for plotting.
        if plot < 5:
            df_accession = df[h]
            corr_4plot = df_accession.corr()
            outf_corr = dosage.replace('.csv', '_corr_'+j+'.csv')
            corr_4plot.to_csv(outf_corr)
            
        # Output pairwise correlations
        if len(h) == 3:
            corr = df[h[0]].corr(df[h[1]])
            outp.write(','.join([h[0], h[1], str(corr)]) + '\n')
            corr = df[h[0]].corr(df[h[2]])
            outp.write(','.join([h[0], h[2], str(corr)]) + '\n')
            corr = df[h[1]].corr(df[h[2]])
            outp.write(','.join([h[1], h[2], str(corr)]) + '\n')
        elif len(h) == 2:
            corr = df[h[0]].corr(df[h[1]])
            outp.write(','.join([h[0], h[1], str(corr)]) + '\n')
        else:
            print('No replicates ', j, h)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="This script is generated based on Alfalfa Project 9")
    
    parser.add_argument('dosage',
                        help='updog dosage with missing data denoted as NA')

    args=parser.parse_args()

    cal_rep_corr(args.dosage)
