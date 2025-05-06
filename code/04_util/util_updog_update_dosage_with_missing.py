#!/usr/bin/python3

import pandas as pd
import numpy as np

def get_totalDepth_perInd_perMarker(read_counts_vcf):
    inp = open(read_counts_vcf)
    line = inp.readline()
    totalDepth_lut = {}
    while line:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                line_array = line.strip().split('\t')
                sample_names = line_array[9:]
            else:
                pass
        else:
            line_array = line.strip().split('\t')
            markerID = line_array[0] + '_' + line_array[1].zfill(9)
            totalDepth_lut[markerID] = []
            index = 9
            while index < len(line_array):
                dp = int(line_array[index].split(':')[0])
                # DP:RA:AD	137:83:83,54
                totalDepth_lut[markerID].append(dp)
                index += 1
        line = inp.readline()
    inp.close()
    return(sample_names, totalDepth_lut)


def update_dosage(dose, sample_names, totalDepth_lut, threshold):
    df_dose = pd.read_csv(dose, index_col=0)
    df_dose = df_dose[sample_names]
    
    df_counts = pd.DataFrame.from_dict(totalDepth_lut, orient='index', columns=sample_names)
    # Make df_counts the same shape as df_dose
    # keep only columns and rows that exist in df_dose
    df_counts = df_counts[df_dose.columns]
    df_counts = df_counts[df_counts.index.isin(df_dose.index)]
      
    """
        Update df_dose with NaN where df_counts values are less than threshold
        Includes error checking and column validation
    """
    try:
        # Check if DataFrames have same shape
        if df_dose.shape != df_counts.shape:
            raise ValueError("DataFrames must have same shape")
    
        # Check if columns match
        if not all(df_dose.columns == df_counts.columns):
            raise ValueError("DataFrame columns must match")
    
        # Update values
        mask = df_counts < int(threshold)
        df_dose[mask] = "."
        
        # Count the number of data points
        total_points = df_dose.size
        
        # Count the number of missing values
        miss_count = (df_dose == '.').sum().sum()
        miss_fraction = round(miss_count / total_points * 100, 2)
        print(f"# Total number of datapoints: {total_points}")
        print(f"# Number of datapoints set to missing: {miss_count}")
        print(f"# Fraction of missing datapoints: {miss_fraction}%")
        
        # Add index names to df_dose
        df_dose.index.name = 'AlleleID'
        outp_report = dose.replace('.csv', '_updateNA.csv')
        df_dose.to_csv(outp_report)
        return df_dose

    except Exception as e:
        print(f"Error: {e}")
        return None



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('read_counts_vcf',
                        help='VCF files with read counts, generated from hap2snpvcf_dose_2_vcf.py')
    
    parser.add_argument('dose',
                        help='updog dosage in csv format')
    
    parser.add_argument('threshold',
                        help='Threshold to determine missing data')

    args=parser.parse_args()
    
    sample_names, totalDepth_lut = get_totalDepth_perInd_perMarker(args.read_counts_vcf)
    
    update_dosage(args.dose, sample_names, totalDepth_lut, args.threshold)
