#!/usr/bin/python3

import pandas as pd
import numpy as np


def calculate_marker_stats(species, species_dose_df, total_alleles, species_samples, species_ploidy, outf, stats, total_markers, results_all, ind_maf_output):
    '''                
                            chr01_000449004 chr01_000449027 ... M6_chr11_40262586_000000129
    S_WDP041__Halterman1_G3                   1               0               .
    S_WDP041__Halterman1_G4                   1               .               0
    '''
    # Vectorize operations instead of using loops for better performance
    # Convert to numeric once
    numeric_species_dose_df = species_dose_df.replace('.', np.nan).astype(float)
    
    # Vectorized calculation of miss_count for all markers at once
    miss_counts = np.sum(np.isnan(numeric_species_dose_df), axis=1)
    ref_alleles = np.nansum(numeric_species_dose_df, axis=1)
    effective_alleles = total_alleles - miss_counts * int(species_ploidy)
    effective_samples = len(species_samples) - miss_counts
    
    # Vectorized frequency calculation
    with np.errstate(divide='ignore', invalid='ignore'):
        ref_freqs = np.where(effective_alleles != 0, ref_alleles/effective_alleles, np.nan)
    alt_freqs = 1 - ref_freqs
    mafs = np.minimum(ref_freqs, alt_freqs)
    
    # Vectorized heterozygosity calculations
    het_counts = np.sum((numeric_species_dose_df > 0) & (numeric_species_dose_df < int(species_ploidy)), axis=1)
    het_freqs = het_counts / effective_samples
    
    # Create results DataFrame directly
    result_df = pd.DataFrame(
        {'Species': species,
         'Marker_ID': species_dose_df.index,
         'Total_samples': len(species_samples),
         'Miss_samples': miss_counts,
         'Total_effective_samples': effective_samples,
         'Total_alleles': total_alleles,
         'Total_effective_alleles': effective_alleles,
         'Ref_allele': ref_alleles,
         'Ref_freq': np.round(ref_freqs, 2),
         'AAF': np.round(alt_freqs, 2),
         'MAF': np.round(mafs, 2),
         'Het_count': het_counts,
         'obHet': np.round(het_freqs, 2)
         })

    # Append results to the all results DataFrame
    results_all = pd.concat([results_all, result_df], ignore_index=True, axis=0)
    
    if ind_maf_output == 'Y':
    # Save results to a CSV file    
        result_df.to_csv(outf, index=False)
    else:
        pass

    # Summary statistics
    # get the number of markers with heterozygosity
    # has to use "&" instead of "and" because we are working with pandas Series
    het_count = (result_df['obHet'] > 0).sum()
    maf_count = ((result_df['MAF'] > 0.05) & (result_df['MAF'] <= 0.5)).sum()
    het_percentage = round((het_count / total_markers) * 100, 2)
    maf_percentage = round((maf_count / total_markers) * 100, 2)
    if 'Species' not in stats:
        stats['Species'] = [species]
        stats['Samples'] = [len(species_samples)]
        stats['Total_markers'] = [total_markers]
        stats['maf_count_005i05'] = [maf_count]
        stats['maf_percentage_markers'] = [maf_percentage]
        stats['het_count'] = [het_count]
        stats['het_percentage_markers'] = [het_percentage]        
    else:
        stats['Species'].append(species)
        stats['Samples'].append(len(species_samples))
        stats['Total_markers'].append(total_markers)
        stats['maf_count_005i05'].append(maf_count)
        stats['maf_percentage_markers'].append(maf_percentage)
        stats['het_count'].append(het_count)
        stats['het_percentage_markers'].append(het_percentage)
        
    # clean up large objects when no longer needed
    del numeric_species_dose_df, miss_counts, ref_alleles, effective_alleles, effective_samples, ref_freqs, alt_freqs, mafs, het_counts, het_freqs
    import gc
    gc.collect()
    return(results_all, stats)
    
    
def cal_allele_status_by_species(dosage, passport, ploidy, ind_maf):
    '''
    # updog dosages represent doses of reference alleles
    Calculate allele status by species from dosage
    - Alternative allele frequency
    - Minor allele frequency
    - Observed heterozygosity (Ho)
    '''
    stats = {}
    # Read the dosage file
    dose_df = pd.read_csv(dosage, index_col=0, low_memory=False)
    dose_df.index = dose_df.index.str.replace('"', '', regex=True)
    total_markers = len(dose_df.index)
    
    passport_df = pd.read_csv(passport)
    
    # Create species mapping from passport file
    species_mapping = dict(zip(passport_df['Sample_ID'], passport_df['Species']))
    
    # Add species row using column names
    species_row = pd.Series(dose_df.columns.map(species_mapping), index=dose_df.columns, name='Species')
    dose_df = pd.concat([pd.DataFrame([species_row]), dose_df])

    # Initialize results DataFrame
    results_all = pd.DataFrame(columns=['Species', 'Marker_ID', 'Total_samples', 'Miss_samples', 'Total_effective_samples', 'Total_alleles', 'Total_effective_alleles', 'Ref_allele', 'Ref_freq', 'AAF', 'MAF', 'Het_count', 'obHet'])
    for species in species_row.unique():
        if pd.isna(species):
            continue
            
        print(f'\n\n# Processing species: {species}')
        
        # Get samples for current species
        species_samples = species_row[species_row == species].index.tolist()
        total_alleles = len(species_samples) * int(ploidy)
        
        # Filter the dosage DataFrame to only include the samples for the current species
        species_dose_df = dose_df[species_samples]

        # Remove the first row (species row) from the dosage DataFrame
        species_dose_df = species_dose_df.iloc[1:]
        print(species_dose_df)
        
        outf = dosage.replace('.csv', f'_{species}_maf.csv')
        
        # calculate allele status
        results_all, stats = calculate_marker_stats(species, species_dose_df, total_alleles, species_samples, ploidy, outf, stats, total_markers, results_all, ind_maf)

    # Save results to a CSV file
    results_all.to_csv(dosage.replace('.csv', '_maf_all.csv'), index=False)
    
    # Save the summary to a CSV file
    stats_df = pd.DataFrame(stats)
    summary = dosage.replace('.csv', '_maf_het_summary.csv')
    stats_df.to_csv(summary, index=False)
    outp = open(summary, 'a')
    now = pd.Timestamp.now()
    outp.write('\n\nProcessed on ' + str(now) + '\n')
    outp.write('*het_count: Number of markers with heterozygosity >0\n')
    outp.write('*het_percentage_markers: Percentage of markers with heterozygosity >0\n')
    outp.write('*maf_count_005i05: Number of markers with MAF between 0.05 and 0.5\n') 
    outp.write('*maf_percentage_markers: Percentage of markers with MAF between 0.05 and 0.5\n')



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")
    
    parser.add_argument('dosage',
                        help='updog dosage with missing data denoted as NA')
    
    parser.add_argument('passport',
                        help='Passport file with species information')
    
    parser.add_argument('ploidy',
                        help='Ploidy of the species')
    
    parser.add_argument('ind_maf_output',
                        help='Y to output individual maf and het statistics, N to not output')

    args=parser.parse_args()

    cal_allele_status_by_species(args.dosage, args.passport, args.ploidy, args.ind_maf_output)
