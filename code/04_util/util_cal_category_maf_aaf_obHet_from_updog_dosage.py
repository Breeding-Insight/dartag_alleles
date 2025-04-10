#!/usr/bin/python3

import pandas as pd
import numpy as np

def calculate_marker_stats(species, dose_df_species, total_alleles, species_samples, species_ploidy, outf, stats, total_markers, results_all):
    results = []
    for markers in dose_df_species.columns:
        ref_allele = dose_df_species[markers].sum()
        
        # calculate frequencies
        ref_freq = ref_allele / total_alleles
        alt_freq = 1 - ref_freq
        
        # MAF is the lower frequency of the two alleles
        maf = min(ref_freq, alt_freq)
        
        # Heterozygosity is the probability that two alleles are different
        het_count = sum((dose_df_species[markers] > 0) & (dose_df_species[markers] < int(species_ploidy)))
        het_freq = float(het_count) / len(species_samples)

        results.append({
            'Species': species,
            'Marker_ID': markers,
            'Total_alleles': total_alleles,
            'Ref_allele': ref_allele,
            'Ref_freq': round(ref_freq,2),
            'AAF': round(alt_freq, 2),
            'MAF': round(maf, 2),
            'Samples': len(species_samples),
            'Het_count': het_count,
            'obHet': round(het_freq, 2)
        })

        results_all.append(
            {'Species'    : species, 'Marker_ID': markers, 'Total_alleles': total_alleles, 'Ref_allele': ref_allele,
                'Ref_freq': round(ref_freq, 2), 'AAF': round(alt_freq, 2), 'MAF': round(maf, 2),
                'Samples' : len(species_samples), 'Het_count': het_count, 'obHet': round(het_freq, 2)})
    
    # Save results to a CSV file
    result_df = pd.DataFrame(results)
    result_df.to_csv(outf, index=False)
    
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
        stats['het_count'] = [het_count]
        stats['maf_count_005i05'] = [maf_count]
        stats['het_percentage_markers'] = [het_percentage]
        stats['maf_percentage_markers'] = [maf_percentage]
    else:
        stats['Species'].append(species)
        stats['Samples'].append(len(species_samples))
        stats['Total_markers'].append(total_markers)
        stats['het_count'].append(het_count)
        stats['maf_count_005i05'].append(maf_count)
        stats['het_percentage_markers'].append(het_percentage)
        stats['maf_percentage_markers'].append(maf_percentage)
    return(results_all, stats)


# With additional details
def detailed_count(df):
    total_markers = len(df.columns)
    between_counts = ((df > 0) & (df < 1)).sum()
    percentage = (between_counts / total_markers) * 100
    
    return pd.DataFrame({'Total_Markers': total_markers, 'Between_0_1': between_counts, 'Percentage': percentage})

    
def cal_allele_status_by_species(dosage, passport, ploidy):
    '''
    # updog dosages represent doses of reference alleles
    Calculate allele status by species from dosage
    - Alternative allele frequency
    - Minor allele frequency
    - Observed heterozygosity (Ho)

    :param dosage: 
    :param passport: 
    :return: 
    '''
    stats = {}
    passport_df = pd.read_csv(passport, index_col='Sample_ID')
    passport_columns = passport_df.columns.tolist()
    dose_df = pd.read_csv(dosage, index_col=0)
    dose_df = dose_df.transpose()
    total_markers = len(dose_df.columns)
    
    # Concatenate the two dataframes along the columns
    # The join() function performs a left join by default, so each of the indexes in the first DataFrame are kept.
    df_join = dose_df.join(passport_df)
    results_all = []
    
    for species in df_join['Species'].unique():
        # Get the sample IDs for the current species
        species_samples = df_join[df_join['Species'] == species].index.tolist()
        total_alleles = len(species_samples) * int(ploidy)
        
        # Filter the dosage DataFrame to only include the samples for the current species
        dose_df_species = df_join.loc[species_samples]
        # Drop the passport columns from the dosage DataFrame
        dose_df_species = dose_df_species.drop(columns=passport_columns)

        # calculate allele status
        outf = dosage.replace('.csv', '_' + species.replace(' ', '_') + '_maf.csv')
        results_all, stats = calculate_marker_stats(species, dose_df_species, total_alleles, species_samples, ploidy, outf, stats, total_markers, results_all)

    # Save results to a CSV file
    result_df = pd.DataFrame(results_all)
    result_df.to_csv(dosage.replace('.csv', '_maf_all.csv'), index=False)
    
    # Save the summary to a CSV file
    stats_df = pd.DataFrame(stats)
    summary = dosage.replace('.csv', '_maf_het_summary.csv')
    stats_df.to_csv(summary, index=False)
    outp = open(summary, 'a')
    now = pd.Timestamp.now()
    outp.write('\n\nProcessed on ' + str(now) + '\n')
    outp.write('*het_count: Number of markers with heterozygosity >0\n')
    outp.write('*maf_count_005i045: Number of markers with MAF between 0.05 and 0.5\n')
    outp.write('*het_percentage_markers: Percentage of markers with heterozygosity >0\n')
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

    args=parser.parse_args()

    cal_allele_status_by_species(args.dosage, args.passport, args.ploidy)
