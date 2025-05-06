#!/usr/bin/python3

import pandas as pd
import numpy as np
from scipy import stats


def calculate_marker_informativeness(dosage_file, passport_file):
	"""
	Analyze marker informativeness between parent species
	"""
	# Read data
	dose_df = pd.read_csv(dosage_file, index_col=0, low_memory=False)
	passport_df = pd.read_csv(passport_file)
	
	# Get parent species samples
	micro_samples = passport_df[passport_df['Species'] == 'Solanum microdontum']['Sample_ID'].tolist()
	tuber_samples = passport_df[passport_df['Species'] == 'Solanum tuberosum']['Sample_ID'].tolist()
	
	results = []
	
	for marker in dose_df.index:
		if marker == 'Species':
			continue
		
		# Get allele frequencies for each species
		micro_freq = np.mean(pd.to_numeric(dose_df.loc[marker, micro_samples], errors='coerce'))
		tuber_freq = np.mean(pd.to_numeric(dose_df.loc[marker, tuber_samples], errors='coerce'))
		
		# Calculate frequency difference
		freq_diff = abs(micro_freq - tuber_freq)
		
		# Simple FST calculation
		# This is a simplified version, you might want to use a more sophisticated FST calculation
		total_freq = np.mean([micro_freq, tuber_freq])
		if total_freq > 0 and total_freq < 1:
			fst = freq_diff / (total_freq * (1 - total_freq))
		else:
			fst = np.nan
		
		results.append(
			{'Marker': marker, 'Micro_freq': micro_freq, 'Tuber_freq': tuber_freq, 'Freq_diff': freq_diff, 'FST': fst})
	
	# Convert to DataFrame
	results_df = pd.DataFrame(results)
	
	# Add marker categories
	results_df['Diagnostic'] = results_df['FST'] > 0.7  # You can adjust this threshold
	
	return results_df


def save_results(results, dose):
	"""
	Save results to CSV
	"""

	outf = dose.replace('.csv', '_marker_informativeness.csv')
	# Save results
	results.to_csv(outf, index=False)
	
	# Summary statistics
	print("\nSummary of marker informativeness:")
	print(f"Total markers analyzed: {len(results)}")
	print(f"Highly differentiated markers (FST > 0.7): {sum(results['Diagnostic'])}")
	print("\nFST distribution:")
	print(results['FST'].describe())
	
	# You might want to plot the results
	import matplotlib.pyplot as plt
	
	plt.figure(figsize=(10, 6))
	plt.hist(results['FST'].dropna(), bins=50)
	plt.xlabel('FST')
	plt.ylabel('Count')
	plt.title('Distribution of FST values between parent species')
	plt.savefig('fst_distribution.png')
	plt.close()
	
	# Plot frequency comparison
	plt.figure(figsize=(10, 10))
	plt.scatter(results['Micro_freq'], results['Tuber_freq'], alpha=0.5)
	plt.xlabel('S. microdontum frequency')
	plt.ylabel('S. tuberosum frequency')
	plt.title('Allele frequency comparison between parent species')
	plt.plot([0, 1], [0, 1], 'r--')  # diagonal line
	plt.savefig('frequency_comparison.png')
	plt.close()


if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="")
	
	parser.add_argument('dose',
						help='Updog dosage file with missing data marked as NA')
	
	parser.add_argument('passport',
						help='Passport file with sample IDs and species information')
	
	args = parser.parse_args()
	
	# Usage:
	results = calculate_marker_informativeness(args.dose, args.passport)
	
	save_results(results, args.dose)
