#!/usr/bin/python3

import pandas as pd

def calculate_coverage(madc):
	# Load the CSV file into a pandas DataFrame
	df = pd.read_csv(madc)

	# Ensure the columns where the samples begin
	# Let's assume the sample columns start at column index 3 onward (0-based index)
	sample_columns = df.columns[3:]

	# Calculate marker-level coverage
	# For each CloneID, the coverage is the sum of all read depths across all alleles (rows with the same CloneID) and all samples. Essentially, it tells how well this marker was sequenced across all alleles and samples collectively.
	# For each sample, the coverage is the sum of all the read depths across all alleles for all markers. This provides information on how much sequencing depth or coverage a specific sample received.
	# Group rows by `CloneID` and calculate the sum of read depths for each marker
	marker_coverage = df.groupby('CloneID')[sample_columns].sum()
	marker_coverage['Total_coverage'] = marker_coverage.sum(axis=1)
	marker_coverage['Mean_sample_coverage'] = marker_coverage[sample_columns].mean(axis=1)
	marker_coverage = marker_coverage.reset_index()

	# Calculate sample-level coverage
	# Sum read depths across all rows for each sample
	sample_coverage = df[sample_columns].sum(axis=0)
	sample_coverage['Mean_marker_coverage'] = df[sample_columns].mean(axis=0)

	# Output the results
	marker_coverage_file = madc.replace('.csv', '_marker_coverage.csv')
	marker_coverage.to_csv(marker_coverage_file, index=False)
	print("Number of samples:", len(sample_columns))
	print("Marker-Level Coverage:")
	print(marker_coverage)
	
	sample_coverage_file = madc.replace('.csv', '_sample_coverage.csv')
	sample_coverage.to_csv(sample_coverage_file, header=['Total_Coverage'])
	print("Number of markers:", len(marker_coverage))
	print("\nSample-Level Coverage:")
	print(sample_coverage)


if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="")
	
	parser.add_argument('madc', help='MADC report with alleles assigned fixed IDs')
	
	args = parser.parse_args()
	
	calculate_coverage(args.madc)
