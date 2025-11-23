#!/usr/bin/python3

import pandas as pd

def aggregate_hapdb_statistics(input_paths_file, output_prefix):
	# Initialize lists to store extracted data for all files
	chromosome_data = []
	allele_length_stats = []
	allele_stats_per_locus = []
	general_summary = []
	
	# Read all file paths from the text file
	with open(input_paths_file, "r") as f:
		csv_file_paths = [line.strip() for line in f.readlines()]  # Read and clean each line
	
	# Process each CSV file
	for file_path in csv_file_paths:
		# /Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v050_summary.csv
		# Extract the first word of the file name (for naming columns)
		identifier = file_path.split("/")[-1].split("_")[0]  # Use the file name (e.g., "DB1" from "DB1_summary.csv")
		# Get the db version
		version = file_path.split("/")[-1].split("_")[-2]  # Use the second part of the file name (e.g., "v1" from "DB1_v1_summary.csv")
		
		# Read the file line by line
		with open(file_path, "r") as file:
			lines = file.readlines()
		
		# Extract general statistics
		total_sequences = None
		for line in lines:
			if "Total sequences" in line:
				total_sequences = int(line.split(",")[-1])
				general_summary.append({"Database": identifier, "Version": version, "Total Sequences": total_sequences})
		
		# Extract "Distribution per chromosome"
		start_index = None
		for i, line in enumerate(lines):
			if "Distribution per chromosome" in line:
				start_index = i + 2  # Skip column headers
				break
		if start_index is not None:
			for line in lines[start_index:]:
				if line.strip() == "" or "Allele Length Statistics" in line:
					# Stop when the next section begins or empty line is reached
					break
				parts = line.split(',')
				if len(parts) == 3:  # Chromosome, Number of Loci, Number of microhaplotypes
					chromosome_data.append({"Database": identifier, "Chromosome": parts[0], "Number of Loci": int(parts[1]),
						"Number of Microhaplotypes"   : int(parts[2])})
		
		# Extract "Allele Length Statistics (in base pairs)"
		start_index = None
		for i, line in enumerate(lines):
			if "Allele Length Statistics" in line:
				start_index = i + 2  # Skip column headers
				break
		if start_index is not None:
			for line in lines[start_index:]:
				if line.strip() == "" or "Allele statistics per locus" in line:
					# Stop when the next section begins
					break
				parts = line.split()
				allele_length_stats.append({"Database": identifier, "Statistic": ' '.join(parts[:-1]), "Value": int(float(parts[-1]))})
		
		# Extract "Allele statistics per locus"
		start_index = None
		for i, line in enumerate(lines):
			if "Allele statistics per locus" in line:
				start_index = i + 2  # Skip column headers
				break
		if start_index is not None:
			for line in lines[start_index:]:
				if line.strip() == "" or "Loci statistics per chromosome" in line:
					# Stop when the next section begins
					break
				parts = line.split(',')
				if len(parts) == 2:  # Statistic and Value
					allele_stats_per_locus.append({"Database": identifier, "Statistic": parts[0], "Value": float(parts[1])})
	
	# Create DataFrames for each section
	chromosome_df = pd.DataFrame(chromosome_data)
	allele_length_stats_df = pd.DataFrame(allele_length_stats)
	allele_stats_per_locus_df = pd.DataFrame(allele_stats_per_locus)
	general_summary_df = pd.DataFrame(general_summary)
	
	# Save the DataFrames to CSV files
	chromosome_df.to_csv(output_prefix + "_chromosome_distribution.csv", index=False)
	allele_length_stats_df.to_csv(output_prefix + "_allele_length_statistics.csv", index=False)
	allele_stats_per_locus_df.to_csv(output_prefix + "_allele_stats_per_locus.csv", index=False)
	general_summary_df.to_csv(output_prefix + "_general_summary.csv", index=False)


if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="")
	
	parser.add_argument('file_list', help='A file containing paths to CSV files to aggregate statistics from.')
	
	parser.add_argument('output_prefix', help='Output prefix for the aggregated CSV files.')
	
	args = parser.parse_args()
	
	aggregate_hapdb_statistics(args.file_list, args.output_prefix)
