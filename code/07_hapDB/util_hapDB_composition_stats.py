#!/usr/bin/python3

import os
import csv
from collections import defaultdict
from datetime import datetime
import pandas as pd

# Function to parse the file
def parse_fasta(file):
	total_sequences = 0
	chromosomes = defaultdict(int)
	loci_by_chromosome = defaultdict(set)  # To store unique loci for each chromosome
	loci_alleles = defaultdict(set)
	allele_lengths = [] # To store lengths of alleles
	
	with open(file, 'r') as f:
		current_allele_seq = None
		for line in f:
			if line.startswith('>'):  # Only process the header lines
				total_sequences += 1
				if current_allele_seq is not None:
					allele_lengths.append(len(current_allele_seq)) # Store length of the previous allele
				header = line.strip()[1:]  # Remove the '>' symbol
				chrom_locus, allele = header.split('|')  # Split into chrom_locus and allele
				chromosome = '_'.join(chrom_locus.split('_')[:-1])  # Extract chromosome (e.g., "chr8.1")
				chromosomes[chromosome] += 1  # Increment chromosome count
				loci_by_chromosome[chromosome].add(chrom_locus)  # Add locus to the set for its chromosome
				loci_alleles[chrom_locus].add(allele)  # Add allele to its respective locus
				current_allele_seq = ''  # Reset for the new allele sequence
			else:
				if current_allele_seq is None:
					current_allele_seq = line.strip()
				else:
					current_allele_seq += line.strip() # Handle multi-line sequences
		# Don't forget to add the last allele sequence length
	if current_allele_seq is not None:
		allele_lengths.append(len(current_allele_seq))
		
	return total_sequences, chromosomes, loci_by_chromosome, loci_alleles, allele_lengths


# Function to calculate allele statistics using `pandas`
def calculate_allele_length_stats_pandas(allele_lengths):
    if not allele_lengths:
        return pd.DataFrame({"Stat": [], "Value": []})  # Return an empty DataFrame if no data exists

    # Create a series for allele lengths
    length_series = pd.Series(allele_lengths)

    # Use `pandas.describe()` and add specific calculations for std deviation and median
    stats = {
        "Total Alleles": len(length_series),
        "Min Length": length_series.min(),
        "Max Length": length_series.max(),
        "Mean Length": length_series.mean(),
        "Median Length": length_series.median(),
        "Standard Deviation": length_series.std()
    }

    # Convert the stats dictionary into a DataFrame for flexibility in output
    stats_df = pd.DataFrame(list(stats.items()), columns=["Stat", "Value"])
    return stats_df



# Function to generate Mb bins
def get_mb_bin(position):
	position_int = int(position)
	mb_bin = position_int // 1_000_000  # Determine the megabase bin
	return f"{mb_bin}"


# Function to write stats to a CSV
def write_csv(loci_alleles, file_wo_ext):
	output_csv = f"{file_wo_ext}_allele_stats.csv"
	
	# Collect all the data into a list
	data = []
	for locus, alleles in loci_alleles.items():
		chromosome = '_'.join(locus.split('_')[:-1])  # Extract chromosome (e.g., "chr8.1")
		position = locus.split('_')[-1]  # Extract position (e.g., "086468647")
		mb_bin = get_mb_bin(position)
		
		# Append the row data as a tuple (chromosome, position, and other fields)
		data.append({'Locus': locus, 'Chromosome': chromosome, 'Position': position, 'Mb_Bin': mb_bin,'Number_Alleles': len(alleles)})
	
	# Convert to DataFrame and save as CSV after sorting
	df = pd.DataFrame(data)
	df = df.sort_values(by=['Chromosome', 'Position']) # sort by chromosome and position
	df.to_csv(output_csv, index=False)
	print(f"CSV file '{output_csv}' written successfully.")


# Main function to calculate statistics and generate CSV
def calculate_statistics(fasta_file, file_wo_ext):
	total_sequences, chromosomes, loci_by_chromosome, loci_alleles, allele_lengths = parse_fasta(fasta_file)
	
	# Calcuate allele length statistics
	allele_length_stats_df = calculate_allele_length_stats_pandas(allele_lengths)
	
	# Generate allele data as a DataFrame
	allele_data = pd.DataFrame({'Locus': loci_alleles.keys(),
								'Number_Alleles': [len(alleles) for alleles in loci_alleles.values()]})
	
	# Summary statistics for alleles per locus
	allele_stats = allele_data['Number_Alleles'].describe()
	
	# Generate chromosome data
	chromosome_data = pd.DataFrame({'Chromosome': loci_by_chromosome.keys(),
									'Number_Loci': [len(loci) for loci in loci_by_chromosome.values()]})
	
	# Summary statistics for loci per chromosome
	chromosome_stats = chromosome_data['Number_Loci'].describe()
	
	time_now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	outf = file_wo_ext + "_summary.csv"
	
	with open(outf, 'w') as outp:
		# Basic information
		outp.write(f"Statistics generated on, {time_now}\n\n")
		outp.write(f"Total sequences, {total_sequences}\n\n")
		print(f"Total sequences: {total_sequences}")
		
		# Chromosome-wise distribution
		outp.write("Distribution per chromosome:\n")
		outp.write(f'Chromosome, Number of loci, Number of microhaplotypes\n')
		print("\nDistribution per chromosome:")
		for chromosome in sorted(chromosomes.keys()):
			num_loci = len(loci_by_chromosome[chromosome])
			num_sequences = chromosomes[chromosome]
			outp.write(f"  {chromosome}, {num_loci}, {num_sequences}\n")
			print(f"  {chromosome}: {num_loci} \t {num_sequences}")
			
		# Write allele length statistics
		outp.write("\nAllele Length Statistics (in base pairs):\n")
		outp.write(allele_length_stats_df.to_string(index=False, header=True) + "\n")
			
		# Allele statistics
		outp.write("\nAllele statistics per locus:\n")
		outp.write(allele_stats.to_csv(header=True) + '\n')
		
		# Chromosome statistics
		outp.write("\nLoci statistics per chromosome:\n")
		outp.write(chromosome_stats.to_csv(header=True) + '\n')
	
	print("\nSummary written to:", outf)
	print("\nSummary statistics:\n")
	print("Allele Statistics (Per Locus):\n", allele_stats)
	print("Chromosome Statistics (Loci Count):\n", chromosome_stats)
	
	# Write details to CSV
	write_csv(loci_alleles, file_wo_ext)


if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="")
	
	parser.add_argument('fasta',
						help='Multi-fasta sequence file with headers in the format: >chr8.1_086468647|allele1')
	
	args = parser.parse_args()
	
	file_wo_ext = os.path.splitext(args.fasta)[0]
	
	calculate_statistics(args.fasta, file_wo_ext)
