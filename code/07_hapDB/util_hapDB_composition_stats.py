#!/usr/bin/python3

import os
import csv
from collections import defaultdict
from datetime import datetime

# Function to parse the file
def parse_fasta(file):
	total_sequences = 0
	chromosomes = defaultdict(int)
	loci_by_chromosome = defaultdict(set)  # To store unique loci for each chromosome
	loci_alleles = defaultdict(set)
	
	with open(file, 'r') as f:
		for line in f:
			if line.startswith('>'):  # Only process the header lines
				total_sequences += 1
				header = line.strip()[1:]  # Remove the '>' symbol
				chrom_locus, allele = header.split('|')  # Split into chrom_locus and allele
				chromosome = '_'.join(chrom_locus.split('_')[:-1])  # Extract chromosome (e.g., "chr8.1")
				chromosomes[chromosome] += 1  # Increment chromosome count
				loci_by_chromosome[chromosome].add(chrom_locus)  # Add locus to the set for its chromosome
				loci_alleles[chrom_locus].add(allele)  # Add allele to its respective locus
	return total_sequences, chromosomes, loci_by_chromosome, loci_alleles


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
	
	# Sort the data primarily by chromosome (alphanumeric), then by position (numeric)
	sorted_data = sorted(data, key=lambda x: (x['Chromosome'], int(x['Position'])))
	
	# Write the sorted data to a CSV file
	with open(output_csv, 'w', newline='') as csvfile:
		fieldnames = ['Locus', 'Chromosome', 'Position', 'Mb_Bin', 'Number_Alleles']
		writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
		writer.writeheader()
		
		# Write each row in the sorted order
		writer.writerows(sorted_data)


# Main function to calculate statistics and generate CSV
def calculate_statistics(fasta_file, file_wo_ext):
	total_sequences, chromosomes, loci_by_chromosome, loci_alleles = parse_fasta(fasta_file)
	time_now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	
	outp = open(file_wo_ext + "_summary.txt", 'w')
	outp.write(f"Statistics generated on, {time_now}\n\n")
	outp.write(f"Total sequences, {total_sequences}\n\n")
	outp.write("Distribution per chromosome:\n")
	print(f"Total sequences: {total_sequences}")
	print("\nDistribution per chromosome:")
	for chromosome in sorted(chromosomes.keys()):
		count = chromosomes[chromosome]
		loci_count = len(loci_by_chromosome[chromosome])  # Count unique loci for this chromosome
		outp.write(f"  {chromosome}, {loci_count}, {count}\n")
		print(f"  {chromosome}: {loci_count} \t {count}")
	
	print("\nWriting output to CSV file...")
	write_csv(loci_alleles, file_wo_ext)



if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="")
	
	parser.add_argument('fasta',
						help='Multi-fasta sequence file with headers in the format: >chr8.1_086468647|allele1')
	
	args = parser.parse_args()
	
	file_wo_ext = os.path.splitext(args.fasta)[0]
	
	calculate_statistics(args.fasta, file_wo_ext)
