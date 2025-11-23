#!/usr/bin/python3

def generate_allele_frequency_table(hapPA_file):
	import pandas as pd
	
	# Load the dataset
	df = pd.read_csv(hapPA_file, sep=",")
	
	# Extract the allele type from AlleleID (Ref, Alt, RefMatch, AltMatch, etc.)
	# Extract part after "|"
	df["AlleleType"] = df["AlleleID"].apply(lambda x: x.split("|")[1].split('_')[0])
	
	# Keep only the specified categories: Ref, Alt, RefMatch, and AltMatch
	df = df[df["AlleleType"].isin(["Ref", "Alt", "RefMatch", "AltMatch"])]
	
	# Group data by CloneID and AlleleType
	allele_frequencies = df.groupby(["CloneID", "AlleleType"])["Presence_count"].sum().reset_index()
	
	# Pivot the table for easy visualization
	# Pivot tables in pandas can produce columns in arbitrary orders
	allele_freq_table = allele_frequencies.pivot(index="CloneID", columns="AlleleType", values="Presence_count").fillna(0)
	
	# Reorder columns explicitly
	column_order = ["Ref", "Alt", "RefMatch", "AltMatch"]
	allele_freq_table = allele_freq_table[column_order]
	
	# Add a Total_Count column for overall presence counts
	allele_freq_table["Total_Count"] = allele_freq_table.sum(axis=1)
	
	# Rename columns for clarity
	allele_freq_table = allele_freq_table.rename(
		columns={"Ref": "Ref_Count", "Alt": "Alt_Count", "RefMatch": "RefMatch_Count", "AltMatch": "AltMatch_Count"})
	
	# Reset index and write to file if needed
	allele_freq_table = allele_freq_table.reset_index()
	print(allele_freq_table)
	
	# Save the allele frequency table
	outf = hapPA_file.replace('.csv', '_alleleFreq.csv')
	allele_freq_table.to_csv(outf, index=False)


if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="")
	
	parser.add_argument('hapPA_file', help='Microhaplotype presence/absence matrix file')
	
	args = parser.parse_args()
	
	generate_allele_frequency_table(args.hapPA_file)
