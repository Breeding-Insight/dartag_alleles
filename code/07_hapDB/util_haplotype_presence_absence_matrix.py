#!/usr/bin/python3

def get_project_haplotypes(madc, project_description, read_threshold):
	import pandas as pd
	df = pd.read_csv(madc, index_col='AlleleID')
	df = df.drop(['CloneID', 'AlleleSequence'], axis=1)
	# Select only numeric columns
	df = df.apply(pd.to_numeric, errors='coerce')
	df[df < int(read_threshold)] = 0
	df[df >= int(read_threshold)] = 1
	outf = madc.replace('.csv', '_hapStatus.csv')
	df.to_csv(outf, index_label='AlleleID')

	# Create a new DataFrame with rows where any value is greater than the threshold
	df_new = df[df.sum(axis=1) >= 1]
	print(madc.replace('.csv', '_hapPA.csv'))
	outf2 = open(madc.replace('.csv', '_hapPA.csv'), 'w')
	haps = df_new.index.tolist()
	outf2.write(','.join(["AlleleID", project_description]) + '\n')
	for i in haps:
		outf2.write(','.join([i, "1"]) + '\n')


if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="")
	
	parser.add_argument('madc',
						help='MADC with alleles assigned fixed IDs')
	
	parser.add_argument('project_description',
						help='Name of the PI to be added to the output file')
	
	parser.add_argument('--read_threshold',
						type=int,
						default=10,
						help='A read count threshold for considering allele presence/absence')	
	
	args = parser.parse_args()
	
	get_project_haplotypes(args.madc, args.project_description, args.read_threshold)
