#!/usr/bin/python3



def generate_plots(file):
	import pandas as pd
	import matplotlib.pyplot as plt
	
	# Read CSV file
	df = pd.read_csv('/Users/dz359/PycharmProjects/BI/blueberry_dartag_validation/data/Report-DBlue21-6097/DBlue21-6097_Counts_missing_allele_discovery_rename_updatedSeq_db_vs_projectCnt.csv')
	
	# Group the data by chromosome
	grouped = df.groupby('Chromosome')
	
	# Create a bar graph for each chromosome
	for name, group in grouped:
		plt.bar(group['Chromosome'], group.iloc[:, 4], label='db_cnt')
		plt.bar(group['Chromosome'], group.iloc[:, 5], label='project_cnt')
		plt.title(f'Chromosome {name}')
		plt.xlabel('Chromosome')
		plt.ylabel('Data')
		plt.legend()
		plt.show()


if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="")
	
	parser.add_argument('file', help='')
	
	args = parser.parse_args()
	
	generate_plots('/Users/dz359/PycharmProjects/BI/blueberry_dartag_validation/data/Report-DBlue21-6097/DBlue21-6097_Counts_missing_allele_discovery_rename_updatedSeq_db_vs_projectCnt.csv')
