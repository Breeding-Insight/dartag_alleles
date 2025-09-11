#!/usr/bin/python3

# Plan:
# Read both CSV files using pandas
# Identify shared and unique markers between the files
# For shared markers, combine (add) the read counts
# Keep unique markers from both files
# Add a source column to indicate where unique markers came from
# Combine everything into a final dataframe

import pandas as pd


def sort_based_on_alleleID(df):
	# Sort the alleles from Ref, Alt, RefMatch, AltMatch
	df = df.set_index('AlleleID')
	df_groupby = df.groupby('CloneID')
	df_ordered = pd.DataFrame()
	for cloneID, clone_df in df_groupby:
		ref = cloneID + '|Ref_0001'
		alt = cloneID + '|Alt_0002'
		reindex = [ref, alt]
		idx_sorted = sorted(clone_df.index.to_list())
		for i in idx_sorted:
			if 'RefMatch' in i:
				reindex.append(i)
			else:
				pass
		
		for i in idx_sorted:
			if 'AltMatch' in i:
				reindex.append(i)
			else:
				pass
		clone_df = clone_df.reindex(reindex)
		df_ordered = pd.concat([df_ordered, clone_df], axis=0)
	return df_ordered
		
		
def combine_marker_counts(file1, file2, first_bp, outf):
	# 1. Read the CSV files
	df1 = pd.read_csv(file1)
	df2 = pd.read_csv(file2)
	
	# Assuming first column contains marker names
	# Store it separately and use rest for counts
	markers1 = df1.iloc[:, 0]
	markers2 = df2.iloc[:, 0]
	counts1 = df1.iloc[:, 3:]
	counts2 = df2.iloc[:, 3:]
	
	# 2. Identify shared and unique markers
	shared_markers = set(markers1) & set(markers2)
	unique_markers1 = set(markers1) - set(markers2)
	unique_markers2 = set(markers2) - set(markers1)
	
	# 3. Combine shared markers
	shared_df = pd.DataFrame()
	clone_ids = []
	allele_seqs = []
	for marker in shared_markers:
		idx1 = markers1[markers1 == marker].index[0]
		idx2 = markers2[markers2 == marker].index[0]
		# Int64Index([1], dtype='int64')
		# [0]: Takes the first (and should be only) index value
		
		# Store the 'CloneID' and AlleleSequence' values
		clone_ids.append(df1.iloc[idx1]['CloneID'])
		allele_seqs.append(df1.iloc[idx1]['AlleleSequence'])
		
		combined_counts = counts1.iloc[idx1] + counts2.iloc[idx2]
		shared_df = pd.concat([shared_df, pd.DataFrame([combined_counts])], ignore_index=True)
	
	# Add marker names to shared_df
	shared_df.insert(0, 'AlleleID', list(shared_markers))
	shared_df.insert(1, column='CloneID', value=clone_ids)
	shared_df.insert(2, column='AlleleSequence', value=allele_seqs)
	outf_shared = outf.replace('.csv', '_shared.csv')
	shared_df_ordered = sort_based_on_alleleID(shared_df)
	shared_df_ordered.to_csv(outf_shared, index=True)
	
	
	# 4. Get unique markers from file1
	# cols_to_use = [0] + list(range(3, len(df1.columns)))
	unique1_df = df1[df1.iloc[:, 0].isin(unique_markers1)].copy()
	
	
	# 5. Get unique markers from file2
	unique2_df = df2[df2.iloc[:, 0].isin(unique_markers2)].copy()
	
	
	# 6. Combine all dataframes
	final_df = pd.concat([shared_df, unique1_df, unique2_df], ignore_index=True)
	final_df = final_df.sort_values('AlleleID')
	df_ordered = sort_based_on_alleleID(final_df)
	df_ordered.to_csv(outf, index=True)
	
	
	
	# 7. Add source column to indicate where the markers came from
	file1_base = file1.split('/')[-1].replace('.csv', '')
	file2_base = file2.split('/')[-1].replace('.csv', '')
	shared_df['Source'] = 'Both'
	unique1_df['Source'] = file1_base
	unique2_df['Source'] = file2_base
	
	
	# 8. Combine all dataframes with source column
	final_df = pd.concat([shared_df, unique1_df, unique2_df], ignore_index=True)
	final_df = final_df.sort_values('AlleleID')
	final_df['AlleleSequence'] = final_df['AlleleSequence'].str.upper()
	final_df['AlleleSequence_rmF20bp'] = final_df['AlleleSequence'].str[int(first_bp):54]
	outf2 = outf.replace('.csv', '_with_source.csv')
	df_ordered = sort_based_on_alleleID(final_df)
	df_ordered.to_csv(outf2, index=True)
	
	# 9. Print out summary statistics
	print("# Summary of combined marker counts:")
	print(f"  # Common markers between {file1_base} and {file2_base}: {len(shared_markers)}")
	print(f"  # Unique markers in file1: {len(unique_markers1)}")
	print(f"  # Unique markers in file2: {len(unique_markers2)}")



if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="")
	
	parser.add_argument('file1', help='First MADC file with fixed allele IDs')
	
	parser.add_argument('file2', help='Second MADC file with fixed allele IDs')
	
	parser.add_argument('first_bp', help='First bp to use for allele sequence extraction (e.g., 20)')
	
	parser.add_argument('outf', help='Output file name for combined MADC counts')
	
	args = parser.parse_args()
	
	combine_marker_counts(args.file1, args.file2, args.first_bp, args.outf)
