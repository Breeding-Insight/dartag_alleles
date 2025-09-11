#!/usr/bin/python3


def collect_ref_alt_sequences(fasta1):
	"""
	Collects reference and alternative sequences from a multi-FASTA file.
	Assumes the reference sequence is labeled as 'Ref_0001' and the alternative sequence as 'Alt_0002'.
	Returns a dictionary with sequence IDs as keys and sequences as values.
	"""
	seq_id = None
	seq_fa = {}
	with open(fasta1, 'r') as f:
		for line in f:
			if line.startswith('>'):
				if seq_id:
					if 'Ref_0001' in seq_id or 'Alt_0002' in seq_id:
						seq_fa[seq_id] = seq
					else:
						pass
				seq_id = line.strip()[1:]  # Remove '>'
				seq = ''
			else:
				seq += line.strip()
		
		if 'Ref_0001' in seq_id or 'Alt_0002' in seq_id:
			seq_fa[seq_id] = seq
		else:
			pass
	f.close()
	return seq_fa



def update_ref_alt_sequences(fasta2, seq_fa, ref_len):
	outp = open(fasta2.replace('.fa', '_' + ref_len + 'bp.fa'), 'w')
	seq_id = None
	seq = ''
	cnt = 0
	present = []
	with open(fasta2, 'r') as f:
		for line in f:
			if line.startswith('>'):
				# Save previous sequence
				if seq_id:
					if seq_id in seq_fa:
						outp.write('>' + seq_id + '\n')
						outp.write(seq_fa[seq_id] + '\n')
						cnt += 1
					else:
						outp.write('>' + seq_id + '\n')
						outp.write(seq + '\n')
					present.append(seq_id)
				# Start new sequence
				seq_id = line.strip()[1:]  # Remove '>'
				seq = ''
			else:
				seq += line.strip()
		
		# Don't forget last sequence
		if seq_id in seq_fa:
			outp.write('>' + seq_id + '\n')
			outp.write(seq_fa[seq_id] + '\n')
			cnt += 1
		else:
			outp.write('>' + seq_id + '\n')
			outp.write(seq + '\n')
		present.append(seq_id)
	# Check for missing sequences
	new_cnt = 0
	for seq_id in seq_fa:
		if seq_id not in present:
			outp.write('>' + seq_id + '\n')
			outp.write(seq_fa[seq_id] + '\n')
			new_cnt += 1
	f.close()
	outp.close()
	print('# Ref and Alt sequences updated to', ref_len, ':', cnt)
	print('# Ref and Alt sequences not present in', fasta2, ':', new_cnt, '\n')



if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="Update Ref and Alt sequences to certain length")
	
	parser.add_argument('fasta1', help='Multi-fasta file containing longer allele sequences.')
	
	parser.add_argument('fasta2',
						help='Multi-fasta file containing shorter allele sequences to be updated with longer Ref and Alt sequences.')
	
	parser.add_argument('ref_len', help='Length of the reference sequence to be used for updating.')
	
	args = parser.parse_args()
	
	seq_fa = collect_ref_alt_sequences(args.fasta1)
	
	update_ref_alt_sequences(args.fasta2, seq_fa, args.ref_len)
