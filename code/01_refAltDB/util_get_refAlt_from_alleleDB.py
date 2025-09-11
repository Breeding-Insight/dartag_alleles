#!/usr/bin/python3

def extract_ref_alt_sequences(fasta_file):
	seq_id = None
	seq_fa = {}
	seq_len = []
	seq = ''
	with open(fasta_file, 'r') as f:
		for line in f:
			if line.startswith('>'):
				# Save previous sequence
				if seq_id:
					if 'Ref_0001' in seq_id or 'Alt_0002' in seq_id:
						seq_fa[seq_id] = seq
						seq_len.append(len(seq))
					else:
						pass
				# Start new sequence
				seq_id = line.strip()[1:]  # Remove '>'
				seq = ''
			else:
				seq += line.strip()
		
		# Don't forget last sequence
		if 'Ref_0001' in seq_id or 'Alt_0002' in seq_id:
			seq_fa[seq_id] = seq
			seq_len.append(len(seq))
		else:
			pass
	if len(list(set(seq_len))) == 1:
		cnt = 0
		len_suffix = str(seq_len[0])
		# blueberry_allele_db_v019.fa
		outf = fasta_file.rsplit('_', 1)[0] + '_v000_refAlt_' + len_suffix + 'bp.fa'
		outp = open(outf, 'w')
		for seq_id, seq in seq_fa.items():
			outp.write('>' + seq_id + '\n')
			outp.write(seq + '\n')
			cnt += 1
		print(f'# From {fasta_file}')
		print(f'# Extracted {cnt} Ref and Alt sequences of length {len_suffix} bp')
	else:
		print('# Error: Sequences have different lengths, cannot extract Ref and Alt sequences.')



if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="Extract Ref and Alt sequences from a FASTA file containing allele sequences. The script assumes that the reference sequence is labeled as 'Ref_0001' and the alternative sequence as 'Alt_0002'. It will output a new FASTA file with these sequences if they are of the same length.")
	
	parser.add_argument('fasta', help='Multi-fasta file containing allele sequences to extract Ref and Alt sequences from. The reference sequence should be labeled as "Ref_0001" and the alternative sequence as "Alt_0002".')
	
	args = parser.parse_args()
	
	extract_ref_alt_sequences(args.fasta)
