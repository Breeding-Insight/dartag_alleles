#!/usr/bin/python3

from io import StringIO
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align.Applications import MuscleCommandline
from collections import defaultdict
import subprocess
import os


def check_first_20bp(alignment, marker_id, outf_summary, outf_exceptions, first_bp):
	results = {}
	
	outp_summary = open(outf_summary, 'a')
	outp_exceptions = open(outf_exceptions, 'a')
	# Get first 20 positions
	first_xbp = [record.seq[0:int(first_bp)] for record in alignment]
	
	# Check if all sequences are identical in this region
	is_identical = len(set(first_xbp)) == 1
	
	# Store result
	results[marker_id] = {'identical': is_identical, 'sequences': first_bp if not is_identical else first_bp[0]}
	
	# Print results for this alignment
	outp_summary.write(f"{marker_id},{is_identical}\n")
	print(f"\nMarker: {marker_id}")
	print(f"First 20bp identical: {is_identical}")
	if not is_identical:
		print("Variations found:")
		for i, seq in enumerate(first_bp):
			print(f"{alignment[i].id}: {seq}")
			outp_exceptions.write(f"{marker_id},{alignment[i].id},{seq}\n")
	
	return results


def find_first_snp(alignment, marker_id, outf_summary):
	'''
	Alignment with 6 rows and 81 columns
	CAAGACAAAGAGGTGGAGAGGGTTGCGTTGTTTGTTAGAAATAC...--- chr1.1_023060630|AltMatch_0003
	CAAGACAAAGAGGTGGAGAGGGTTGAGTTGTTTGTTAGAAATAC...--- chr1.1_023060630|AltMatch_0001
	CAAGACAAAGAGGTGGAGAGGATTGAGTTGTTTGTTAGAAATAT...--- chr1.1_023060630|AltMatch_0002
	CAAGACAAAGAGGTGGAGAGGGTTTACTTGTTTGTTAGAAATAT...--- chr1.1_023060630|RefMatch_0001
	CAAGACAAAGAGGTGGAGAGGGTTGAGTTGTTTGTTAGAAATAT...GGT chr1.1_023060630|Alt_0002
	CAAGACAAAGAGGTGGAGAGGGTTGACTTGTTTGTTAGAAATAT...GGT chr1.1_023060630|Ref_0001
	'''
	outp = open(outf_summary, 'a')
	
	# Get length of alignment
	aln_length = alignment.get_alignment_length()
	
	# Check each position
	for pos in range(aln_length):
		# Get all bases at this position
		bases = set(record.seq[pos] for record in alignment if record.seq[pos] != '-')

		# If more than one base type exists, it's a SNP
		if len(bases) > 1:
			# Determine type of variation
			if '-' in bases:
				var_type = "indel"
			else:
				var_type = "SNP"
			outp.write(f"{marker_id},{pos + 1},{var_type},','.join{bases}\n")
			print(f"First SNP found at position {pos + 1}")
			for record in alignment:
				print(f"{record.id}: {record.seq[pos]}")
			return pos + 1
	
	print("No SNPs found")
	outp.close()
	return None


def run_msa_for_marker(sequences, marker_id, output_dir, outf_summary, outf_exceptions, aln_suffix, first_bp):
	"""Run MSA for one marker's sequences"""
	# Create marker-specific fasta
	temp_fasta = os.path.join(output_dir, f"{marker_id}_{aln_suffix}_input.fasta")
	temp_aln = os.path.join(output_dir, f"{marker_id}_{aln_suffix}.aln")
	
	# Write sequences for this marker
	SeqIO.write(sequences, temp_fasta, "fasta")
	
	muscle_cline = MuscleCommandline(input=temp_fasta, out=temp_aln, clw=True, maxiters=2)
	subprocess.run(str(muscle_cline).split())
	
	# Write sequences for this marker
	temp_fasta = StringIO()
	SeqIO.write(sequences, temp_fasta, "fasta")
	fasta_data = temp_fasta.getvalue()
	
	# Run MUSCLE and capture output
	muscle_cline = MuscleCommandline()
	child = subprocess.Popen(str(muscle_cline),
							 stdin=subprocess.PIPE,
							 stdout=subprocess.PIPE,
							 stderr=subprocess.PIPE,
							 universal_newlines=True)
	
	# Get alignment from MUSCLE
	align_data, error = child.communicate(input=fasta_data)
	
	# Parse alignment
	alignment = AlignIO.read(StringIO(align_data), "fasta")
	
	#check_first_20bp(alignment, marker_id, outf_summary, outf_exceptions, first_bp)
	
	find_first_snp(alignment, marker_id, outf_summary)


def process_each_marker(marker_seqs, output_dir, outf_summary, outf_exceptions, aln_suffix, first_bp):
	# Process each marker
	print(f"Processing {len(marker_seqs)} markers...")
	for marker_id, sequences in marker_seqs.items():
		if len(sequences) > 1:  # Only align if multiple sequences exist
			#print(f"Aligning marker {marker_id} ({len(sequences)} sequences)")
			try:
				run_msa_for_marker(sequences, marker_id, output_dir, outf_summary, outf_exceptions, aln_suffix, first_bp)
			except Exception as e:
				print(f"Error processing marker {marker_id}: {e}")
		else:
			print(f"Skipping marker {marker_id}: only one sequence")
			
			

def generate_seqRecord(line_array, marker_seqs, allele_idx):
	cloneID = line_array[1]
	# Check if cloneID already exists for this marker
	if cloneID in marker_seqs:
		allele_idx += 1
	else:
		allele_idx = 1
	
	if allele_idx <= 10:
		seqID = f"s{allele_idx:04d}_{line_array[0]}"
		sequence = Seq(line_array[2])
		
		# Create a SeqRecord object
		record = SeqIO.SeqRecord(seq=sequence, id=seqID)
		
		marker_seqs[cloneID].append(record)
	else:
		pass
	return marker_seqs, allele_idx
	
	
	
def parse_fasta_by_marker(madc, ref_alt_only='N'):
	"""Group sequences by marker locus"""
	marker_seqs = defaultdict(list)
	
	# Parse MADC and group by marker
	inp = open(madc)
	line = inp.readline()
	allele_idx = 1
	while line:
		if line.startswith('*') or line.startswith('AlleleID'):
			pass
		else:
			line_array = line.strip().split(',')
			if ref_alt_only == 'Y':
				if line_array[0].endswith('|Ref') or line_array[0].endswith('|Alt'):
					marker_seqs, allele_idx = generate_seqRecord(line_array, marker_seqs, allele_idx)
				else:
					pass
			else:
				marker_seqs, allele_idx = generate_seqRecord(line_array, marker_seqs, allele_idx)
		line = inp.readline()
	
	return marker_seqs



			
			
			
if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="")
	
	parser.add_argument('madc', help='MADC file containing allele sequences')
	
	parser.add_argument('output_dir',
						help='Directory to save multiple sequence alignment output files')
	
	parser.add_argument('outf_prefix',
						help='Prefix for output files')
	
	parser.add_argument('-ref_alt_only', type=str, default='N')
	
	parser.add_argument('-aln_suffix',
						type=str,
						default='allAlleles')
	
	parser.add_argument('-first_bp',
						type=int,
						default=20,
						help='First xx base pair to check for identical sequences')
	
	args = parser.parse_args()
	
	# Parse sequences by marker
	marker_seqs = parse_fasta_by_marker(args.madc, args.ref_alt_only)
	
	# Process each marker
	outf_summary = args.outf_prefix + "_summary.csv"
	outf_exceptions = args.outf_prefix + "_exceptions.txt"
	process_each_marker(marker_seqs, args.output_dir, outf_summary, outf_exceptions, args.aln_suffix, args.first_bp)
