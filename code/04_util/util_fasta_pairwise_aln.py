#!/usr/bin/python3
# Updated on 1/25/2024: Ns in some allele sequences, added a condition to not consider those bases

import sys
import pandas as pd

def rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-': '-'}
    reverse_seq = seq[::-1] # reverse the string
    complement_seq = ''.join([complement[base] for base in reverse_seq])
    return(complement_seq)


def format_alignment_with_numbers(alignment):
    """Format alignment with base numbering"""
    seq1_aln, match_line, seq2_aln = alignment[:3]
    
    # Initialize variables
    line_length = 50
    total_length = len(seq1_aln)
    
    for i in range(0, total_length, line_length):
        # Get sequence chunks
        chunk_end = min(i + line_length, total_length)
        seq1_chunk = seq1_aln[i:chunk_end]
        match_chunk = match_line[i:chunk_end]
        seq2_chunk = seq2_aln[i:chunk_end]
        
        # Calculate actual base positions (excluding gaps)
        seq1_pos = len(seq1_aln[:chunk_end].replace('-', ''))
        seq2_pos = len(seq2_aln[:chunk_end].replace('-', ''))
        
        # Format output with numbers
        print(f"\nPosition {i + 1}-{chunk_end}:")
        print(f"Seq1  {seq1_chunk} {seq1_pos:4d}")
        print(f"      {match_chunk}")
        print(f"Seq2  {seq2_chunk} {seq2_pos:4d}")


def compare(markerID, seq1, seq2):
    try:
        from Bio import pairwise2
        from Bio.pairwise2 import format_alignment
        
        # (seqA, seqB, match_point, mismatch_point, gap_point, gap_extend_point)
        alignments = pairwise2.align.globalms(seq1, seq2, 1, -.5, -1.5, -.1)

        if not alignments:
            return None, "No alignments found"
        
        # Only look at the first alignment even if there are multiple alignments with the same score
        aln = alignments[0]
        aln_list = format_alignment(*aln).split("\n")
        # ['GAGGAAAAACACACACTGTATGATTTTGGAAACTCGACATAGGCCTATTGGAGG', '||||||||||||||||||||||||||.|||||||||||||||||||||||||||', 'GAGGAAAAACACACACTGTATGATTTCGGAAACTCGACATAGGCCTATTGGAGG', '  Score=52.5', '']
        '''
        seqA: --GAGGAAAAACACACACTGTATGATTTTGGAAACTCG-CATAGGCCTATTGGAGG
        anno:   |||||||||||||||| .|||||||||||||||||| |||||||||||||||||
        seqB: AAGAGGAAAAACACACAC-ATATGATTTTGGAAACTCGACATAGGCCTATTGGAGG
        In sequence, "-" denotes deletions
        In match-line, "|" denotes a match; "." denotes a mismatch; " " denotes a gap.
        '''

        # Get score
        score = float(aln_list[3].split("=")[1])
        indel_pos = 'None'
        if '-' in aln_list[0]:
            # If there is an indel, find the first position of the indel
            indel_pos = aln_list[0].find('-') # Provide the anchor base position
        elif '-' in aln_list[2]:
            if aln_list[2].count('-') == 1:
                indel_pos = aln_list[2].find('-') + 1 # for single-base deletion in the variant
                
            else:
                indel_pos = aln_list[2].find('-') # provide anchor base position for multi-base deletion
        else:
            indel_pos = 'None'
            
        # Print alignment info
        print("-" * 60)
        print(f"Marker ID: {markerID}")
        print(f"Alignment Score: {score}")
        print(f"Sequence 1 length: {len(seq1)}")
        print(f"Sequence 2 length: {len(seq2)}")
        print(f"First indel position: {indel_pos}")
        
        # Format and print alignment with numbers
        format_alignment_with_numbers(aln_list)
        print("-" * 60)
        
    except Exception as e:
        print(f"Error in alignment: {str(e)}")
        return None
 


def read_fasta_sequences(fasta):
    """
    Reads a FASTA file and returns a dictionary of records
    """
    from Bio import SeqIO
    records = {}
    for record in SeqIO.parse(fasta, "fasta"):
        allele_id = record.id
        clone_id = allele_id.split('|')[0]  # Extract clone ID from allele ID
        allele_sequence = str(record.seq)
        if clone_id in records:
            records[clone_id].append(allele_sequence)
        else:
            records[clone_id] = [allele_sequence]

    for key, value in records.items():
        if len(value) == 2:
            compare(key, value[0], value[1])
        else:
            print(f"Clone ID {key} has {len(value)} sequences, expected 2 for comparison.")
    



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate pairwsie alignment of alleles in a fasta file")

    parser.add_argument('fasta',
                        help='Multi-fasta file containing allele sequences')
    
    args=parser.parse_args()

    read_fasta_sequences(args.fasta)
