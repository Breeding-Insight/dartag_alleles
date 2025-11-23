#!/usr/bin/python3
# Things polyRAD changes in the VCF output after running readDArTag
# Amplicons from the bottom strand are reverse complemented
# The positions are changed from target SNP to the beginning of amplicons

from Bio import SeqIO
from Bio.Seq import Seq
import vcf

def create_seq_to_ID_dict(madc):
    # Create a dictionary with sequences as keys and IDs as values
    seq_to_ID = {}
    inp = open(madc, 'r')
    line = inp.readline()  # Skip the first line (header)
    line = inp.readline()  # Read the first sequence line
    while line:
        # Get original sequence and its reverse complement
        line_array = line.strip().split(',')
        orig_seq = line_array[2].upper()
        rev_comp = str(Seq(orig_seq).reverse_complement())
        # Store both original and reverse complement sequences with the same ID
        seq_to_ID[orig_seq] = line_array[0]
        seq_to_ID[rev_comp] = line_array[0]
        line = inp.readline()  # Read the next sequence line
    inp.close()
    return(seq_to_ID)


def get_bottom_loci(botloci):
    # Get bottom strand loci
    bottom_loci = []
    with open(botloci, 'r') as inp:
        for line in inp:
            bottom_loci.append(line.strip())
    return(bottom_loci)


def get_sub_sequence(sequence, bottom_loci, targetSNP_ID, first_bp):
    if targetSNP_ID not in bottom_loci:
        # Remove the first xx bp from REF sequence
        # Can't trim the end of the sequence because the variants might be in bases beyond 54 bp
        sequence = sequence[int(first_bp):]
    else:
        # on bottom strand, the REF sequence is reverse complemented
        # No need to +1 because python list is half-open interval
        end_bp = len(sequence) - int(first_bp)
        sequence = sequence[:end_bp]
    return sequence


def update_vcf_info(input_vcf, output_vcf, seq_to_ID, bottom_loci, first_bp):
    # Update VCF file by adding sequence IDs to INFO field
    vcf_reader = vcf.Reader(open(input_vcf, 'r'))

    # Remove SAMPLE entries from the Reader’s metadata dictionary
    vcf_reader.metadata.pop('SAMPLE', None)
    
    # Add new INFO field to the header
    # Number can be
    # A: one value per alt allele
    # R: one value per allele including ref
    # A fixed number: specifies how many values can be associated with the INFO tag
    vcf_reader.infos['targetSNP'] = vcf.parser._Info(id='targetSNP', num='1', type='String',
        desc='Target SNP position in the reference genome', source=None, version=None)
    
    vcf_reader.infos['DArTstrand'] = vcf.parser._Info(id='DArTstrand', num='1', type='String',
        desc='DArTag amplicon strand', source=None, version=None)
    
    vcf_reader.infos['hapID'] = vcf.parser._Info(id='hapID', num='R', type='String',
        desc='Haplotype IDs separated by comma', source=None, version=None)
        
    # Create VCF writer with updated header
    vcf_writer = vcf.Writer(open(output_vcf, 'w'), vcf_reader)
    for record in vcf_reader:
        allele_IDs = []
        targetSNP_ID = '.'
        
        # The coordinates in the VCF file are for the starting positions of the alleles on the top strand of the reference genome
        # These are computed by polyRAD based on target SNP positions and strand information of the alleles
        if str(record.REF) in seq_to_ID:
            allele_ID = seq_to_ID[str(record.REF)]
            targetSNP_ID = allele_ID.split('|')[0]  # Extract target SNP ID from allele ID
            record.ID = targetSNP_ID
            allele_IDs.append(seq_to_ID[str(record.REF)])
            record.REF = get_sub_sequence(str(record.REF), bottom_loci, targetSNP_ID, first_bp)
        else:
            print(f"Warning: REF sequence '{record.REF}' not found in the database. Using '.' as allele ID. {record.CHROM}:{record.POS}")
        
        # Adjust ALT sequences
        new_alts = []
        if record.ALT:
            #new_alt = get_sub_sequence(str(record.ALT[0]), bottom_loci, targetSNP_ID, first_bp)
            #new_alts.append(new_alt)
            #allele_IDs.append(targetSNP_ID + '|Alt_0002')
            index = 0
            while index < len(record.ALT):
                match_seq = str(record.ALT[index]).upper()
                if match_seq in seq_to_ID:
                    allele_IDs.append(seq_to_ID[match_seq])
                else:
                    allele_IDs.append('.')
                index += 1
                match_seq = get_sub_sequence(match_seq, bottom_loci, targetSNP_ID, first_bp)
                new_alts.append(match_seq)
            record.ALT = new_alts

        # Add all information to INFO field
        # Adjust position for top strand
        # No need to adjust position for bottom strand because the trimming is done on the 3' end of the bottom-strand sequence
        record.INFO['targetSNP'] = targetSNP_ID
        if targetSNP_ID in bottom_loci:
            record.INFO['DArTstrand'] = '-'
        else:
            record.INFO['DArTstrand'] = '+'
            record.POS += int(first_bp)  # Adjust position for top strand
        record.INFO['hapID'] = ','.join(allele_IDs)
        vcf_writer.write_record(record)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('MADC',
                        help='MADC file with fixed allele IDs and sequences')

    parser.add_argument('botloci',
                        help='A file containing the marker loci from the bottom strand')
    
    parser.add_argument('readdartag_vcf',
                        help='A readme file to add change information')

    parser.add_argument('--first_bp', type=int, default=0,
                        help='The first base pair to use for allele sequence extraction (default: 0)')

    args=parser.parse_args()

    # Use sequences as keys and allele IDs as values
    # Also get reverse complement sequences to capture those alleles on the bottom strand
    seq_to_ID = create_seq_to_ID_dict(args.MADC)

    # Get bottom strand loci
    bottom_loci = get_bottom_loci(args.botloci)
    
    # Update VCF file
    outp_vcf = args.readdartag_vcf.replace('.vcf', '_hapIDs.vcf')
    update_vcf_info(args.readdartag_vcf, outp_vcf, seq_to_ID, bottom_loci, args.first_bp)
