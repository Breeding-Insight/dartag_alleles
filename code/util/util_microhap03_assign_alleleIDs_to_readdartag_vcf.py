#!/usr/bin/python3
# Things polyRAD changes in the VCF output after running readDArTag
# Amplicons from the bottom strand are reverse complemented
# The positions are changed from target SNP to the beginning of amplicons

from Bio import SeqIO
from Bio.Seq import Seq
import vcf

def create_seq_to_ID_dict(db_fasta):
    # Create a dictionary with sequences as keys and IDs as values
    seq_to_ID = {}
    for record in SeqIO.parse(db_fasta, 'fasta'):
        # Get original sequence and its reverse complement
        orig_seq = str(record.seq).upper()
        rev_comp = str(Seq(orig_seq).reverse_complement())
        # Store both original and reverse complement sequences with the same ID
        seq_to_ID[orig_seq] = record.id
        seq_to_ID[rev_comp] = record.id 
    return(seq_to_ID)


def get_bottom_loci(botloci):
    # Get bottom strand loci
    bottom_loci = []
    with open(botloci, 'r') as inp:
        for line in inp:
            bottom_loci.append(line.strip())
    return(bottom_loci)


def update_vcf_info(input_vcf, output_vcf, seq_to_ID, bottom_loci):
    # Update VCF file by adding sequence IDs to INFO field
    vcf_reader = vcf.Reader(open(input_vcf, 'r'))
    
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
        
        # Check reference allele
        ref_seq = str(record.REF).upper()
        if ref_seq in seq_to_ID:
            allele_IDs.append(seq_to_ID[ref_seq])
            targetSNP_ID = seq_to_ID[ref_seq].split('|')[0]
        else:
            allele_IDs.append('.')
        
        # Check alternate alleles
        if record.ALT:
            for alt in record.ALT:
                alt_seq = str(alt).upper()
                if alt_seq in seq_to_ID:
                    allele_IDs.append(seq_to_ID[alt_seq])
                else:
                    allele_IDs.append('.')

        # Add all information to INFO field
        record.INFO['targetSNP'] = targetSNP_ID
        if targetSNP_ID in bottom_loci:
            record.INFO['DArTstrand'] = '-'
        else:
            record.INFO['DArTstrand'] = '+'
        record.INFO['hapID'] = ','.join(allele_IDs)
        vcf_writer.write_record(record)    


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('db_fasta',
                        help='Micorhaplotype db fasta file')

    parser.add_argument('botloci',
                        help='A file containing the marker loci from the bottom strand')
    
    parser.add_argument('readdartag_vcf',
                        help='A readme file to add change information')

    args=parser.parse_args()

    # Get reverse complement sequences to capture those alleles on the bottom strand
    seq_to_ID = create_seq_to_ID_dict(args.db_fasta)

    # Get bottom strand loci
    bottom_loci = get_bottom_loci(args.botloci)
    
    # Update VCF file
    outp_vcf = args.readdartag_vcf.replace('.vcf', '_hapIDs.vcf')
    update_vcf_info(args.readdartag_vcf, outp_vcf, seq_to_ID, bottom_loci)
