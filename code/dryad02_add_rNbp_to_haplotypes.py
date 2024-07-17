#!/usr/bin/python3
# Generate a key file for sfetch, then run the following
# esl-sfetch -Cf [options] seqfile subseq-coord-file
# (retrieve multiple subsequences using file of keys and coords)
# /programs/hmmer-3.1b1-linux-intel-x86_64/binaries/esl-sfetch -Cf input.fasta sfetch_keyfile.txt

'''
Extend alleles from 54 bp to 81 bp
Ref alleles: locate beginning of each allele in 180 bp flanking sequences and extract 81 bp
Alt, RefMatch, and AltMatch alleles: locating the ending of each allele in 180 bp flanking sequences and extract 27 bp on the right hand side of the allele
'''

def sort_alleles_and_output(extend81bp_allele_seq_remove, outf):
    outp = open(outf, 'w')
    # Sort the fasta sequences by ID; e.g., Ref_0001, Alt_0002, RefMatch_0001, AltMatch_0001,...
    locus_haplotypes = {}
    for key in sorted(extend81bp_allele_seq_remove):
        # chr1.1_000194324|AltMatch_0004
        marker = key.split('|')[0]
        haplotype = key.split('|')[1]
        if marker in locus_haplotypes:
            locus_haplotypes[marker].append(haplotype)
        else:
            locus_haplotypes[marker] = [haplotype]
            
    for marker in sorted(locus_haplotypes):
        ref = marker + '|Ref_0001'
        alt = marker + '|Alt_0002'
        outp.write(ref + '\n' + extend81bp_allele_seq_remove[ref] + '\n')
        outp.write(alt + '\n' + extend81bp_allele_seq_remove[alt] + '\n')
        locus_haplotypes[marker].remove('Ref_0001')
        locus_haplotypes[marker].remove('Alt_0002')
        for i in sorted(locus_haplotypes[marker]):
            match = marker + '|' + i
            outp.write(match + '\n' + extend81bp_allele_seq_remove[match] + '\n')
    outp.close()
    print('## Number of alleles written to the output: ', len(extend81bp_allele_seq_remove))

def remove_alleles_with_low_coverage(extend81bp_allele_seq, remove_alleles):
    extend81bp_allele_seq_remove = extend81bp_allele_seq
    inp = open(remove_alleles)
    cnt = 0
    for line in inp.readlines():
        del extend81bp_allele_seq_remove['>' + line.strip()]
        cnt += 1
    inp.close()
    print('## Number of alleles removed due to low sequence coverage: ', cnt)
    return(extend81bp_allele_seq_remove)
    
    
    
def add_Nbp_to_54bp_alleles(current_allele_seq, rNbp):
    inp = open(rNbp)
    line = inp.readline()
    '''
    >alfalfaRep2vsXJDY1_shared_1000570|AltMatch_0001_Nbp
    TCAAATCTTAGTCCTGTGTAGGACATC
    '''
    cnt = 0
    first_seq = 'true'
    extend81bp = ''
    extend81bp_allele_seq = current_allele_seq
    while line:
        if line.startswith('>'):
            if first_seq == 'true':
                first_seq = 'false'
            else:
                if fastaID in current_allele_seq:
                    extend81bp_allele_seq[fastaID] = extend81bp
                    cnt += 1
                else:
                    print('  This allele does not exist in the 54 bp db: ', fastaID)
            fastaID = line.strip().replace('_rNbp', '')
            extend81bp = ''
        else:
            line = line.strip()
            extend81bp = current_allele_seq[fastaID] + line
        line = inp.readline()
        
    # add last sequence
    if fastaID in current_allele_seq:
        current_allele_seq[fastaID] = extend81bp
        cnt += 1
    else:
        print('  This allele does not exist in the 54 bp db: ', fastaID)
    inp.close()
    import datetime
    now = datetime.datetime.now()
    nowf = now.strftime("%Y-%m-%d, %H:%M:%S")
    print(nowf)
    print('# Running dryad02_add_rNbp_to_haplotypes.py')
    print('## Number of current alleles: ', len(current_allele_seq))
    print('## Number of alleles with sequences extended to 81 bp: ', cnt)
    return(extend81bp_allele_seq)


def add_current_alleles_to_dict(current_alleles):
    inp = open(current_alleles)
    line = inp.readline()
    current_allele_seq = {}
    first = 'True'
    while line:
        if line.startswith('>'):
            if first == 'True':
                first = 'False'
            else:
                current_allele_seq[fastaID] = seq
            fastaID = line.strip()
            seq = ''
        else:
            seq += line.strip()
        line = inp.readline()
    # last sequence
    current_allele_seq[fastaID] = seq
    inp.close()
    return(current_allele_seq)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract fasta sequences of 3K DArTag panel")

    parser.add_argument('alleles_54bp',
                        help='DArTag alleles in 54 bp')

    parser.add_argument('rightNbp',
                        help='N bp on the right hand side of the alleles')

    parser.add_argument('remove_alleles', help='Alleles need to be removed from microhaplotype db due to low coverage')

    args=parser.parse_args()

    current_allele_seq = add_current_alleles_to_dict(args.alleles_54bp)

    extend81bp_allele_seq = add_Nbp_to_54bp_alleles(current_allele_seq, args.rightNbp)

    extend81bp_allele_seq_remove = remove_alleles_with_low_coverage(extend81bp_allele_seq, args.remove_alleles)

    outf = args.alleles_54bp.replace('.fa', '_81bp.fa')
    sort_alleles_and_output(extend81bp_allele_seq_remove, outf)
