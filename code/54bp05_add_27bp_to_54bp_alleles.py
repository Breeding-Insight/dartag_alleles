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

def add_27bp_to_54bp_alleles(allele_54bp_seq, right27bp, outf):
    inp = open(right27bp)
    line = inp.readline()
    '''
    >alfalfaRep2vsXJDY1_shared_1000570|AltMatch_0001_27bp
    TCAAATCTTAGTCCTGTGTAGGACATC
    '''
    outp = open(outf, 'w')
    cnt = 0
    while line:
        if line.startswith('>'):
            fastaID = line.strip().replace('_27bp', '')
            if fastaID in allele_54bp_seq:
                r27bp = inp.readline()
                seq = allele_54bp_seq[fastaID] + r27bp
                outp.write(fastaID + '\n' + seq)
                cnt += 1
            else:
                print('  This allele does not exist in the 54 bp db: ', fastaID)
        else:
            print('Check this line', line)
        line = inp.readline()
    inp.close()
    outp.close()
    print('  # Number of initial RefMatch and AltMatch alleles: ', len(allele_54bp_seq))
    print('  # Note that some alleles were removed because of low coverage/identity during BLAST search against the f180bp sequences')
    print('  # Number of RefMatch and AltMatch alleles written out: ', cnt)


def add_54bp_alleles_to_dict(alleles_54bp):
    inp = open(alleles_54bp)
    line = inp.readline()
    allele_54bp_seq = {}
    first = 'True'
    while line:
        if line.startswith('>'):
            if first == 'True':
                first = 'False'
            else:
                allele_54bp_seq[fastaID] = seq
            fastaID = line.strip()
            seq = ''
        else:
            seq += line.strip()
        line = inp.readline()
    # last sequence
    allele_54bp_seq[fastaID] = seq
    inp.close()
    return(allele_54bp_seq)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract fasta sequences of 3K DArTag panel")

    parser.add_argument('alleles_54bp',
                        help='DArTag alleles in 54 bp')

    parser.add_argument('right27bp',
                        help='27 bp on the right hand side of the 54 bp alleles')

    args=parser.parse_args()

    allele_54bp_seq = add_54bp_alleles_to_dict(args.alleles_54bp)

    outf = args.alleles_54bp.replace('.fa', '_81bp.fa')

    add_27bp_to_54bp_alleles(allele_54bp_seq, args.right27bp, outf)