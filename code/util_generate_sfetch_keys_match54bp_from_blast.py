#!/usr/bin/python3
# Generate a key file for sfetch, then run the following
# esl-sfetch -Cf [options] seqfile subseq-coord-file
# (retrieve multiple subsequences using file of keys and coords)
# /programs/hmmer-3.1b1-linux-intel-x86_64/binaries/esl-sfetch -Cf input.fasta sfetch_keyfile.txt


'''
-C Multiple subsequence retrieval mode, with -f option (required).
Specifies that the second command line argument
is to be parsed as a subsequence coordinate file, consisting of lines containing four whitespace-delimited fields:
new_name, from, to, name/accession. For each such line, sequence
name/accession is found, a subsequence from..to is extracted,
and the subsequence is renamed new_name before being output.
Any other fields after the first four are ignored. Blank lines
and lines beginning with # are ignored.

In retrieving subsequences listed in a file (-C -f, or just -Cf), each line of the file
  is in GDF format: <newname> <from> <to> <source seqname>, space/tab delimited.
'''

'''
Generate keys for sfetch
RefMatch, and AltMatch alleles: locating the ending of each allele in 180 bp flanking sequences and extract 27 bp on the right hand side of the allele
'''
def get_sfetch_keys(blast):
    '''
    VaccDscaff11_000042737|RefMatch_0001    54      1       54      VaccDscaff11_000042737|Ref      361     152     205     54      100     98.148  8.85e-23
    VaccDscaff11_000042737|RefMatch_0001    54      1       54      VaccDscaff11_000042737|Alt      361     152     205     54      100     96.296  2.16e-20
    '''
    inp = open(blast)
    line = inp.readline()
    outp = open(blast + '_sfetchKeys.txt', 'w')
    alleles = []
    while line:
        line_array = line.strip().split()
        query_alleleID_array = line_array[0].split('|')
        subject_alleleID_array = line_array[4].split('|')
        if query_alleleID_array[0] == subject_alleleID_array[0]:
            # Confirm plus strand alignment
            if int(line_array[6]) < int(line_array[7]):
                if int(line_array[3]) == int(line_array[1]):
                    start_from = int(line_array[7]) + 1
                else:
                    start_from = int(line_array[7]) + (int(line_array[1]) - int(line_array[3])) + 1
                    print('Alt: This allele does not align to end of sequence: ', line_array)
                newname = line_array[0] + '_27bp'
                end_to = start_from + 26
                if end_to > 361:
                    end_to = 361

                if newname not in alleles:
                    outp.write('\t'.join([newname, str(start_from), str(end_to), line_array[4]]) + '\n')
                    alleles.append(newname)
                else:
                    print('This alt allele seems to be duplicated: ', newname)

                if int(line_array[2]) != 1:
                    print('Alt: This allele does not start alignment from 1: ', line_array)
            else:
                print('Alt: Align on minus strand: ', line_array)
        else:
            print('This allele did not align to the right marker locus: ', line_array)
        line = inp.readline()
    outp.close()
    print('  # Number of RefMatch and AltMatch alleles written out: ', len(alleles))



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract fasta sequences of 3K DArTag panel")

    parser.add_argument('blast',
                        help='BLAST results of alleles against the f180-bp flanking sequences with the same orientation')

    args=parser.parse_args()

    get_sfetch_keys(args.blast)