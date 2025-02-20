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
Ref alleles: locate beginning of each allele in 180 bp flanking sequences and extract 81 bp
Alt, RefMatch, and AltMatch alleles: locating the ending of each allele in 180 bp flanking sequences and extract 27 bp on the right hand side of the allele
'''
def get_sfetch_keys(blast):
    inp = open(blast)
    line = inp.readline()
    outp_ref = open(blast + '_sfetchKeys_ref.txt', 'w')
    outp_other = open(blast + '_sfetchKeys_match.txt', 'w')
    ref = []
    other = []
    while line:
        line_array = line.strip().split()
        query_alleleID_array = line_array[0].split('|')
        subject_alleleID_array = line_array[4].split('|')
        if query_alleleID_array[0] == subject_alleleID_array[0]:
            if line_array[0].endswith('Ref_0001'):
                if line_array[0] == 'VaccDscaff7_034027370|Ref_0001':
                    print(line_array)
                    break
                if not line_array[4].endswith('Ref'):
                    print('Ref: This allele did not align to Ref sequence: ', line_array)
                else:
                    pass

                # Confirm plus strand alignment
                if int(line_array[6]) < int(line_array[7]):
                    if int(line_array[2]) == 1:
                        newname = line_array[0]
                        start_from = int(line_array[6])
                        end_to = int(line_array[6]) + 80
                        if end_to > 361:
                            end_to = 361

                        if newname not in ref:
                            # sftech_key: <newname> <from> <to> <source seqname>
                            outp_ref.write('\t'.join([newname, str(start_from), str(end_to), line_array[4]]) + '\n')
                            ref.append(newname)
                        else:
                            print('This ref allele seems to be duplicated: ', newname)

                        # Just print this out for warning
                        if int(line_array[3]) != int(line_array[1]):
                            print('Ref: This allele does not align to end of sequence: ', line_array)
                    else:
                        print('Ref: This allele does not start alignment from 1: ', line_array)
                else:
                    print('Ref: Align on minus strand: ', line_array)
            else:
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

                    if newname not in other:
                        outp_other.write('\t'.join([newname, str(start_from), str(end_to), line_array[4]]) + '\n')
                        other.append(newname)
                    else:
                        print('This alt allele seems to be duplicated: ', newname)

                    if int(line_array[2]) != 1:
                        print('Alt: This allele does not start alignment from 1: ', line_array)
                else:
                    print('Alt: Align on minus strand: ', line_array)
        else:
            print('This allele did not align to the right marker locus: ', line_array)
        line = inp.readline()
    outp_ref.close()
    print('  # Number of ref alleles written out: ', len(ref))
    print('  # Number of other alleles written out: ', len(other))



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract fasta sequences of 3K DArTag panel")

    parser.add_argument('blast',
                        help='BLAST results of alleles against the f180-bp flanking sequences with the same orientation')

    args=parser.parse_args()

    get_sfetch_keys(args.blast)