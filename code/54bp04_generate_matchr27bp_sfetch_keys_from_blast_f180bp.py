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
def get_sfetch_keys(blast, blast_unique):
    outp = open(blast.replace('f180bp.bn.best', 'f180bp.bn.best_r27bp_sfetchKeys.txt'), 'w')
    for key, line_array in blast_unique.items():
        if int(line_array[3]) == int(line_array[1]):
            start_from = int(line_array[7]) + 1
        else:
            start_from = int(line_array[7]) + (int(line_array[1]) - int(line_array[3])) + 1
            #print('Alt: This allele does not align to end of sequence: ', line_array)

        newname = line_array[0] + '_27bp'
        end_to = start_from + 26
        if end_to > 361:
            end_to = 361
        outp.write('\t'.join([newname, str(start_from), str(end_to), line_array[4]]) + '\n')
    outp.close()


def get_query_unique_hits(blast):
    # VaccDscaff11_000042737|RefMatch_0001    54      1       54      VaccDscaff11_000042737|Ref      361     152     205     54      100     98.148  8.85e-23
    inp = open(blast)
    line = inp.readline()
    blast_unique = {}
    while line:
        line_array = line.strip().split()
        query = line_array[0].split('|')[0]
        subject = line_array[4].split('|')[0]
        if query == subject:
            if ('RefMatch' in line_array[0] and 'Ref' in line_array[4]) or ('AltMatch' in line_array[0] and 'Alt' in line_array[4]):
                # Add information to a dictionary and make sure every subject sequence only appear once
                cov = (int(line_array[3]) - int(line_array[2]) + 1)/float(line_array[1]) * 100
                if line_array[0] not in blast_unique:
                    if cov >= 90.0:
                        blast_unique[line_array[0]] = line_array
                    else:
                        # Remove alleles if coverage is lower than 90%
                        pass
                else:
                    # Compare length of coverage
                    cov_pre = (int(blast_unique[line_array[0]][3]) - int(blast_unique[line_array[0]][2]) + 1)/float(blast_unique[line_array[0]][1]) * 100
                    if cov > cov_pre:
                        print('  # Update this query because there is alignment with a higher COVERAGE: ')
                        print(cov, cov_pre)
                        print('    # current record: ', line_array)
                        print('    # existing record: ', blast_unique[line_array[0]])
                        blast_unique[line_array[0]] = line_array
                    elif cov == cov_pre:
                        # Compare alignment identity
                        if float(line_array[10]) > float(blast_unique[line_array[0]][10]):
                            print('  # Update this query because there is alignment with a higher IDENTITY: ', line_array)
                            blast_unique[line_array[0]] = line_array
                        else:
                            pass
                    else:
                        pass
            else:
                pass
        else:
            pass
        line = inp.readline()
    inp.close()

    # Sanity check
    for key, line_array in blast_unique.items():
        # Alignment should be all on the plus strand, check it here
        if int(line_array[6]) > int(line_array[7]):
            print('  # Alignment on MINUS strand: ', line_array)
        else:
            pass
    print('  # Number of RefMatch and AltMatch alleles: ', len(blast_unique))
    return(blast_unique)


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract fasta sequences of 3K DArTag panel")

    parser.add_argument('blast',
                        help='BLAST results of alleles against the f180-bp flanking sequences with the same orientation')

    args=parser.parse_args()

    blast_unique = get_query_unique_hits(args.blast)

    get_sfetch_keys(args.blast, blast_unique)
