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
Ref blast_unique: locate beginning of each allele in 180 bp flanking sequences and extract 81 bp
Alt, RefMatch, and AltMatch blast_unique: locating the ending of each allele in 180 bp flanking sequences and extract 27 bp on the right hand side of the allele
'''
def get_sfetch_keys(blast, blast_unique):
    outp = open(blast + '_sfetchKeys.txt', 'w')
    out_count = 0
    for key, line_array in blast_unique.items():
        # Confirm plus strand alignment
        if int(line_array[6]) < int(line_array[7]):
            if int(line_array[2]) == 1:
                start_from = int(line_array[6])
                end_to = int(line_array[6]) + 80
                if end_to > 361:
                    end_to = 361
                # sftech_key: <newname> <from> <to> <source seqname>
                outp.write('\t'.join([line_array[0], str(start_from), str(end_to), line_array[4]]) + '\n')
                out_count += 1
            else:
                print('  # This allele does not start alignment from 1, which is probably due to IUPAC codes, will extend x-bp based on the alignment: ', line_array)
                start_from = int(line_array[6]) - int(line_array[2]) + 1
                end_to = start_from + 80
                if end_to > 361:
                    end_to = 361
                outp.write('\t'.join([line_array[0], str(start_from), str(end_to), line_array[4]]) + '\n')
                out_count += 1
        else:
            print('  # Align on minus strand: ', line_array)
    outp.close()
    print('  # Total records written out: ', out_count)


def get_query_unique_hits(blast):
    # VaccDscaff11_000042737|Ref_0001    54      1    54      VaccDscaff11_000042737|Ref      361     210     157     54      100     100.000 3.63e-25
    inp = open(blast)
    line = inp.readline()
    blast_unique = {}
    while line:
        line_array = line.strip().split()
        query = line_array[0].rsplit('_', 1)[0]
        subject = line_array[4]
        if query == subject:
            if (line_array[0].endswith('Ref_0001') and line_array[4].endswith('Ref')) or (line_array[0].endswith('Alt_0002') and line_array[4].endswith('Alt')):
                # Add information to a dictionary and make sure every subject sequence only appear once
                if line_array[0] not in blast_unique:
                    blast_unique[line_array[0]] = line_array 
                else:
                    # Some Ref/Alt alleles contain IUPAC code, therefore, there may not be 100% identity and coverage BLAST hits
                    # Compare length of coverage
                    query_cov_inDict = abs(int(blast_unique[line_array[0]][3]) - int(blast_unique[line_array[0]][2])) + 1
                    query_cov = abs(int(line_array[3]) - int(line_array[2])) + 1
                    if query_cov > query_cov_inDict:
                        print('  # Update this query: \n', blast_unique[line_array[0]])
                        print('  # With this one:', line_array)
                        blast_unique[line_array[0]] = line_array
                    elif query_cov == query_cov_inDict:
                        # Compare alignment identity
                        if float(line_array[10]) > float(blast_unique[line_array[0]][10]):
                            print('  # Update this query: \n', blast_unique[line_array[0]])
                            print('  # With this one:', line_array)
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
    ref_count = alt_count = other_count = 0
    for key, line_array in blast_unique.items():
        # Alignment should be all on the plus strand, check it here
        if int(line_array[6]) > int(line_array[7]):
            print('  # Alignment on MINUS strand: ', line_array)
        else:
            pass

        if key.endswith('Ref_0001'):
            ref_count += 1
        elif key.endswith('Alt_0002'):
            alt_count += 1
        else:
            other_count += 1
    print('  # Extract unique hits for queries')
    print('     # Number of ref blast_unique: ', ref_count)
    print('     # Number of alt blast_unique: ', alt_count)
    return(blast_unique)


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract fasta sequences of 3K DArTag panel")

    parser.add_argument('blast',
                        help='BLAST results of blast_unique against the f180-bp flanking sequences with the same orientation')

    args=parser.parse_args()
    
    print('  # Running db07_generate_ref_alt_sfetch_keys_from_blast.py on', args.blast)

    blast_unique = get_query_unique_hits(args.blast)

    get_sfetch_keys(args.blast, blast_unique)
