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


def get_chr_len(chr_len_file):
    chr_len = {}
    inp = open(chr_len_file)
    line = inp.readline()
    while line:
        line_array = line.strip().split()
        chr_len[line_array[0]] = line_array[1]
        line = inp.readline()
    inp.close()
    return(chr_len)


def get_sfetch_keys(markers, chr_len, flankBP):
    outp = open(markers.replace('.csv', '') + '_f' + flankBP + 'bp_sfetchKeys.txt', 'w')
    out_count = 0
    inp = open(markers)
    line = inp.readline()
    # solcap_snp_c2_51460,chr01_000449027,chr01,449027,A,G
    while line:
        line_array = line.strip().split(',')
        if int(line_array[3]) < int(flankBP):
            start_from = 1
            end_to = int(flankBP)
            if line_array[2] in chr_len:
                if end_to > int(chr_len[line_array[2]]):
                    end_to = int(chr_len[line_array[2]])
                else:
                    end_to = int(line_array[3]) + int(flankBP) - 1
                # sftech_key: <newname> <from> <to> <source seqname>
                outp.write('\t'.join([line_array[1], str(start_from), str(end_to), line_array[2]]) + '\n')
                out_count += 1
        else:
            start_from = int(line_array[3]) - int(flankBP)
            end_to = int(line_array[3]) + int(flankBP)
            if line_array[2] in chr_len:
                if end_to > int(chr_len[line_array[2]]):
                    end_to = int(chr_len[line_array[2]])
                else:
                    end_to = int(line_array[3]) + int(flankBP)
                outp.write('\t'.join([line_array[1], str(start_from), str(end_to), line_array[2]]) + '\n')
                out_count += 1
        line = inp.readline()
    inp.close()
    outp.close()
    print('  # Total records written out: ', out_count)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract fasta sequences of 3K DArTag panel")

    parser.add_argument('markerIDs',
                        help='Marker ID in chrxx_000000000 format')

    parser.add_argument('chr_len',
                        help='Chromosome length file')
    
    parser.add_argument('flankBP',
                        help='Length of sequences to be extracted')

    args=parser.parse_args()
    
    chr_len = get_chr_len(args.chr_len)
    
    get_sfetch_keys(args.markerIDs, chr_len, args.flankBP)
