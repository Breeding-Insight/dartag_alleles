#!/usr/bin/python3
# Generate a key file for sfetch, then run the following
# esl-sfetch -Cf [options] seqfile subseq-coord-file
# (retrieve multiple subsequences using file of keys and coords)
# /programs/hmmer-3.1b1-linux-intel-x86_64/binaries/esl-sfetch -Cf input.fasta sfetch_keyfile.txt


'''
-C Multiple subsequence retrieval mode, with -f option (required).
Specifies that the second command line argument
is to be parsed as a subsequence coordinate file, consisting of lines containing four whitespace-delimited fields:
new_name, from_bp, to, name/accession. For each such line, sequence
name/accession is found, a subsequence from_bp..to is extracted,
and the subsequence is renamed new_name before being output.
Any other fields after the first four are ignored. Blank lines
and lines beginning with # are ignored.

In retrieving subsequences listed in a file (-C -f, or just -Cf), each line of the file
  is in GDF format: <newname> <from_bp> <to> <source seqname>, space/tab delimited.
'''

def get_sequence_length(len_file):
    inp = open(len_file)
    line = inp.readline()
    seq_len = {}
    while line:
        line_array = line.strip().split()
        seq_len[line_array[0]] = int(line_array[1])
        line = inp.readline()
    inp.close()
    return(seq_len)


def get_sfetch_keys(chrPos, seq_len, f_len):
    # chr_1D_000510214
    inp = open(chrPos)
    line = inp.readline()
    outp = open(chrPos.replace('.csv', '_f' + f_len + 'bp_sfetchKeys.txt'), 'w')
    cnt = 0
    while line:
        line_array = line.strip().split(',')
        chrom = line_array[0]
        position = int(line_array[1])
        if position - int(f_len) > 0:
            from_bp = position - int(f_len)
        else:
            from_bp = 0
            
        if position + int(f_len) < seq_len[chrom]:
            to_bp = position + int(f_len)
        else:
            to_bp = seq_len[chrom]
        outp.write(' '.join([line_array[0] + '_' + line_array[1], str(from_bp), str(to_bp), chrom]) + '\n')
        cnt += 1
        line = inp.readline()
    inp.close()
    outp.close()
    print('# Number of keys written out: ', cnt)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Prepare sfetch keys")

    parser.add_argument('len_file',
                        help='File containing the length of the chromosomes')

    parser.add_argument('chrPos',
                        help='Marker position ')

    parser.add_argument('f_len',
                        help='Length of the flanking sequences to fetch')
    
    args=parser.parse_args()

    seq_len = get_sequence_length(args.len_file)
    
    get_sfetch_keys(args.chrPos, seq_len, args.f_len)
