#!/usr/bin/python3

def cal_sequence_length(fasta):
    inp = open(fasta)
    out = open(fasta + '.len', 'w')
    line = inp.readline()
    seq_len = {}
    seq_ID = ''
    seq = ''
    while line:
        line = line.strip()
        if line.startswith('>'):
            if seq != '':
                seq_len[seq_ID] = len(seq)
                out.write(seq_ID + '\t' + str(len(seq)) + '\n')
            else:
                pass
            seq_ID = line[1:]
            seq = ''
        else:
            seq += line
        line = inp.readline()
    inp.close()
    
    # Last sequence
    seq_len[seq_ID] = len(seq)
    out.write(seq_ID + '\t' + str(len(seq)) + '\n')
    out.close()
    return(seq_len)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('fasta',
                        help='A multiple-sequences fasta file')

    args=parser.parse_args()

    cal_sequence_length(args.fasta)
