#!/usr/bin/python3

def collect_fasta_sequences(fasta):
    inp = open(fasta)
    line = inp.readline()
    fasta_seq = {}
    seq = ''
    while line:
        if line.startswith('>'):
            if seq != '':
                fasta_seq[marker] = seq
            marker = line.strip().split()[0][1:]
            seq = ''
        else:
            seq += line.strip()
        line = inp.readline()
        
        # Last sequence
        if seq != '':
            fasta_seq[marker] = seq
    inp.close()
    print('  \n# Number of sequences in the fasta file: ', fasta, len(fasta_seq))
    return(fasta_seq)


def update_db_fasta_sequences(new_fasta_seq, db_fasta_seq, outf):
    outp = open(outf, 'w')
    cnt = 0
    for marker in db_fasta_seq:
        if marker in new_fasta_seq:
            outp.write('>' + marker + '\n')
            outp.write(new_fasta_seq[marker] + '\n')
            cnt += 1
        else:
            outp.write('>' + marker + '\n')
            outp.write(db_fasta_seq[marker] + '\n')
    print('  # Number of sequences updated in db_fasta:', cnt)
    outp.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('new_fasta',
                        help='A file containing updated fasta sequences')

    parser.add_argument('db_fasta',
                        help='Database fasta file that need to be updated')

    args=parser.parse_args()

    new_fasta_seq = collect_fasta_sequences(args.new_fasta)
    
    db_fasta_seq = collect_fasta_sequences(args.db_fasta)

    outf = args.db_fasta.replace('.fa', '_109bpRefAlt.fa')
    update_db_fasta_sequences(new_fasta_seq, db_fasta_seq, outf)
