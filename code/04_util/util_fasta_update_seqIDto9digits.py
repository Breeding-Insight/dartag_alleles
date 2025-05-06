#!/usr/bin/python3

def update_fasta_seqID(fasta):
    inp = open(fasta)
    line = inp.readline()
    seq = ''
    cnt = 0
    outp = open(fasta.replace('.fa', '_seqID9digits.fa'), 'w')
    while line:
        if line.startswith('>'):
            if seq != '':
                outp.write('>' + marker_new + '\n' + seq + '\n')
                cnt += 1
            else:
                pass
            marker = line.strip().split('|')[0][1:]
            allele = line.strip().split('|')[1]
            chr = marker.rsplit('_', 1)[0]
            pos = marker.rsplit('_', 1)[1]
            marker_new = chr + '_' + pos.zfill(9) + '|' + allele
            seq = ''
        else:
            seq += line.strip()
        line = inp.readline()
        
    # Last sequence
    if seq != '':
        cnt += 1
        outp.write('>' + marker_new + '\n' + seq + '\n')
        
    inp.close()
    print('# Number of sequences in the fasta file: ', cnt)
    outp.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('db_fasta',
                        help='Database fasta file that need to be updated')

    args=parser.parse_args()

    update_fasta_seqID(args.db_fasta)
