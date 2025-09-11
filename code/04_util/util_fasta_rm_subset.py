#!/usr/bin/python3


def get_subset_seqIDs(subset):
    subsetIDs = []
    inp = open(subset)
    line = inp.readline()
    while line:
        line = line.strip()
        subsetIDs.append(line)
        line = inp.readline()
    inp.close()
    print(f'# Number of sequences needing removal: {len(subsetIDs)}')
    return subsetIDs


    
def rm_subset_seq_from_fasta(fasta, subsetIDs, outf):
    inp = open(fasta)
    outp = open(outf, 'w')
    line = inp.readline()
    remove = 'false'
    cnt = 0
    while line:
        if line.startswith('>'):
            line_array = line.strip().split()
            if line_array[0][1:] in subsetIDs:
                remove = 'true'
                cnt += 1
            else:
                outp.write(line)
                remove = 'false'
        else:
            if remove == 'false':
                outp.write(line)
            else:
                pass
        line = inp.readline()
    print(f'# Number of sequences removed from fasta file: {cnt}')
    inp.close()
    outp.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract a subset of fasta sequences")

    parser.add_argument('subset',
                        help='File with sequence IDs to extract')
    
    parser.add_argument('fasta',
                        help='Fasta sequence file')
    
    parser.add_argument('outf', help='Output file name')
    
    args=parser.parse_args()

    subsetIDs = get_subset_seqIDs(args.subset)

    rm_subset_seq_from_fasta(args.fasta, subsetIDs, args.outf)
