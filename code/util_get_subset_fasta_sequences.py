#!/usr/bin/python3

def collect_seqIDs(targets):
    inp = open(targets)
    line = inp.readline()
    seqIDs = []
    while line:
        line_array = line.strip().split()
        if line_array[0] not in seqIDs:
            seqIDs.append(line_array[0])
        else:
            print('This sequence ID is duplicated.')
        line = inp.readline()
    inp.close()
    print('Number of sequences: ', len(seqIDs))
    return(seqIDs)



def get_fasta(input, seqIDs, outf):
    inp = open(input)
    outp = open(outf, 'w')
    line = inp.readline()
    match = 'false'
    while line:
        if line.startswith('>'):
            line_array = line.strip().split()
            if line_array[0][1:] in seqIDs:
                outp.write(line)
                match = 'true'
            else:
                match = 'false'
        else:
            if match == 'true':
                outp.write(line)
            else:
                pass
        line = inp.readline()
    inp.close()
    outp.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract a subset of fasta sequences")

    parser.add_argument('targets',
                        help='A file containing target seqIDs')

    parser.add_argument('input',
                        help='Fasta sequence file')

    parser.add_argument('outf', help='Output file name')

    args=parser.parse_args()

    seqIDs = collect_seqIDs(args.targets)
    
    get_fasta(args.input, seqIDs, args.outf)
