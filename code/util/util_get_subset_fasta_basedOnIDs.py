#!/usr/bin/python3


def get_fasta(input, seqID_pattern, outf):
    inp = open(input)
    outp = open(outf, 'w')
    line = inp.readline()
    match = 'false'
    while line:
        if line.startswith('>'):
            line_array = line.strip().split()
            if seqID_pattern in line_array[0][1:]:
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

    parser.add_argument('input',
                        help='Fasta sequence file')
    
    parser.add_argument('outf', help='Output file name')

    parser.add_argument('pattern',
                        help='Pattern to match in sequence IDs')
    
    args=parser.parse_args()
    
    get_fasta(args.input, args.pattern, args.outf)
