#!/usr/bin/python3


def check_iupac_code(seq):
    iupac_code = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
    for base in seq:
        if base in iupac_code:
            return True
    return False
    
    
def check_fasta_sequences(fasta):
    inp = open(fasta)
    line = inp.readline()
    seq = ''
    cnt = 0
    while line:
        if line.startswith('>'):
            if seq != '':
                iupac = check_iupac_code(seq)
                if iupac == True:
                    cnt += 1
            marker = line.strip().split()[0][1:]
            seq = ''
        else:
            seq += line.strip()
        line = inp.readline()
        
    # Last sequence
    if seq != '':
        iupac = check_iupac_code(seq)
        if iupac == True:
            cnt += 1
    inp.close()
    print('  \n# Number of sequences with IUPAC bases: ', cnt)


def check_madc_sequences(madc):
    inp = open(madc)
    line = inp.readline()
    total = 0
    cnt = 0
    markers = {}
    while line:
        if line.startswith('AlleleID') or line.startswith('*') or line.startswith(' '):
            pass
        else:
            total += 1
            line_array = line.strip().split(',')
            iupac = check_iupac_code(line_array[2])
            if iupac == True:
                cnt += 1
                if line_array[1] not in markers:
                    markers[line_array[1]] = 1
                else:
                    markers[line_array[1]] += 1
        line = inp.readline()
    inp.close()
    print('  \n# Total number of sequences: ', total)
    print('  # Number of marker loci with IUPAC bases: ', len(markers))
    print('  # Number of microhaplotypes with IUPAC bases: ', cnt)
    #print(markers)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('input_file',
                        help='A multi-fasta file or a MADC file')
    
    parser.add_argument('fasta_or_madc',
                        help='Type of input file: fasta or madc')
    
    args=parser.parse_args()

    if args.fasta_or_madc == 'madc':
        check_madc_sequences(args.input_file)
    else:
        check_fasta_sequences(args.input_file)
