#!/usr/bin/python3

def collect_ref_alt_sequence(hapDB):
    inp = open(hapDB)
    line = inp.readline()
    seq_dict = {}
    seq = ''
    while line:
        line = line.strip()
        if line.startswith('>'):
            if seq != '':
                seq_dict[id] = seq
            seq = ''
            id = line[1:]
        else:
            seq += line
        line = inp.readline()
    # Add last sequence to dict
    seq_dict[id] = seq
    inp.close()
    return(seq_dict)


def add_missing_ref_alt_to_madc(seq_dict, report):
    inp = open(report)
    outp = open(report.replace('.csv', '_addZero.csv'), 'w')
    line = inp.readline()
    alleles = {}
    alleles_full = {}
    while line:
        line_array = line.strip().split(',')
        if line_array[0] == '' or line_array[0] == '*' or line_array[0] == 'AlleleID':
            outp.write(line)
        else:
            cols = len(line_array)
            # Put all haplotypes into dict
            # Use cloneID as key, and a list of lists for all haplotypes
            if line_array[1] in alleles:
                alleles[line_array[1]].append(line_array[0])
                alleles_full[line_array[1]].append(line)
            else:
                alleles[line_array[1]] = [line_array[0]]
                alleles_full[line_array[1]] = [line]
        line = inp.readline()
    inp.close()

    # Process the dict
    for key, value in alleles.items():
        if key + '|Ref' not in value:
            print('Marker loci without Ref haplotype: ', key)
            alleleID = key + '|Ref'
            cloneID = key
            seq = seq_dict[alleleID + '_0001']
            ref_record = ','.join([alleleID, cloneID, seq, 'NA'] + ['0'] * (cols - 4)) + '\n'
            alleles_full[cloneID].insert(0, ref_record)
        else:
            pass
        
        if key + '|Alt' not in value:
            print('Marker loci without Alt haplotype: ', key)
            alleleID = key + '|Alt'
            cloneID = key
            seq = seq_dict[alleleID + '_0002']
            alt_record = ','.join([alleleID, cloneID, seq, 'NA'] + ['0'] * (cols - 4)) + '\n'
            alleles_full[cloneID].insert(1, alt_record)
        else:
            pass

    for key, value in alleles_full.items():
        outp.write(''.join(value))
    outp.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Update allele sequences in report with alleles assigned temporary names")

    parser.add_argument('hapDB',
                        help='Haplotype db v0001 with only Ref and Alt haplotypes')

    parser.add_argument('report',
                        help='DArTag report')

    args=parser.parse_args()

    seq_dict = collect_ref_alt_sequence(args.hapDB)

    add_missing_ref_alt_to_madc(seq_dict, args.report)
