#!/usr/bin/python3

def collect_81bp_microhaplotypes_from_microhapDB(microhap):
    inp = open(microhap)
    line = inp.readline()
    microhap_81bp = {}
    seq = ''
    while line:
        if line.startswith('>'):
            if seq != '':
                microhap_81bp[seqID] = seq
            else:
                pass
            seqID = line.strip()[1:]
            seq = ''
        else:
            seq += line.strip()
        line = inp.readline()

    # Last sequence
    microhap_81bp[seqID] = seq
    inp.close()
    print('Number of microhaplotypes in the DB: ', len(microhap_81bp))
    #print(microhap_81bp)
    return(microhap_81bp)


def update_madc_microhap_to_81bp(report, microhap_81bp):
    inp = open(report)
    outp = open(report.replace('.csv', '_81bp.csv'), 'w')
    header = inp.readline()
    # AlleleID,CloneID,AlleleSequence,sample1,
    outp.write(header)
    line = inp.readline()
    while line:
        line_array = line.strip().split(',')
        if line_array[0] in microhap_81bp:
            outp.write(','.join(line_array[:2] + [microhap_81bp[line_array[0]]] + line_array[3:]) + '\n')
        else:
            outp.write(line)
            print('This microhaplotype is not present in the DB:', line_array[0])
        line = inp.readline()
    inp.close()
    outp.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('microhap',
                        help='Microhaplotype db extended to 81 bp')

    parser.add_argument('report',
                        help='MADC report with alleles assigned fixed IDs')

    args=parser.parse_args()

    microhap_81bp = collect_81bp_microhaplotypes_from_microhapDB(args.microhap)

    update_madc_microhap_to_81bp(args.report, microhap_81bp)
