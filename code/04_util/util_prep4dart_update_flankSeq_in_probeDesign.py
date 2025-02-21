#!/usr/bin/python3

def collect_new_flankSeq(flankSeq):
    inp = open(flankSeq)
    line = inp.readline()
    flankSeq_dict = {}
    seq = ''
    while line:
        if line.startswith('>'):
            if seq != '':
                flankSeq_dict[seqID] = seq
            else:
                pass
            seqID = line.strip()[1:]
            seq = ''
        else:
            seq = seq + line.strip()
        line = inp.readline()
    # Last sequence
    flankSeq_dict[seqID] = seq
    inp.close()
    print('Number of sequences added to the dictionary: ', len(flankSeq_dict))
    return(flankSeq_dict)


def update_flankSeq_in_probe(probe, flankSeq_dict, flankSeq_len):
    inp = open(probe)
    outp = open(probe.replace('.csv', '_f300bp.csv'), 'w')
    line = inp.readline() # header
    outp.write(line)
    line = inp.readline()
    while line:
        line_array = line.strip().split(',')
        if line_array[0] in flankSeq_dict:
            target_seq = flankSeq_dict[line_array[0]][0:int(flankSeq_len)] + line_array[5] + flankSeq_dict[line_array[0]][int(flankSeq_len) + 1:]
            #print(line_array[1])
            #print(target_seq, '\n')
            outp.write(','.join([line_array[0], target_seq] + line_array[2:]) + '\n')
        else:
            outp.write(line)
        line = inp.readline()
    inp.close()
    outp.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('flankSeq',
                        help='Flanking sequences of markers to be used to update probe design')

    parser.add_argument('probe',
                        help='Probe design file provided to DArT for QC')
    
    parser.add_argument('flankSeq_len',
                        help='Length of flanking sequence to the left of the markers')

    args=parser.parse_args()

    flankSeq_dict = collect_new_flankSeq(args.flankSeq)

    update_flankSeq_in_probe(args.probe, flankSeq_dict, args.flankSeq_len)
