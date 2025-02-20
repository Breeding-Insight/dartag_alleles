#!/usr/bin/python3

def get_ref_alt_bases(markerIDs):
    '''
    Get reference and alternate bases for each markerID
    snpST00073,chr11_002499502,chr11,2499502,T,C
    [   0         1             2       3    4 5]
    '''
    inp = open(markerIDs, 'r')
    line = inp.readline()
    ref_alt_bases={}
    while line:
        line_array = line.strip().split(',')
        ref_alt_bases[line_array[1]] = [line_array[4], line_array[5]]
        line = inp.readline()
    inp.close()
    return ref_alt_bases


def pre_ref_alt_flankSeq(flankSeq, ref_alt_bases):
    inp = open(flankSeq, 'r')
    line = inp.readline()
    flank = {}
    seq = ''
    while line:
        if line.startswith('>'):
            if seq != '':
                flank[markerID] = seq
            else:
                pass
            markerID = line.strip()[1:]
            seq = ''
        else:
            seq += line.strip()
        line = inp.readline()
    # Last sequence
    flank[markerID] = seq
    outp = open(flankSeq.replace('.fa', '_ref_alt.fa'), 'w')
    for key, value in flank.items():
        if key in ref_alt_bases:
            ref = value[0]
            alt = value[1]
            if value[180] == ref:
                outp.write('>' + key + '|Ref\n' + value + '\n')
                outp.write('>' + key + '|Alt\n' + value[:180] + alt + value[181:] + '\n')
            else:
                print(value[180], ref)
    inp.close()
    outp.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract f180 bp sequences (Ref and Alt) of 3K DArTag panel and generate snpID lookup table")

    parser.add_argument('markerIDs',
                        help='Marker ID look up table with ref and alt bases')
    
    parser.add_argument('flankSeq',
                        help='Flanking sequences of the markers in fasta format')

    args=parser.parse_args()

    ref_alt_bases = get_ref_alt_bases(args.markerIDs)

    pre_ref_alt_flankSeq(args.flankSeq, ref_alt_bases)
