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

def get_flankSeq(flankSeq):
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
    inp.close()
    return(flank)


def pre_ref_alt_flankSeq(flank, ref_alt_bases, flank_len, outf):
    outp = open(outf, 'w')
    for key, value in flank.items():
        if key in ref_alt_bases:
            ref = ref_alt_bases[key][0]
            alt = ref_alt_bases[key][1]
            ref_in_flank = value[int(flank_len)]
            if ref == '-':
                # print('match', ref_in_flank, ref)
                outp.write('>' + key + '|Ref\n' + value + '\n')
                outp.write('>' + key + '|Alt\n' + value[:int(flank_len)] + alt + value[int(flank_len) + 1:] + '\n')
            elif alt == '-':
                indel_len = len(ref)
                outp.write('>' + key + '|Ref\n' + value + '\n')
                outp.write('>' + key + '|Alt\n' + value[:int(flank_len)] + value[int(flank_len) + indel_len:] + '\n')
            else:
                if ref_in_flank == ref:
                    #print('match', ref_in_flank, ref)
                    outp.write('>' + key + '|Ref\n' + value + '\n')
                    outp.write('>' + key + '|Alt\n' + value[:int(flank_len)] + alt + value[int(flank_len) + 1:] + '\n')
                elif ref_in_flank == alt:
                    # There is a special case in the potato panel where the ref and alt bases are swapped
                    # For the microhap db, will use the probe design file as the reference
                    outp.write('>' + key + '|Ref\n' + value[:int(flank_len)] + alt + value[int(flank_len) + 1:] + '\n')
                    outp.write('>' + key + '|Alt\n' + value + '\n')
                    print('\nRef and alt swapped:', key, ref_in_flank, ref)
                    print(key, ref_in_flank, ref_alt_bases[key])
                else:
                    print('\nmiss', key, ref_in_flank, ref)
                    print(key, ref_in_flank, ref_alt_bases[key])
    outp.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract f180 bp sequences (Ref and Alt) of 3K DArTag panel and generate snpID lookup table")

    parser.add_argument('markerIDs',
                        help='Marker ID look up table with ref and alt bases')
    
    parser.add_argument('flankSeq',
                        help='Flanking sequences of the markers in fasta format')
    
    parser.add_argument('flank_len',
                        help='Length of the left-side flanking sequence')

    args=parser.parse_args()

    ref_alt_bases = get_ref_alt_bases(args.markerIDs)

    flank = get_flankSeq(args.flankSeq)

    outf = args.flankSeq.replace('.fa', '_ref_alt.fa')
    pre_ref_alt_flankSeq(flank, ref_alt_bases, args.flank_len, outf)
