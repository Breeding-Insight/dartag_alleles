#!/usr/bin/python3

def get_ref_alt_bases(markerIDs):
    '''
    Added this function to consider the case of Indels
    Get reference and alternate bases for each markerID
    # OFP20_M6_CDS_75	M6_chr10_48867893_000000225	M6_chr10_48867893	225	    -	   AGC	Indel
    # OFP20_M6_CDS_290	M6_chr10_48867893_000000441	M6_chr10_48867893	441	TCACGATGT	-	Indel
    [       0                               1             2               3    4        5      6]
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
                # OFP20_M6_CDS_75	M6_chr10_48867893_000000225	M6_chr10_48867893	225	    -	   AGC	Indel
                outp.write('>' + key + '|Ref\n' + value + '\n')
                outp.write('>' + key + '|Alt\n' + value[:int(flank_len)] + alt + value[int(flank_len):] + '\n')
            elif alt == '-':
                # OFP20_M6_CDS_290	M6_chr10_48867893_000000441	M6_chr10_48867893	441	TCACGATGT	-	Indel
                del_in_alt = len(ref)
                outp.write('>' + key + '|Ref\n' + value + '\n')
                outp.write('>' + key + '|Alt\n' + value[:int(flank_len)] + value[int(flank_len) + del_in_alt:] + '\n')
            else:
                if ref_in_flank == ref:
                    #print('match', ref_in_flank, ref)
                    outp.write('>' + key + '|Ref\n' + value + '\n')
                    outp.write('>' + key + '|Alt\n' + value[:int(flank_len)] + alt + value[int(flank_len) + 1:] + '\n')
                elif ref_in_flank == alt:
                    # There is a special case in the potato panel where the ref and alt bases are swapped
                    # For the microhap db, will use the probe design file as the reference
                    outp.write('>' + key + '|Ref\n' + value[:int(flank_len)] + ref + value[int(flank_len) + 1:] + '\n')
                    outp.write('>' + key + '|Alt\n' + value + '\n')
                    print('\n# Ref and alt swapped:', key)
                    print('  # Ref in flank sequence:', ref_in_flank)
                    print('  # Ref in probe design file:', ref, ref_alt_bases[key])
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
