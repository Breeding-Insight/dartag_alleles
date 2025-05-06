#!/usr/bin/python3

def update_seqID(lut):
    inp = open(lut)
    line = inp.readline()
    cnt = 0
    outp = open(lut.replace('.txt', '_seqID9digits.txt'), 'w')
    while line:
        # Chr01_00085461|RefMatch	0	0
        line_array = line.strip().split('\t')
        marker = line_array[0].split('|')[0]
        allele = line_array[0].split('|')[1]
        chr = marker.rsplit('_', 1)[0]
        pos = marker.rsplit('_', 1)[1]
        marker_new = chr + '_' + pos.zfill(9) + '|' + allele
        outp.write('\t'.join([marker_new] + line_array[1:]) + '\n')
        cnt += 1
        line = inp.readline()
    inp.close()
    print('# Number of entries written out: ', cnt)
    outp.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('lut',
                        help='Match count LUT file to be updated with 9 digit seqID')

    args=parser.parse_args()

    update_seqID(args.lut)
